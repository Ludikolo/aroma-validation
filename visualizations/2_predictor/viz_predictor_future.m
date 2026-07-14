% VIZ_PREDICTOR_FUTURE  Predictor check 1 of 2: a 4 h-ahead forecast on a normal day.
%
% At every step I forecast the delivered heat 4 h ahead, two ways, and lay both on
% top of what the plant actually does: the direct multi-step Koopman map (one fit per
% horizon) and the iterated one-step map (A,B) rolled forward. The direct map stays on
% the truth; the iterated map falls behind as its one-step error compounds over the horizon.
% Shown for two high-demand consumers. (How the error grows with the horizon is
% quantified in viz_predictor_error.)
%
% Runs live: loads only the trained models, generates the day here, computes every
% forecast here; nothing is loaded from a results .mat.

clear; clc;

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(root); startup();

p = params(); p.Ts = 900;              % all plant/predictor constants; Ts = 900 s (15-min sample)
pred_dir = fullfile(root, 'predictor_open_loop', 'results');

% the two trained models I compare here, both saved by the open-loop fit step:
%   fits     = direct multi-step V_seq, one fit per (horizon, consumer); from fit_vseq.m.
%              fits.cp{h,i} weights the lift z, fits.dp{h,i} weights [controls; demand].
%   fit_iter = single one-step Koopman map z_{k+1} = A z_k + B [u_k; d_{k+1}], rolled forward.
F  = load(fullfile(pred_dir, 'vseq_fits_full.mat'));  fits     = F.fits;
G  = load(fullfile(pred_dir, 'iterated_AB.mat'));      fit_iter = G.fit_iter;
H      = fits.horizons(:)';
has_dV = isfield(fits, 'includes_d_in_V') && fits.includes_d_in_V;  % did the fit put demand in V_seq? (yes, matches fit_vseq)
h0  = 16;  hi0 = find(H == h0);         % 4 h ahead = 16 steps at Ts = 900 s; hi0 = its row in the fits

% --- generate a normal operating day: nominal flow, full diurnal demand ---
% this is the deployment regime (smooth diurnal demand at the design flow),
% NOT the PRBS training set the fits were learned on, so it's an honest check.
[net, z0_cold] = build_plant(p);       % network + cold start state
ei = edge_user_index(net);  mdotE = net.mdotEdges(:);   % mdotE = nominal edge flows, held fixed below
F0_idx = find(strcmp({net.Nodes.name}, 'F0'));
R0_idx = find(strcmp({net.Nodes.name}, 'R0'));
% hold flow at nominal (r_q_fun = mdotE), no source-temp offset (u = 0), full demand (scale 1.0), 26 h run
p_run = p; p_run.t_offset = 0; p_run.r_q_fun = @(t) mdotE;
res = simulate_plant(net, z0_cold, p_run, @(t) 0, @(t) 1.0, 26 * 3600);

traj.t = res.t; traj.t_offset = 0;
traj.T_0s = p.Tin_nom + res.u;  traj.r_q = res.r_q;  traj.T_ir = res.T_r_i;
% theta^(i) = (T_s - T_r)*q per consumer [K.kg/s]; multiply by cp later to get heat
traj.d = res.d_i;  traj.theta = (res.T_s_i - res.T_r_i) .* res.q_users;
traj.T_is = res.T_s_i;  traj.T_0r = res.Tout(R0_idx, :);  traj.T_F0 = res.Tout(F0_idx, :);
traj.q_users = res.q_users;  traj.q_edges = res.q_edges;  traj.Tout = res.Tout;
[Z, ~, meta] = candidate_library(traj, p);   % Z = 47-feature Koopman lift (36 base + 11 exergy); meta.idx.theta locates the 5 delivered-heat features

% U stacks the inputs the same way the fits expect: 35 rows = 1 source-supply-temp
% offset T_0s, then 29 edge flow setpoints r_q, then 5 consumer return-temp setpoints T_ir.
% D = per-consumer heat demand (the known forecast that lives in V_seq).
U = [traj.T_0s(:)'; traj.r_q; traj.T_ir];  D = traj.d;
n_u = size(U, 1);  th_ix = meta.idx.theta;    % n_u = 35; th_ix(i) = row of theta^(i) inside the lift
% window: skip the lift warm-up (valid_start), stop h0 steps early so target k+h0 stays in range
N = size(Z, 2);  ks = max(1, meta.valid_start);  ke = N - h0;  kk = ks:ke;
t_h = (kk + h0 - 1) * p.Ts / 3600;            % x-axis in hours (sample index -> h)
r2 = @(pred, act) 1 - sum((pred - act).^2) / max(sum((act - mean(act)).^2), eps);

cons = [2, 4];                          % two high-demand consumers
figure('Color', 'w');
for pi = 1:numel(cons)
    i0 = cons(pi);
    c_act = zeros(size(kk)); c_multi = zeros(size(kk)); c_one = zeros(size(kk));
    for n = 1:numel(kk)
        % z0 = lift now; Vu = the h0 applied controls; Vd = the demand forecast over the same window
        k = kk(n); z0 = Z(:, k);  Vu = U(:, k : k+h0-1);  Vd = D(:, k+1 : k+h0);
        % direct V_seq head at h0, then cp converts theta [K.kg/s] -> heat [W]
        c_multi(n) = p.cp * head_apply(fits.cp{hi0, i0}, fits.dp{hi0, i0}, z0, Vu, Vd, n_u, h0, has_dV);
        % iterated one-step rolled h0 steps: each step compounds its error
        z_it = z0;
        for s = 0:h0-1
            z_it = fit_iter.A * z_it + fit_iter.B * [U(:, k+s); D(:, k+1+s)];
        end
        c_one(n) = p.cp * z_it(th_ix(i0));
        c_act(n) = p.cp * traj.theta(i0, k + h0);   % ground truth heat the plant actually delivered
    end
    subplot(2, 1, pi);
    plot(t_h, c_act/1e3, 'k', 'LineWidth', 1.8); hold on;
    plot(t_h, c_multi/1e3, 'b', 'LineWidth', 1.3);
    plot(t_h, c_one/1e3, 'r', 'LineWidth', 1.1);
    grid on; xlabel('time [h]'); ylabel('delivered heat [kW]');
    legend('actual plant', 'direct multi-step', 'iterated one-step', 'Location', 'best');
    title(sprintf('4 h-ahead forecast, normal operation, consumer C%d: delivered heat (R^2 %.2f vs %.2f)', ...
          i0, r2(c_multi, c_act), r2(c_one, c_act)));
    fprintf('  C%d delivered-heat R2: multi %.3f   iter %.3f\n', i0, r2(c_multi, c_act), r2(c_one, c_act));
end

exportgraphics(gcf, fullfile(here, 'viz_predictor_future.pdf'), 'ContentType', 'vector');
fprintf('Saved viz_predictor_future.pdf  (4 h-ahead, normal operation)\n');


% one direct V_seq forecast: theta_{k+h} = cp'*z0 + dp'*[controls; demand].
% dp is split: the first n_u*h entries weight the h applied control vectors,
% the rest weight the h demand-forecast vectors (only present when has_dV).
function y = head_apply(cp_h, dp_h, z0, Vu, Vd, n_u, h, has_dV)
if has_dV
    n_v_u = n_u * h;                                  % length of the control block in dp
    y = cp_h(:)' * z0 + dp_h(1:n_v_u)' * Vu(:) + dp_h(n_v_u+1:end)' * Vd(:);
else
    y = cp_h(:)' * z0 + dp_h(:)' * Vu(:);             % no demand block: controls only
end
end
