% VIZ_PREDICTOR_ERROR  Forecast skill vs horizon.
%
% On a held-out operating day (full diurnal demand at nominal flow, the regime the
% controller runs in, at a phase the fit never saw) I forecast every horizon from 1 to 16 steps (15 min .. 4 h)
% two ways and score each against the plant: the direct multi-step Koopman map
% (V_seq) and the iterated one-step Koopman map (A,B). Skill is the R^2 on delivered
% heat, averaged over the five consumers. The direct map keeps R^2 near 1 across the
% whole horizon; the iterated map starts level at one step and then falls away as its
% rollout compounds the per-step error. That growing gap is the open-loop result
% behind the direct multi-step choice, and the same ranking on the held-out PRBS test
% is checked as an invariant in validate_predictor.
%
% Runs live: loads only the trained models, generates the day here, scores here.

clear; clc;

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(root); startup();

p = params();
pred_dir = fullfile(root, 'predictor_open_loop', 'results');

% the two trained predictors, both from predictor_open_loop/results:
%   fits     = direct multi-step V_seq, one (cp, dp) per (horizon, consumer)
%              from fit_vseq; cp hits the lift z, dp hits the input/demand block
%   fit_iter = the one-step Koopman map A,B from fit_iterated, rolled out below
F = load(fullfile(pred_dir, 'vseq_fits_full.mat'));  fits     = F.fits;
G = load(fullfile(pred_dir, 'iterated_AB.mat'));      fit_iter = G.fit_iter;
H      = fits.horizons(:)';   % the horizons the fits cover: 1..16 (= 1:p.o1.H_max, 15 min .. 4 h)
nH     = numel(H);
n_user = size(fits.cp, 2);    % consumers = number of columns in the per-(h,i) fit grid (5)
% true when the V_seq fit folds the demand forecast into its input block (it does, see fit_vseq);
% tells head_apply below to split dp into a control part and a demand part
has_dV = isfield(fits, 'includes_d_in_V') && fits.includes_d_in_V;

% --- generate a normal operating day: nominal flow, full diurnal demand ---
% This is the deployment regime (smooth diurnal day at design flow) rather than the
% PRBS the fits were mostly trained on. The operational days in the training set sit
% on a 3 h phase grid (t_offset = 0, 3, ... 21 h, see data/generate_operational.m), so
% I score on a day at 10.5 h, half way between two training phases. That day is not in
% the fit, which is what makes this a held-out horizon test rather than a replay of a
% training day.
t_off_h = 10.5;                % phase of the day scored here [h], off the training grid
[net, z0] = build_plant(p);
mdotE  = net.mdotEdges(:);     % the network's nominal per-edge mass flows [kg/s] = the design operating point
F0_idx = find(strcmp({net.Nodes.name}, 'F0'));   % source feed node, for T_F0 below
R0_idx = find(strcmp({net.Nodes.name}, 'R0'));   % source return node, for T_0r below
% r_q_fun pins flow at the nominal edge setpoints (drives q as a state instead of the
% scalar exogenous w); t_offset shifts where in the diurnal cycle the day starts
p_run = p; p_run.t_offset = t_off_h * 3600; p_run.r_q_fun = @(t) mdotE;
% u_fun = @(t) 0  -> zero source-supply offset, so supply sits at p.Tin_nom all day
% w_fun = @(t) 1  -> exogenous flow scale, ignored here since r_q_fun drives the flow
% 26 h so there is a full day of valid samples after the 2 h warm-up is dropped
res = simulate_plant(net, z0, p_run, @(t) 0, @(t) 1.0, 26 * 3600);

traj.t = res.t; traj.t_offset = p_run.t_offset;
traj.T_0s = p.Tin_nom + res.u;  traj.r_q = res.r_q;  traj.T_ir = res.T_r_i;
traj.d = res.d_i;  traj.theta = (res.T_s_i - res.T_r_i) .* res.q_users;
traj.T_is = res.T_s_i;  traj.T_0r = res.Tout(R0_idx, :);  traj.T_F0 = res.Tout(F0_idx, :);
traj.q_users = res.q_users;  traj.q_edges = res.q_edges;  traj.Tout = res.Tout;
% Z = the Koopman lift z_k of this day (n_z x N); 47 features = 36 base + 11 exergy
% bilinears in production. meta.idx says where each feature sits in z.
[Z, ~, meta] = candidate_library(traj, p);
% U = the input vector the predictors expect, stacked the same way as in the fit:
%   T_0s (1 source supply temp), r_q (29 edge flow setpoints), T_ir (5 consumer
%   return-temp setpoints) -> n_u = 35 rows. D = per-consumer heat demand [W] (the forecast).
U = [traj.T_0s(:)'; traj.r_q; traj.T_ir];  D = traj.d;
% th_ix = rows of z holding theta^i = (Ts-Tr)*q per consumer; n_u read off U for head_apply's dp split
th_ix = meta.idx.theta;  n_u = size(U, 1);

% --- R^2 of each predictor at each horizon, per consumer ---
r2 = @(pred, act) 1 - sum((pred - act).^2) / max(sum((act - mean(act)).^2), eps);
R2_v = zeros(nH, n_user);  R2_i = zeros(nH, n_user);
for hi = 1:nH
    h  = H(hi);
    ks = max(1, meta.valid_start);
    ke = size(Z, 2) - h;
    pv = zeros(n_user, ke-ks+1);  pii = zeros(n_user, ke-ks+1);  act = zeros(n_user, ke-ks+1);
    col = 0;
    for k = ks:ke
        col = col + 1;
        % iterated leg: start at the true lift z_k and step the one-step map h times.
        % each step's input carries the applied controls U plus the next-step demand
        % D(k+s+1), the same info the V_seq fit gets, so this compares structure not inputs.
        z_it = Z(:, k);
        for s = 0:h-1, z_it = fit_iter.A * z_it + fit_iter.B * [U(:, k+s); D(:, k+s+1)]; end
        % direct leg: the V_seq input window is the h applied controls and the h-step
        % demand forecast, exactly the regressor block fit_vseq trained on
        Vu = U(:, k : k+h-1);  Vd = D(:, k+1 : k+h);
        for i = 1:n_user
            % all three are delivered heat c = p.cp * theta [W]; p.cp converts theta [K.kg/s] to W
            pv(i, col)  = p.cp * head_apply(fits.cp{hi, i}, fits.dp{hi, i}, Z(:, k), Vu, Vd, n_u, h, has_dV);
            pii(i, col) = p.cp * z_it(th_ix(i));   % read theta out of the rolled lift
            act(i, col) = p.cp * traj.theta(i, k + h);   % plant truth at k+h
        end
    end
    for i = 1:n_user
        R2_v(hi, i) = r2(pv(i, :),  act(i, :));
        R2_i(hi, i) = r2(pii(i, :), act(i, :));
    end
end
r2_v = mean(R2_v, 2)';   % suite mean across consumers
r2_i = mean(R2_i, 2)';

% --- figure: skill vs horizon ---
figure('Color', 'w');
plot(H, r2_i, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
plot(H, r2_v, 'b-o', 'LineWidth', 1.8, 'MarkerSize', 5);
grid on; xlabel('forecast horizon [15-min steps]'); ylabel('R^2 (delivered heat)');
xlim([1 H(end)]); ylim([floor(min(r2_i)*20)/20, 1]);
legend('iterated one-step', 'direct multi-step', 'Location', 'southwest');
title('Forecast skill vs horizon, held-out operating day (all consumers)');

% WHAT YOU SEE: direct multi-step stays near R^2 = 1 the whole way out, while the
% iterated one-step is equal at h = 1 and then degrades as its rollout compounds.

exportgraphics(gcf, fullfile(here, 'viz_predictor_error.pdf'), 'ContentType', 'vector');
fprintf('Saved viz_predictor_error.pdf\n');
fprintf('  R^2 at %d steps:  iterated %.3f   direct multi-step %.3f\n', ...
        H(end), r2_i(end), r2_v(end));


function y = head_apply(cp_h, dp_h, z0, Vu, Vd, n_u, h, has_dV)
% Evaluate one V_seq head: theta_{k+h} = cp_h' z_k + dp_h' [controls; demand].
% dp_h is stacked control-block-then-demand-block, so the first n_u*h entries
% multiply the control window and the rest multiply the demand forecast.
if has_dV
    n_v_u = n_u * h;            % length of the control block (35 inputs x h steps)
    y = cp_h(:)' * z0 + dp_h(1:n_v_u)' * Vu(:) + dp_h(n_v_u+1:end)' * Vd(:);
else
    y = cp_h(:)' * z0 + dp_h(:)' * Vu(:);   % legacy fit with no demand in V
end
end
