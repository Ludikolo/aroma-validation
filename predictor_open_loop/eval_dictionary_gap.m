% EVAL_DICTIONARY_GAP  Score the spatial lift against the 47 production lift and
%   write the gap the dictionary study reports.
%
% The question the study asks: does adding the in-network (supply-ring) block to
% the lift predict delivered heat better than the 47 features I ship? Both lifts
% are fitted by fit_vseq_spatial under the same protocol, so the comparison is
% like for like:
%   'none' = the production 47 lift, refit here  -> vseq_fits_spatial_none.mat
%   'full' = the 47 plus the supply-ring block   -> vseq_fits_spatial_full.mat
% I score both on three regimes and report gap = RMSE(spatial) - RMSE(47) in W,
% so a positive bar means the 47 predicts better.
%
%   PRBS        the two held-out PRBS days (val_01, val_02), the identification regime
%   op-day      a held-out operating day, phase 10.5 h, off the 0:3:21 h training grid
%   on-policy   the trajectory the closed-loop KPC controller actually visits, which
%               is the only regime the controller is ever asked to predict in
%
% The on-policy leg matters most and is the reason the 47 stays: a lift can win on
% the identification data and still lose where the controller drives the plant.
%
% Writes results/dict_gap.mat; plot_dictionary_gap.m draws it. Runs in about
% 5 min (the on-policy leg is a 24 h closed-loop run).

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root); startup();

p   = params();
net = build_plant(p);
res_dir = fullfile(here, 'results');

S_full = load(fullfile(res_dir, 'vseq_fits_spatial_full.mat'));  f_full = S_full.fits;
S_none = load(fullfile(res_dir, 'vseq_fits_spatial_none.mat'));  f_none = S_none.fits;
H = f_none.horizons(:)';
assert(isequal(H, f_full.horizons(:)'), 'the two fits cover different horizons');
fprintf('lifts: 47 (n_z = %d) vs spatial (n_z = %d), horizons %d..%d\n', ...
        size(f_none.cp{1,1},1), size(f_full.cp{1,1},1), H(1), H(end));

%% --- the three regimes ---
trajs = struct('name', {}, 'traj', {});

% PRBS: the two held-out validation days
for j = 1:2
    S = load(fullfile(res_dir, sprintf('val_%02d.mat', j)));
    trajs(end+1) = struct('name', sprintf('PRBS val_%02d', j), 'traj', S.traj); %#ok<SAGROW>
end

% operational: a held-out day, same construction as the training days but at a
% phase halfway between two of them, so the fit has not seen it
fprintf('simulating the held-out operating day...\n');
trajs(end+1) = struct('name', 'op-day (held out)', 'traj', make_opday(p, net, 10.5)); %#ok<SAGROW>

% on-policy: the stressed benchmark trajectory (mdot_scale = 0.35), run in both
% directions. A lift can look good on a trajectory its own controller produced,
% so I score both lifts on the closed loop driven by the 47 AND on the closed
% loop driven by the spatial lift. Neither gets to pick its own operating point.
fprintf('running the closed loop driven by the 47 lift...\n');
trajs(end+1) = struct('name', 'on-policy (47-driven)', ...
                      'traj', make_onpolicy(p, 'none')); %#ok<SAGROW>
fprintf('running the closed loop driven by the spatial lift...\n');
trajs(end+1) = struct('name', 'on-policy (spatial-driven)', ...
                      'traj', make_onpolicy(p, 'full')); %#ok<SAGROW>

%% --- score both lifts on every regime and horizon ---
nR = numel(trajs);  nH = numel(H);
rmse_full = nan(nR, nH);  rmse_none = nan(nR, nH);

for r = 1:nR
    tj = trajs(r).traj;
    [Zf, ~, mf] = spatial_library(tj, p, net, 'full');
    [Zn, ~, mn] = spatial_library(tj, p, net, 'none');
    U = [tj.T_0s(:)'; tj.r_q; tj.T_ir];  D = tj.d;
    for hi = 1:nH
        rmse_full(r, hi) = score(Zf, U, D, tj.theta, mf, f_full, hi, H(hi), p);
        rmse_none(r, hi) = score(Zn, U, D, tj.theta, mn, f_none, hi, H(hi), p);
    end
    fprintf('  %-24s scored\n', trajs(r).name);
end

gap = rmse_full - rmse_none;         % > 0  ->  the 47 predicts better

%% --- report ---
fprintf('\nRMSE(spatial) - RMSE(47), delivered heat [W].  positive = the 47 wins\n');
fprintf('%-26s', 'regime');  fprintf('%9s', 'h=1', 'h=4', 'h=8', 'h=12', 'h=16');  fprintf('\n');
show = [1 4 8 12 16];
for r = 1:nR
    fprintf('%-26s', trajs(r).name);
    for h = show
        hi = find(H == h);
        if isempty(hi), fprintf('%9s', '-'); else, fprintf('%+9.0f', gap(r, hi)); end
    end
    fprintf('\n');
end
fprintf('\nfor scale, the 47 lift''s own RMSE [W] on the same cells:\n');
for r = 1:nR
    fprintf('%-26s', trajs(r).name);
    for h = show
        hi = find(H == h);
        if isempty(hi), fprintf('%9s', '-'); else, fprintf('%9.0f', rmse_none(r, hi)); end
    end
    fprintf('\n');
end

names = {trajs.name};
save(fullfile(res_dir, 'dict_gap.mat'), 'gap', 'rmse_full', 'rmse_none', 'H', 'names');
fprintf('\nsaved results/dict_gap.mat\n');


%% ---------- helpers ----------
function rm = score(Z, U, D, theta, meta, fits, hi, h, p)
% RMSE of predicted delivered heat at horizon h, pooled over the five consumers.
n_z = size(Z, 1);
ks  = max(1, meta.valid_start);
ke  = size(Z, 2) - h;
err = [];
for i = 1:size(theta, 1)
    cp = fits.cp{hi, i};  dp = fits.dp{hi, i};
    pred = zeros(1, ke - ks + 1);  act = pred;  col = 0;
    for k = ks:ke
        col = col + 1;
        V  = U(:, k:k+h-1);  Vd = D(:, k+1:k+h);
        phi = [Z(:, k); V(:); Vd(:)];
        pred(col) = cp(:)' * phi(1:n_z) + dp(:)' * phi(n_z+1:end);
        act(col)  = theta(i, k + h);
    end
    err = [err, p.cp * (pred - act)]; %#ok<AGROW>
end
rm = sqrt(mean(err.^2));
end

function traj = make_opday(p, net, t_off_h)
% one operating day at a phase the fit never saw (training grid is 0:3:21 h)
[~, z0] = build_plant(p);
mdotE = net.mdotEdges(:);
R0 = find(strcmp({net.Nodes.name}, 'R0'));  F0 = find(strcmp({net.Nodes.name}, 'F0'));
pr = p; pr.t_offset = t_off_h * 3600; pr.r_q_fun = @(t) mdotE;
res = simulate_plant(net, z0, pr, @(t) 0, @(t) 1.0, 26 * 3600);
traj = struct();
traj.t = res.t;  traj.t_offset = pr.t_offset;
traj.T_0s = p.Tin_nom + res.u;  traj.r_q = res.r_q;  traj.T_ir = res.T_r_i;
traj.T_is = res.T_s_i;  traj.T_0r = res.Tout(R0, :);  traj.T_F0 = res.Tout(F0, :);
traj.q_users = res.q_users;  traj.q_edges = res.q_edges;  traj.Tout = res.Tout;
traj.d = res.d_i;  traj.c = res.c_i;
traj.theta = (res.T_s_i - res.T_r_i) .* res.q_users;
end

function traj = make_onpolicy(p, driver)
% One 24 h closed-loop KPC run on the stressed benchmark scenario
% (mdot_scale = 0.35, the same setup run_comparison uses), driven by the lift
% named in `driver`: 'none' = the production 47, 'full' = the spatial lift.
% kpc_step_loop returns the full trajectory as its second output.
root = fileparts(fileparts(mfilename('fullpath')));
pred_res = fullfile(root, 'predictor_open_loop', 'results');
ctrl_res = fullfile(root, 'controller_closed_loop', 'results');

mdot_scale = 0.35;  T_warm = 30 * 60;  T_sim = 24 * 3600;  sc_start = 0;
[net, z0_cold] = build_plant(p);
net.mdotEdges = mdot_scale * net.mdotEdges;
net.flows = net.mdotEdges;  net.q0 = net.mdotEdges(:);

% the controller's own predictor: the 47 in production, the spatial lift when
% we ask what trajectory that controller would have visited instead
if strcmp(driver, 'full')
    F = load(fullfile(pred_res, 'vseq_fits_spatial_full.mat'));
    lift_fn = @(tj, pp) spatial_library(tj, pp, net, 'full');
else
    F = load(fullfile(pred_res, 'vseq_fits_full.mat'));
    lift_fn = [];
end
Hb = load(fullfile(ctrl_res, 'best_horizons.mat'));
B  = load(fullfile(ctrl_res, 'best_tune.mat'));
A  = load(fullfile(ctrl_res, 'best_alpha.mat'));
pred = truncate_pred_to_Np(build_kpc_v2_matrices(F.fits), Hb.Np_best);

tune = B.best_tune;
tune.dT_0s_max = 7.5;  tune.dr_q_max = 4.5;  tune.dT_ir_max = 15.0;
tune.use_kirchhoff = true;  tune.use_Tr_le_Ts = true;  tune.use_mixing = true;
tune.rho_slack_mix = 1;  tune.Nc = Hb.Nc_best;  tune.alpha_energy = A.alpha_star;

ei = edge_user_index(net);
F0_idx = find(strcmp({net.Nodes.name}, 'F0'));
R0_idx = find(strcmp({net.Nodes.name}, 'R0'));
con_idx = zeros(ei.n_user, 1);
for i = 1:ei.n_user
    con_idx(i) = find(strcmp({net.Nodes.name}, ei.consumers{i}));
end
p_wu = p; p_wu.t_offset = sc_start; p_wu.r_q_fun = @(t) net.mdotEdges(:);
res_wu = simulate_plant(net, z0_cold, p_wu, @(t) 0, @(t) 1.0, T_warm);
if isempty(lift_fn)
    [~, traj] = kpc_step_loop(net, res_wu, p, sc_start, T_warm, T_sim, ...
                              pred, tune, ei, F0_idx, R0_idx, con_idx);
else
    [~, traj] = kpc_step_loop(net, res_wu, p, sc_start, T_warm, T_sim, ...
                              pred, tune, ei, F0_idx, R0_idx, con_idx, lift_fn);
end
end
