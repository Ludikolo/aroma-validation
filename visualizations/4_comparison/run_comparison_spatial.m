% NOTE: reproduces the spatial dictionary comparison in dictionary_study.pdf; not on the run_all_tests path.
% RUN_COMPARISON_SPATIAL  Closed-loop leg of the dictionary study: the two
% Koopman controllers (KPC and the iterated Koopman-LMPC) re-run with the
% spatial full lift, on exactly the scenario, tune and envelope of
% run_comparison. The NMPC and Jacobian-LMPC legs do not use the lift, so
% they are read from the saved comparison.mat as the reference rows instead
% of being re-run. Saves to comparison_spatial.mat; nothing existing is
% overwritten.

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(root); addpath(fullfile(root, 'controller_closed_loop', 'lmpc')); startup();

% same stressed scenario as run_comparison
mdot_scale = 0.35;
sc_start   = 0;
T_warm     = 30 * 60;
T_sim      = 24 * 3600;

p = params();
pred_dir = fullfile(root, 'predictor_open_loop', 'results');
ctrl_dir = fullfile(root, 'controller_closed_loop', 'results');

% spatial fits (with the extras from extend_vseq_extras_spatial) + locked tune
F  = load(fullfile(pred_dir, 'vseq_fits_spatial_full.mat'));
assert(strcmp(F.fits.dict, 'spatial_full'), 'expected spatial full fits');
assert(isfield(F.fits, 'Ts_cp'), 'run extend_vseq_extras_spatial first (extras missing)');
Hb = load(fullfile(ctrl_dir, 'best_horizons.mat'));
Ab = load(fullfile(ctrl_dir, 'best_alpha.mat'));
Bb = load(fullfile(ctrl_dir, 'best_tune.mat'));
pred = truncate_pred_to_Np(build_kpc_v2_matrices(F.fits), Hb.Np_best);
tune = Bb.best_tune;
tune.dT_0s_max = 7.5; tune.dr_q_max = 4.5; tune.dT_ir_max = 15.0;
tune.use_kirchhoff = true; tune.use_Tr_le_Ts = true; tune.use_mixing = true; tune.rho_slack_mix = 1;
tune.Nc = Hb.Nc_best; tune.alpha_energy = Ab.alpha_star;
Np = Hb.Np_best;

S = load(fullfile(pred_dir, 'iterated_AB_spatial.mat'));
fit_iter = S.fit_iter;
assert(strcmp(fit_iter.dict, 'spatial_full'), 'expected the spatial iterated fit');
tune_lmpc = tune; tune_lmpc.use_mixing = false;

% plant at stressed capacity + shared warmup, as in run_comparison
[net, z0_cold] = build_plant(p);
net.mdotEdges = mdot_scale * net.mdotEdges; net.flows = net.mdotEdges; net.q0 = net.mdotEdges(:);
ei = edge_user_index(net); n_user = ei.n_user;
F0_idx = find(strcmp({net.Nodes.name}, 'F0'));
R0_idx = find(strcmp({net.Nodes.name}, 'R0'));
con_idx = arrayfun(@(i) find(strcmp({net.Nodes.name}, ei.consumers{i})), 1:n_user);

p_wu = p; p_wu.t_offset = sc_start; p_wu.r_q_fun = @(t) net.mdotEdges(:);
res_wu = simulate_plant(net, z0_cold, p_wu, @(t) 0, @(t) 1.0, T_warm);

% the spatial lift, injected through the existing lift_fn hook; the fit-time
% and deploy-time lift must be the same variant, checked by a pre-flight lift
lift_fn = @(traj, pp) spatial_library(traj, pp, net, 'full');
Tprobe = load(fullfile(pred_dir, 'test_01.mat'));
[Zp, ~, mp] = lift_fn(Tprobe.traj, p);
assert(size(Zp,1) == pred.n_z && strcmp(mp.dict, F.fits.dict), ...
    'deploy lift does not match the fitted lift');

fprintf('Running KPC (spatial, n_z=%d) ...\n', pred.n_z); t0 = tic;
res_kpc = kpc_step_loop(net, res_wu, p, sc_start, T_warm, T_sim, pred, tune, ei, F0_idx, R0_idx, con_idx, lift_fn);
fprintf('  done (%.1f s wall)\n', toc(t0));

fprintf('Running Koopman-LMPC (spatial) ...\n'); t0 = tic;
res_klmpc = lmpc_step_loop(net, res_wu, p, sc_start, T_warm, T_sim, fit_iter, Np, tune_lmpc, ei, F0_idx, R0_idx, con_idx, lift_fn);
fprintf('  done (%.1f s wall)\n', toc(t0));

% metrics, same definitions as run_comparison; the 47-lift rows and the two
% plant-based baselines come from the saved comparison.mat for reference
C = load(fullfile(here, 'comparison.mat'), 'M');
names = {'KPC-spatial', 'Koopman-LMPC-spatial'};
R = {res_kpc, res_klmpc};
fprintf('\nDictionary study, closed loop (stressed mdot=%.2f, Np=%d):\n', mdot_scale, Np);
fprintf('%-22s  worst-met%%  unmet[kWh]  med-solve  feas\n', 'controller');
M = struct();
for j = 1:numel(R)
    r = R{j};
    pc  = 100 * sum(r.c_i, 2) ./ max(sum(r.d_i, 2), 1e-9);
    unm = sum(max(r.d_i - r.c_i, 0), 'all') * p.Ts / 3.6e6;
    feas = sum(r.exitflag > 0);
    M.(matlab.lang.makeValidName(names{j})) = struct('percons', pc, 'worst', min(pc), ...
        'unmet', unm, 'med', median(r.solve_ms), 'mx', max(r.solve_ms), 'feas', feas, 'n', numel(r.exitflag));
    fprintf('%-22s  %8.3f  %9.2f  %8.0f  %d/%d\n', ...
        names{j}, min(pc), unm, median(r.solve_ms), feas, numel(r.exitflag));
end
fprintf('%-22s  %8.3f  (47-lift, from comparison.mat)\n', 'KPC',          C.M.KPC.worst);
fprintf('%-22s  %8.3f  (47-lift, from comparison.mat)\n', 'Koopman-LMPC', C.M.Koopman_LMPC.worst);
fprintf('%-22s  %8.3f  (lift-independent, from comparison.mat)\n', 'NMPC',          C.M.NMPC.worst);
fprintf('%-22s  %8.3f  (lift-independent, from comparison.mat)\n', 'Jacobian-LMPC', C.M.Jacobian_LMPC.worst);
M.KPC_ref47        = C.M.KPC;
M.Koopman_ref47    = C.M.Koopman_LMPC;
M.NMPC_ref         = C.M.NMPC;
M.Jacobian_ref     = C.M.Jacobian_LMPC;

save(fullfile(here, 'comparison_spatial.mat'), 'res_kpc', 'res_klmpc', ...
     'names', 'M', 'mdot_scale', 'sc_start', 'T_sim', 'Np', 'p', 'ei', '-v7.3');
fprintf('\nSaved comparison_spatial.mat\n');
