% RUN_STRESS_SWEEP  Run the two fast Koopman controllers (KPC direct multi-step
% vs the iterated one-step Koopman-LMPC) at progressively harder flow capacity
% and save the worst-consumer demand-met. As the bind tightens the gap opens up:
% KPC stays above the comfort bar while the iterated model drops below it. That
% is why the head-to-head is run at a hard point. Render with viz_stress_sweep.

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(root); addpath(fullfile(root, 'controller_closed_loop', 'lmpc')); startup();

p = params();
pred_dir = fullfile(root, 'predictor_open_loop', 'results');
ctrl_dir = fullfile(root, 'controller_closed_loop', 'results');
F  = load(fullfile(pred_dir, 'vseq_fits_full.mat'));
Hb = load(fullfile(ctrl_dir, 'best_horizons.mat'));
Ab = load(fullfile(ctrl_dir, 'best_alpha.mat'));
Bb = load(fullfile(ctrl_dir, 'best_tune.mat'));
pred = truncate_pred_to_Np(build_kpc_v2_matrices(F.fits), Hb.Np_best);
tune = Bb.best_tune;
tune.dT_0s_max = 7.5; tune.dr_q_max = 4.5; tune.dT_ir_max = 15.0;
tune.use_kirchhoff = true; tune.use_Tr_le_Ts = true; tune.use_mixing = true; tune.rho_slack_mix = 1;
tune.Nc = Hb.Nc_best; tune.alpha_energy = Ab.alpha_star;
Np = Hb.Np_best;
S = load(fullfile(pred_dir, 'iterated_AB.mat'));
if isfield(S, 'fit_iter'), fit_iter = S.fit_iter; else, fit_iter = S; end
tune_lmpc = tune; tune_lmpc.use_mixing = false;

T_warm = 30*60; T_sim = 24*3600; sc = 0;
mdots = [0.50 0.45 0.40 0.35 0.30];
met_kpc  = zeros(size(mdots));
met_lmpc = zeros(size(mdots));
fprintf('\n mdot |  KPC worst |  iterated worst |  gap\n');
for mi = 1:numel(mdots)
    md = mdots(mi);
    [net, z0] = build_plant(p);
    net.mdotEdges = md*net.mdotEdges; net.flows = net.mdotEdges; net.q0 = net.mdotEdges(:);
    ei = edge_user_index(net); n_user = ei.n_user;
    F0 = find(strcmp({net.Nodes.name}, 'F0')); R0 = find(strcmp({net.Nodes.name}, 'R0'));
    con = arrayfun(@(i) find(strcmp({net.Nodes.name}, ei.consumers{i})), 1:n_user);
    pwu = p; pwu.t_offset = 0; pwu.r_q_fun = @(t) net.mdotEdges(:);
    rwu = simulate_plant(net, z0, pwu, @(t)0, @(t)1.0, T_warm);
    rk = kpc_step_loop(net, rwu, p, sc, T_warm, T_sim, pred, tune, ei, F0, R0, con);
    rl = lmpc_step_loop(net, rwu, p, sc, T_warm, T_sim, fit_iter, Np, tune_lmpc, ei, F0, R0, con);
    met_kpc(mi)  = min(100 * sum(rk.c_i, 2) ./ max(sum(rk.d_i, 2), 1e-9));
    met_lmpc(mi) = min(100 * sum(rl.c_i, 2) ./ max(sum(rl.d_i, 2), 1e-9));
    fprintf(' %.2f | %10.3f | %15.3f | %+6.3f\n', md, met_kpc(mi), met_lmpc(mi), met_kpc(mi) - met_lmpc(mi));
end

save(fullfile(here, 'stress_sweep.mat'), 'mdots', 'met_kpc', 'met_lmpc');
fprintf('\nSaved stress_sweep.mat (render it with viz_stress_sweep)\n');
