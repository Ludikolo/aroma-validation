% DEMO_CONTROLLER  Closed-loop demonstration and validation of the
%   convex-QP MPC controller (KPC v2) built on top of the V_seq
%   Koopman predictor from predictor_open_loop. The controller solves
%   one QP per sample step (Ts = 15 min) to track demand on every
%   consumer while respecting the substation contracts and the
%   network's mass-conservation constraints.
%
%   The locked-in tune (Np = 5, Nc = 1, alpha = 1, rho_slack = 1.195)
%   was selected by a cross-scenario Pareto sweep followed by
%   Bayesian optimisation over the slack penalty. See
%   CONTROLLER_NOTES.md for the formulation and tune-selection details.
%
%   What this demo does (~30 s wall):
%     1. loads the predictor fits (predictor_open_loop/results),
%        the locked tune, and the 6-scenario suite at design capacity.
%     2. builds the plant at stressed capacity (mdot_scale = 0.6),
%        the regime where the controller's advantage opens up.
%     3. runs one 24 h closed loop with KPC v2 + the hold-nominal
%        baseline on a single deterministic seed.
%     4. computes per-consumer demand met %, plots a tracking trace
%        plus a controller-comparison bar chart, and saves demo.mat.
%
%   Validation of seven controller invariants on the saved demo.mat is
%   in tests/validate_controller.m.

clear; clc;
startup;
p = params();

root       = fileparts(mfilename('fullpath'));
res_dir    = fullfile(root, 'results');
figs_dir   = fullfile(root, 'figures');
pred_res   = fullfile(fileparts(root), 'predictor_open_loop', 'results');
if ~exist(figs_dir, 'dir'), mkdir(figs_dir); end

fprintf('\n=== Controller closed-loop demo (KPC v2 vs hold-nominal) ===\n');
fprintf('Sample rate Ts = %g s, locked tune from cross-scenario Pareto + BO\n', p.Ts);

%% Load predictor fits + locked tune
F = load(fullfile(pred_res, 'vseq_fits_full.mat'));
pred_full = build_kpc_v2_matrices(F.fits);
H = load(fullfile(res_dir, 'best_horizons.mat'));
A = load(fullfile(res_dir, 'best_alpha.mat'));
B = load(fullfile(res_dir, 'best_tune.mat'));
pred = truncate_pred_to_Np(pred_full, H.Np_best);

tune = B.best_tune;
tune.dT_0s_max   = 7.5;
tune.dr_q_max    = 4.5;
tune.dT_ir_max   = 15.0;
tune.use_kirchhoff = true;
tune.use_Tr_le_Ts  = true;
tune.use_mixing    = true;
tune.rho_slack_mix = 1;
tune.Nc            = H.Nc_best;
tune.alpha_energy  = A.alpha_star;

fprintf('Locked tune: (Np, Nc, alpha, rho_slack) = (%d, %d, %g, %.2f)\n', ...
        H.Np_best, H.Nc_best, A.alpha_star, tune.rho_slack);

%% Build plant at stressed capacity (mdot_scale = 0.6)
mdot_scale = 0.6;
T_warm     = 30 * 60;
T_sim      = 24 * 3600;
sc_start   = 0;

[net, z0_cold] = build_plant(p);
net.mdotEdges = mdot_scale * net.mdotEdges;
net.flows     = net.mdotEdges;
net.q0        = net.mdotEdges(:);

ei = edge_user_index(net);
n_user = ei.n_user;
F0_idx = find(strcmp({net.Nodes.name}, 'F0'));
R0_idx = find(strcmp({net.Nodes.name}, 'R0'));
con_idx = zeros(n_user, 1);
for i = 1:n_user
    con_idx(i) = find(strcmp({net.Nodes.name}, ei.consumers{i}));
end

fprintf('Plant: 23 nodes, 29 edges, mdot_scale = %.2f (stressed)\n', mdot_scale);

%% 30 min warmup
fprintf('Warmup (%.0f min)...\n', T_warm / 60);
p_wu = p; p_wu.t_offset = sc_start;
res_wu = simulate_plant(net, z0_cold, p_wu, @(t) 0, @(t) 1.0, T_warm);

%% Closed loop: KPC v2
fprintf('Closed loop KPC v2 (%d h, %d steps)...\n', T_sim / 3600, round(T_sim / p.Ts));
t_kpc = tic;
res_kpc = kpc_step_loop(net, res_wu, p, sc_start, T_warm, T_sim, ...
                        pred, tune, ei, F0_idx, R0_idx, con_idx);
fprintf('  done in %.1f s, mean solve %.0f ms, max %.0f ms\n', ...
        toc(t_kpc), mean(res_kpc.solve_ms), max(res_kpc.solve_ms));

%% Hold-nominal baseline
fprintf('Closed loop hold-nominal (same plant, same warmup)...\n');
p_hold = p; p_hold.t_offset = sc_start + T_warm;
res_hold = simulate_plant(net, res_wu.z_final, p_hold, @(t) 0, @(t) 1.0, T_sim);

%% Per-consumer met %
met_kpc  = nan(1, n_user);
met_hold = nan(1, n_user);
for i = 1:n_user
    sd_k = sum(res_kpc.d_i(i, :));   sk_k = sum(res_kpc.c_i(i, :));
    sd_h = sum(res_hold.d_i(i, :));  sk_h = sum(res_hold.c_i(i, :));
    met_kpc(i)  = 100 * sk_k / max(sd_k, 1e-9);
    met_hold(i) = 100 * sk_h / max(sd_h, 1e-9);
end

fprintf('\nPer-consumer met %% (stressed-S2, 24 h, single seed):\n');
fprintf('%-4s %10s %10s\n', 'C', 'KPC', 'hold');
for i = 1:n_user
    fprintf('  %s %10.3f %10.3f\n', ei.consumers{i}, met_kpc(i), met_hold(i));
end
fprintf('%-4s %10.3f %10.3f   (suite mean)\n', '', mean(met_kpc), mean(met_hold));

%% Load 6-scenario suite at design capacity (headline result)
S = load(fullfile(res_dir, 'scenario_suite.mat'));
sc_names = fieldnames(S.results);
n_sc = numel(sc_names);
met_design_kpc  = nan(n_sc, n_user);
met_design_hold = nan(n_sc, n_user);
for si = 1:n_sc
    r = S.results.(sc_names{si});
    for i = 1:n_user
        sd_k = sum(r.kpc.d_i(i, :));   sk_k = sum(r.kpc.c_i(i, :));
        sd_h = sum(r.hold.d_i(i, :));  sk_h = sum(r.hold.c_i(i, :));
        met_design_kpc(si, i)  = 100 * sk_k / max(sd_k, 1e-9);
        met_design_hold(si, i) = 100 * sk_h / max(sd_h, 1e-9);
    end
end

fprintf('\nHeadline at design capacity (6 scenarios x 5 consumers):\n');
fprintf('  KPC  worst cell : %.3f %%   (only controller >= 95 %% on every cell)\n', ...
        min(met_design_kpc(:)));
fprintf('  hold worst cell : %.3f %%\n', min(met_design_hold(:)));

%% Plot 1: tracking trace for the worst-consumer-at-stressed
[~, worst_i] = min(met_kpc);
fig = figure('Visible', 'off', 'Position', [100 100 900 360]);
t_h_kpc  = (0:size(res_kpc.d_i, 2) -1) * p.Ts / 3600;
t_h_hold = (0:size(res_hold.c_i, 2)-1) * (T_sim / size(res_hold.c_i, 2)) / 3600;
plot(t_h_kpc,  res_kpc.d_i(worst_i, :)  / 1000, 'k-',  'LineWidth', 1.4, 'DisplayName', 'demand d^{(i)}');
hold on;
plot(t_h_kpc,  res_kpc.c_i(worst_i, :)  / 1000, 'b-',  'LineWidth', 1.2, 'DisplayName', 'KPC delivered c^{(i)}');
plot(t_h_hold, res_hold.c_i(worst_i, :) / 1000, 'r--', 'LineWidth', 1.0, 'DisplayName', 'hold delivered c^{(i)}');
grid on;
xlabel('Time [h]'); ylabel('Heat [kW]');
title(sprintf('Worst-consumer tracking at stressed capacity (C%d, met%% = %.2f)', ...
              worst_i, met_kpc(worst_i)));
legend('Location', 'best');
print(fig, fullfile(figs_dir, 'tracking.pdf'), '-dpdf', '-bestfit');
print(fig, fullfile(figs_dir, 'tracking.png'), '-dpng', '-r150');
close(fig);

%% Plot 2: met %% bar at stressed + design
fig = figure('Visible', 'off', 'Position', [100 100 900 360]);
subplot(1, 2, 1);
bar([met_design_kpc(:)'; met_design_hold(:)']', 'grouped');
ylim([min([met_design_kpc(:); met_design_hold(:)]) - 1, 101]);
yline(95, 'r--', '95 %% bar');
xlabel('Cell index (scenario x consumer)'); ylabel('met %');
title('Design capacity: 6 scenarios x 5 consumers');
legend({'KPC', 'hold'}, 'Location', 'southwest');
grid on;

subplot(1, 2, 2);
bar([met_kpc; met_hold]', 'grouped');
ylim([0, 102]);
yline(95, 'r--', '95 %% bar');
set(gca, 'XTick', 1:n_user, 'XTickLabel', ei.consumers);
xlabel('Consumer'); ylabel('met %');
title(sprintf('Stressed capacity (mdot = %.1f): 24 h', mdot_scale));
legend({'KPC', 'hold'}, 'Location', 'southwest');
grid on;
print(fig, fullfile(figs_dir, 'met_bar.pdf'), '-dpdf', '-bestfit');
print(fig, fullfile(figs_dir, 'met_bar.png'), '-dpng', '-r150');
close(fig);

%% Save demo artefact for validate_controller.m
save(fullfile(res_dir, 'demo.mat'), ...
     'mdot_scale', 'T_sim', 'T_warm', 'sc_start', ...
     'res_kpc', 'res_hold', 'met_kpc', 'met_hold', ...
     'met_design_kpc', 'met_design_hold', ...
     'tune', '-v7.3');

fprintf('\nSaved %s\n', fullfile(res_dir, 'demo.mat'));
fprintf('Figures: tracking.{pdf,png}, met_bar.{pdf,png}\n');
fprintf('Run validate_controller next to check the seven invariants.\n');
