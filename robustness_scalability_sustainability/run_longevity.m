% RUN_LONGEVITY  90-day continuous closed loop, to check the controller does not drift
% over a long deployment. One single kpc_step_loop call (no day-by-day stitching, so
% no boundary artefacts), then aggregate per day. The same locked KPC controller runs
% for 90 simulated days at the stressed half-capacity (mdot 0.50); that is ~8640
% control steps, so any accumulating error (numerical, model, or feasibility) would
% show up as the per-day demand-met drifting down or the solve creeping up. A flat
% met % line and a flat solve-time line mean it stays stable for long-horizon use.
%
% Output: results/longevity.mat (met_pct_day, suite_per_day, solve_ms_day, drift).

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root); startup();

p = params();
pred_dir = fullfile(root, 'predictor_open_loop', 'results');
ctrl_dir = fullfile(root, 'controller_closed_loop', 'results');

% --- KPC controller (predictor maps + locked tune) ---
F  = load(fullfile(pred_dir, 'vseq_fits_full.mat'));
Hb = load(fullfile(ctrl_dir, 'best_horizons.mat'));
Ab = load(fullfile(ctrl_dir, 'best_alpha.mat'));
Bb = load(fullfile(ctrl_dir, 'best_tune.mat'));
pred = truncate_pred_to_Np(build_kpc_v2_matrices(F.fits), Hb.Np_best);
tune = Bb.best_tune;
tune.dT_0s_max = 7.5; tune.dr_q_max = 4.5; tune.dT_ir_max = 15.0;
tune.use_kirchhoff = true; tune.use_Tr_le_Ts = true; tune.use_mixing = true; tune.rho_slack_mix = 1;
tune.Nc = Hb.Nc_best; tune.alpha_energy = Ab.alpha_star;

% --- plant at stressed capacity (mdot 0.50, the hard case), so a no-drift result
% means stable under binding load, not just at an easy design point ---
T_warm = 30 * 60;
n_days = 90;
T_sim  = n_days * 86400;
mdot_scale = 0.50;
[net, z0_cold] = build_plant(p);
net.mdotEdges = mdot_scale * net.mdotEdges; net.flows = net.mdotEdges; net.q0 = net.mdotEdges(:);
ei = edge_user_index(net); n_user = ei.n_user;
F0_idx = find(strcmp({net.Nodes.name}, 'F0'));
R0_idx = find(strcmp({net.Nodes.name}, 'R0'));
con_idx = arrayfun(@(i) find(strcmp({net.Nodes.name}, ei.consumers{i})), 1:n_user);

p_wu = p; p_wu.t_offset = 0; p_wu.r_q_fun = @(t) net.mdotEdges(:);
res_wu = simulate_plant(net, z0_cold, p_wu, @(t) 0, @(t) 1.0, T_warm);
fprintf('Warmup done. Starting %d-day closed loop ...\n', n_days);

t0 = tic;
% under the extreme bind the QP can very rarely fail to converge; the controller then
% holds its previous solution and recovers (we report that infeasible-step rate below),
% so silence the solver's per-failure warning to keep this long run's log readable
ws = warning('off', 'all');
res = kpc_step_loop(net, res_wu, p, 0, T_warm, T_sim, pred, tune, ei, F0_idx, R0_idx, con_idx);
warning(ws);
fprintf('%d-day closed loop done in %.1f min\n', n_days, toc(t0)/60);

% --- per-day aggregation ---
steps_per_day = round(86400 / p.Ts);
N_cl = size(res.d_i, 2);
n_days_actual = floor(N_cl / steps_per_day);
met_pct_day  = nan(n_days_actual, n_user);
solve_ms_day = nan(n_days_actual, 3);
infeas_day   = zeros(n_days_actual, 1);
for d = 1:n_days_actual
    idx = (d-1)*steps_per_day + 1 : d*steps_per_day;
    for i = 1:n_user
        met_pct_day(d, i) = 100 * sum(res.c_i(i, idx)) / max(sum(res.d_i(i, idx)), 1e-9);
    end
    solve_ms_day(d, :) = [min(res.solve_ms(idx)), median(res.solve_ms(idx)), max(res.solve_ms(idx))];
    infeas_day(d) = sum(res.exitflag(idx) <= 0);
end
worst_per_day = min(met_pct_day, [], 2);
suite_per_day = mean(met_pct_day, 2);
drift_total   = worst_per_day(end) - worst_per_day(1);
infeas_rate   = sum(infeas_day) / (n_days_actual * steps_per_day) * 100;

fprintf('\n=== %d-day summary ===\n', n_days_actual);
fprintf('  worst-consumer met%%  day 1 : %.4f %%\n', worst_per_day(1));
fprintf('  worst-consumer met%%  day %d: %.4f %%\n', n_days_actual, worst_per_day(end));
fprintf('  total drift          : %+.4f pp\n', drift_total);
fprintf('  infeasible-step rate : %.4f %%\n', infeas_rate);
fprintf('  median solve (day 1 / day %d): %.0f / %.0f ms\n', ...
        n_days_actual, solve_ms_day(1, 2), solve_ms_day(end, 2));

save(fullfile(here, 'results', 'longevity.mat'), ...
     'n_days_actual', 'met_pct_day', 'worst_per_day', 'suite_per_day', ...
     'solve_ms_day', 'infeas_day', 'drift_total', 'infeas_rate', '-v7.3');
fprintf('\nSaved results/longevity.mat\n');
