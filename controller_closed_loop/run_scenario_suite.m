% RUN_SCENARIO_SUITE  Six short closed-loop scenarios across the day at design
% capacity. Each runs three controllers on the same warmed-up plant and the same
% demand realisation: hold-nominal, demand-following RBC, and KPC. Saves all runs
% to controller_closed_loop/results/scenario_suite.mat (consumed by demo_controller,
% whose saved demo.mat is what validate_controller then checks).
%
%   S1_morning   07:00 start (residential morning peak)
%   S2_midday    12:30 start (commercial midday peak)
%   S3_evening   18:00 start (residential evening peak)
%   S4_night     02:00 start (low-demand sanity)
%   S5_shop_pm   16:30 start (shop afternoon peak)
%   S6_shoulder  21:00 start (evening-to-night shoulder)
%
% Each loop is 4 h = 16 samples at Ts = 900 s, so it spans a full peak event.

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root); startup();

p = params();
pred_dir = fullfile(root, 'predictor_open_loop', 'results');
ctrl_dir = fullfile(root, 'controller_closed_loop', 'results');

T_sim  = 4 * 3600;     % 4 h closed loop (16 samples at Ts = 900)
T_warm = 30 * 60;

scenarios = struct( ...
    'name',  {'S1_morning', 'S2_midday', 'S3_evening', 'S4_night', 'S5_shop_pm', 'S6_shoulder'}, ...
    'start', { 7*3600,       12.5*3600,   18*3600,      2*3600,     16.5*3600,    21*3600     });

% predictor maps + the locked tune (same build as demo_controller)
F  = load(fullfile(pred_dir, 'vseq_fits_full.mat'));
Hb = load(fullfile(ctrl_dir, 'best_horizons.mat'));
Ab = load(fullfile(ctrl_dir, 'best_alpha.mat'));
Bb = load(fullfile(ctrl_dir, 'best_tune.mat'));
pred = truncate_pred_to_Np(build_kpc_v2_matrices(F.fits), Hb.Np_best);
fprintf('Tune: Np=%d, Nc=%d, alpha=%g\n', Hb.Np_best, Hb.Nc_best, Ab.alpha_star);

tune = Bb.best_tune;
tune.dT_0s_max = 7.5; tune.dr_q_max = 4.5; tune.dT_ir_max = 15.0;   % actuator rate caps
tune.use_kirchhoff = true; tune.use_Tr_le_Ts = true; tune.use_mixing = true; tune.rho_slack_mix = 1;
tune.Nc = Hb.Nc_best; tune.alpha_energy = Ab.alpha_star;

[net, z0_cold] = build_plant(p);
ei = edge_user_index(net); n_user = ei.n_user;
F0_idx = find(strcmp({net.Nodes.name}, 'F0'));
R0_idx = find(strcmp({net.Nodes.name}, 'R0'));
con_idx = arrayfun(@(i) find(strcmp({net.Nodes.name}, ei.consumers{i})), 1:n_user);

results = struct();
for si = 1:numel(scenarios)
    sc = scenarios(si);
    fprintf('\n--- %s (start %.1f h) ---\n', sc.name, sc.start/3600);

    % warm the plant up at this time of day, then start the loop
    p_wu = p; p_wu.t_offset = sc.start;
    res_wu = simulate_plant(net, z0_cold, p_wu, @(t) 0, @(t) 1.0, T_warm);
    cl_start = sc.start + T_warm;

    % hold-nominal: open loop at nominal flow, supply at Tin_nom. this leg keeps
    % simulate_plant's extra boundary sample (one more column than the kpc/rbc
    % legs); demo_controller drops it when scoring, so the legs line up there
    p_hold = p; p_hold.t_offset = cl_start;
    res_hold = simulate_plant(net, res_wu.z_final, p_hold, @(t) 0, @(t) 1.0, T_sim);
    res_hold.T_0s = p.Tin_nom + res_hold.u;          % enrich so all legs carry the same fields
    res_hold.T_0r = res_hold.Tout(R0_idx, :);

    % demand-following rule-based baseline
    p_rbc = p; p_rbc.scenario_start = cl_start;
    res_rbc = run_rbc_demand_following(net, res_wu.z_final, p_rbc, T_sim);

    % KPC closed loop
    res_kpc = kpc_step_loop(net, res_wu, p, sc.start, T_warm, T_sim, ...
                            pred, tune, ei, F0_idx, R0_idx, con_idx);

    results.(sc.name).hold    = res_hold;
    results.(sc.name).rbc     = res_rbc;
    results.(sc.name).kpc     = res_kpc;
    results.(sc.name).start_h = sc.start / 3600;
end

save(fullfile(ctrl_dir, 'scenario_suite.mat'), 'results', 'tune', 'scenarios');
fprintf('\nSaved scenario_suite.mat (Np=%d)\n', Hb.Np_best);

% summary table (suite-mean met %, KPC solve time per scenario)
fprintf('\n%-15s  %-9s  %-9s  %-9s  %s\n', 'scenario', 'hold met%', 'rbc met%', 'kpc met%', 'kpc solve');
for si = 1:numel(scenarios)
    nm = scenarios(si).name;
    fprintf('%-15s  %-9.2f  %-9.2f  %-9.2f  %.0f ms\n', nm, ...
        suite_met(results.(nm).hold, n_user), suite_met(results.(nm).rbc, n_user), ...
        suite_met(results.(nm).kpc, n_user), mean(results.(nm).kpc.solve_ms));
end

function m = suite_met(res, n_user)
% suite-mean met % = total delivered / total demanded over the loop
sd = 0; sc = 0;
for i = 1:n_user
    sd = sd + sum(res.d_i(i, :));
    sc = sc + sum(res.c_i(i, :));
end
m = 100 * sc / max(sd, 1e-9);
end
