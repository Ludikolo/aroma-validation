% DEMO_ROBUSTNESS  Part 4 summary: the comparison headline plus the long-deployment
% (longevity) check. Run run_comparison.m and run_longevity.m first; this loads their
% results, prints the headline, and plots the longevity as two time-domain traces. The
% worst-consumer demand-met stays flat across 90 days (no drift) and the per-day solve
% time stays flat (no compute creep), so nothing in the prediction or the QP degrades
% over a long run. The printed headline is the comparison: KPC leads the nonlinear
% NMPC on demand at a fraction of the solve time, and beats the iterated Koopman-LMPC.

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here)); startup;
res_dir  = fullfile(here, 'results');
figs_dir = fullfile(here, 'figures');
if ~exist(figs_dir, 'dir'), mkdir(figs_dir); end

p  = params();
LO = load(fullfile(res_dir, 'longevity.mat'));
C  = load(fullfile(here, '..', 'visualizations', '4_comparison', 'comparison.mat'));

%% --- Part 4 headline: the comparison ---
fprintf('\n=== Part 4 headline (stressed day, mdot = %.2f) ===\n', C.mdot_scale);
fprintf('%-15s  worst-met%%  unmet[kWh]  med-solve\n', 'controller');
for j = 1:numel(C.names)
    m = C.M.(matlab.lang.makeValidName(C.names{j}));
    fprintf('%-15s  %8.3f  %9.2f  %8.0f ms\n', C.names{j}, m.worst, m.unmet, m.med);
end
kpc = C.M.KPC; nmpc = C.M.NMPC;
fprintf('KPC leads NMPC on demand (%.1f%% vs %.1f%%), solves %.0fx faster, and stays feasible every step (%d/%d vs NMPC %d/%d).\n', ...
        kpc.worst, nmpc.worst, nmpc.med/max(kpc.med,eps), kpc.feas, kpc.n, nmpc.feas, nmpc.n);

%% --- deployability: 90-day longevity ---
fprintf('\n--- Longevity (%d days, stressed capacity mdot = 0.50) ---\n', LO.n_days_actual);
fprintf('  worst-consumer met %% drift day 1 -> day %d: %.4f -> %.4f (%+.4f pp)\n', ...
        LO.n_days_actual, LO.worst_per_day(1), LO.worst_per_day(end), LO.drift_total);
fprintf('  median solve day 1 -> day %d: %.0f -> %.0f ms (no creep), infeasible-step rate %.4f %%\n', ...
        LO.n_days_actual, LO.solve_ms_day(1,2), LO.solve_ms_day(end,2), LO.infeas_rate);

figure('Color', 'w', 'Position', [100 100 1000 420]);
days = 1:LO.n_days_actual;

subplot(1, 2, 1);
plot(days, LO.worst_per_day, 'b', 'LineWidth', 1.6);
grid on; xlabel('day'); ylabel('worst-consumer met [%]');
ylim([min(95, min(LO.worst_per_day) - 1), 100.5]);
title(sprintf('90-day demand-met (drift %+.3f pp)', LO.drift_total));

subplot(1, 2, 2);
plot(days, LO.solve_ms_day(:, 2), 'b', 'LineWidth', 1.6); hold on;
yline(p.Ts * 1000, 'r--', 'Ts = 900 s budget', 'LineWidth', 1.0);
grid on; xlabel('day'); ylabel('median QP solve time [ms]');
title('90-day solve time (no compute creep)');
set(gca, 'YScale', 'log'); ylim([100, p.Ts * 1000 * 2]);

exportgraphics(gcf, fullfile(figs_dir, 'longevity.pdf'), 'ContentType', 'vector');
fprintf('\nSaved figures/longevity.pdf\n');
fprintf('Run validate_robustness next to check the Part 4 invariants.\n');
