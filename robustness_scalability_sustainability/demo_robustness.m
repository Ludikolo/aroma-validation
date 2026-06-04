% DEMO_ROBUSTNESS  Closed-loop validation summary for the locked KPC v2
%   controller across 15 robustness / scalability / sustainability axes.
%   Unlike demo_plant / demo_predictor / demo_controller, this part is
%   a meta-analysis: each individual axis (capacity, sensor noise,
%   long horizon, climate shift, ...) has its own multi-hour MC or
%   sweep that is not feasible to rerun in a 30-second demo. The
%   committed .mat artefacts come from the full thesis matrix; this
%   demo loads them, computes the headline numbers, plots two summary
%   figures, and saves demo.mat for the validation step.
%
%   What this demo does (~5 s wall):
%     1. loads sustainability_headline.mat (4 KPIs: source-energy
%        ratio, T_0r band compliance, low-exergy Delta-T, delivery
%        ratio)
%     2. loads monte_carlo_stressed.mat (30-seed paired Monte Carlo
%        at stressed capacity, mdot = 0.6 x design)
%     3. computes the headline statistics (Cohen's d on unmet,
%        paired-t p-value, Wilcoxon signed-rank p, bootstrap CI)
%     4. plots a 4-bar sustainability summary + a paired-MC boxplot
%     5. saves demo.mat
%
%   Validation of nine invariants on the saved demo.mat is in
%   tests/validate_robustness.m.

clear; clc;
startup;
p = params(); %#ok<NASGU>

root     = fileparts(mfilename('fullpath'));
res_dir  = fullfile(root, 'results');
figs_dir = fullfile(root, 'figures');
if ~exist(figs_dir, 'dir'), mkdir(figs_dir); end

fprintf('\n=== Robustness / scalability / sustainability validation summary ===\n');

%% Sustainability headline (T17 acceptance bars)
SH = load(fullfile(res_dir, 'sustainability_headline.mat'));
fprintf('\n--- Sustainability headline (KPC v2 vs hold-nominal, suite) ---\n');
fprintf('  suite source energy   hold = %7.1f MJ   kpc = %7.1f MJ   saving = %+.2f %%\n', ...
        SH.suite_E_plant_hold, SH.suite_E_plant_kpc, SH.energy_saving_pct);
fprintf('  suite delivered       hold = %7.1f MJ   kpc = %7.1f MJ   ratio  =  %.3f\n', ...
        SH.suite_E_delivered_hold, SH.suite_E_delivered_kpc, SH.delivery_ratio);
fprintf('  suite Delta T source  hold = %7.2f K    kpc = %7.2f K\n', ...
        SH.suite_deltaT_hold, SH.suite_deltaT_kpc);
fprintf('  T_0r band compliance per scenario, kpc:  min = %5.2f %% (bar >= 90 %%)\n', ...
        min(SH.band_kpc));

%% Stressed Monte Carlo headline (T10b acceptance)
MC = load(fullfile(res_dir, 'monte_carlo_stressed.mat'));
suite_kpc  = mean(MC.kpc.met_pct(:));
suite_hold = mean(MC.hold_ctl.met_pct(:));
suite_rbc  = mean(MC.rbc.met_pct(:));

unmet_best_base = min([MC.hold_ctl.unmet_total, MC.rbc.unmet_total], [], 2);
delta_unmet = unmet_best_base - MC.kpc.unmet_total;
cohens_d    = mean(delta_unmet) / max(std(delta_unmet), eps);
[~, p_t]    = ttest(MC.kpc.unmet_total, unmet_best_base, 'Tail', 'left');
p_wilcox    = signrank(MC.kpc.unmet_total, unmet_best_base, 'tail', 'left');
rng_state = rng; rng(42);
d_ci = bootci(1000, @(x) mean(x)/max(std(x),eps), delta_unmet);
rng(rng_state);

fprintf('\n--- Stressed Monte Carlo (mdot = %.2f x design, %d seeds) ---\n', ...
        MC.mdot_scale, MC.n_seeds);
fprintf('  suite-mean met %%:  kpc = %.3f, hold = %.3f, rbc = %.3f\n', ...
        suite_kpc, suite_hold, suite_rbc);
fprintf('  paired t-test  KPC < best(hold, rbc) on unmet : p = %.4g\n', p_t);
fprintf('  Wilcoxon signed-rank p                         : %.4g\n', p_wilcox);
fprintf('  Cohen''s d on unmet (KPC vs best)               : %.3f  CI = [%+.2f, %+.2f]\n', ...
        cohens_d, d_ci(1), d_ci(2));
fprintf('  mean QP solve %.0f ms / max %.0f ms (Ts budget 900 s)\n', ...
        mean(MC.kpc.solve_ms_mean), max(MC.kpc.solve_ms_max));

%% Deployment headline (scalability + longevity)
% These two artefacts come from the larger plant-in-the-loop runs that are
% too slow to redo in a demo (a 50-consumer network and a 90-day loop); the
% saved .mat hold the headline results and the validator re-checks them.
N50 = load(fullfile(res_dir, 'n50_closed_loop.mat'));
n50_worst  = N50.result.worst_cons;
n50_med_ms = N50.result.median_ms;
D90 = load(fullfile(res_dir, 'longevity_90day.mat'));
drift_90 = D90.drift_total;
n90_days = D90.n_days;
fprintf('\n--- Deployment headline ---\n');
fprintf('  n_user = 50 plant-in-the-loop: worst-consumer met %% = %.2f %%, median solve %.0f ms\n', ...
        n50_worst, n50_med_ms);
fprintf('  %g-day closed loop: worst-consumer drift = %.3g pp (zero-drift longevity)\n', ...
        n90_days, drift_90);

%% Plot 1: 4-bar sustainability summary
fig = figure('Visible', 'off', 'Position', [100 100 900 360]);
subplot(1, 2, 1);
bar([SH.suite_E_plant_hold, SH.suite_E_plant_kpc; ...
     SH.suite_E_delivered_hold, SH.suite_E_delivered_kpc]);
set(gca, 'XTick', 1:2, 'XTickLabel', {'E\_plant', 'E\_delivered'});
ylabel('MJ over 6-scenario suite');
title(sprintf('Energy balance (saving %+.2f %%)', SH.energy_saving_pct));
legend('hold', 'KPC', 'Location', 'northeast'); grid on;

subplot(1, 2, 2);
metrics = [SH.suite_deltaT_kpc, SH.delivery_ratio, min(SH.band_kpc) / 100, ...
           (SH.suite_E_plant_hold - SH.suite_E_plant_kpc) / SH.suite_E_plant_hold];
bars  = [15, 0.90, 0.90, 0];
labels = {'\Delta T_{src}', 'delivery', 'band frac', 'E save'};
plot_metric = metrics;
plot_bar    = bars;
b = bar([plot_metric; plot_bar]');
set(gca, 'XTick', 1:4, 'XTickLabel', labels);
ylabel('value');
title('Sustainability bars (T17, KPC)');
legend('value', 'bar', 'Location', 'northeast'); grid on;
print(fig, fullfile(figs_dir, 'sustainability.pdf'), '-dpdf', '-bestfit');
print(fig, fullfile(figs_dir, 'sustainability.png'), '-dpng', '-r150');
close(fig);

%% Plot 2: stressed-MC boxplot KPC vs hold vs RBC
fig = figure('Visible', 'off', 'Position', [100 100 700 360]);
suite_per_seed_kpc  = mean(MC.kpc.met_pct, 2);
suite_per_seed_hold = mean(MC.hold_ctl.met_pct, 2);
suite_per_seed_rbc  = mean(MC.rbc.met_pct, 2);
boxplot([suite_per_seed_kpc, suite_per_seed_rbc, suite_per_seed_hold], ...
        'Labels', {'KPC', 'RBC', 'hold'});
yline(95, 'r--', '95 %% bar');
ylabel('Suite-mean met % per seed');
title(sprintf('Stressed Monte Carlo (mdot = %.2f x design, %d seeds, Cohen''s d = %.2f)', ...
              MC.mdot_scale, MC.n_seeds, cohens_d));
grid on;
print(fig, fullfile(figs_dir, 'mc_boxplot.pdf'), '-dpdf', '-bestfit');
print(fig, fullfile(figs_dir, 'mc_boxplot.png'), '-dpng', '-r150');
close(fig);

%% Save demo artefact
save(fullfile(res_dir, 'demo.mat'), ...
     'SH', 'MC', ...
     'suite_kpc', 'suite_hold', 'suite_rbc', ...
     'cohens_d', 'p_t', 'p_wilcox', 'd_ci', ...
     'delta_unmet', 'n50_worst', 'n50_med_ms', 'drift_90', 'n90_days', '-v7.3');

fprintf('\nSaved %s\n', fullfile(res_dir, 'demo.mat'));
fprintf('Figures: sustainability.{pdf,png}, mc_boxplot.{pdf,png}\n');
fprintf('Run validate_robustness next to check the nine invariants.\n');
