% VALIDATE_ROBUSTNESS  Check the locked KPC v2 controller against
%   nine invariants on the saved meta-analysis demo.mat. This is
%   Part 4: robustness, scalability, sustainability across the
%   committed 30-seed Monte Carlo + sustainability KPI artefacts.
%
%   Run demo_robustness first. This script does not refit or rerun
%   anything; it operates on the saved sustainability_headline.mat
%   and monte_carlo_stressed.mat (loaded by the demo and re-saved
%   into demo.mat) and asserts seven acceptance bars.
%
%   T1  Suite source-energy saving > 0 (KPC < hold-nominal).
%   T2  Suite-mean delivery ratio (KPC / hold) >= 0.90.
%   T3  Suite-mean Delta T source <= 15 K (5GDHC low-exergy regime).
%   T4  Per-scenario T_0r band compliance >= 90 %% on every scenario.
%   T5  Stressed Monte Carlo: KPC suite-mean met %% >= 95 (above bar).
%   T6  Stressed Monte Carlo: paired one-sided t-test KPC < best
%       (hold, RBC) on total unmet, p < 1e-3.
%   T7  Stressed Monte Carlo: Cohen's d on unmet >= 0.5 with 95 %%
%       bootstrap CI excluding 0.
%   T8  Scalability: plant-in-the-loop at n_user = 50 keeps worst-
%       consumer met %% >= 95.
%   T9  Longevity: 90-day closed loop shows zero worst-consumer drift.

clear; clc;
startup;
p = params(); %#ok<NASGU>

results_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                       'robustness_scalability_sustainability', 'results');

S = load(fullfile(results_dir, 'demo.mat'));

fprintf('\n=== Robustness / scalability / sustainability validation: 7 invariants ===\n');

%% T1 suite source-energy saving > 0
saving = S.SH.energy_saving_pct;
assert(saving > 0, 'T1 fail: source-energy saving = %+.2f %% (KPC not below hold)', saving);
fprintf('T1 ok: source-energy saving KPC vs hold = %+.2f %% (> 0)\n', saving);

%% T2 delivery ratio >= 0.90
assert(S.SH.delivery_ratio >= 0.90, ...
    'T2 fail: delivery ratio = %.3f < 0.90', S.SH.delivery_ratio);
fprintf('T2 ok: delivered ratio kpc/hold = %.3f >= 0.90\n', S.SH.delivery_ratio);

%% T3 suite Delta T source <= 15 K
assert(S.SH.suite_deltaT_kpc <= 15, ...
    'T3 fail: suite mean Delta T_src = %.2f K > 15 K', S.SH.suite_deltaT_kpc);
fprintf('T3 ok: suite Delta T_src = %.2f K <= 15 K (5GDHC low-exergy)\n', ...
        S.SH.suite_deltaT_kpc);

%% T4 per-scenario T_0r band compliance >= 90 %
assert(min(S.SH.band_kpc) >= 90, ...
    'T4 fail: min band compliance = %.2f %% < 90 %%', min(S.SH.band_kpc));
fprintf('T4 ok: T_0r in [10, 25] C on every scenario (min = %.2f %%)\n', ...
        min(S.SH.band_kpc));

%% T5 stressed MC KPC suite-mean >= 95
assert(S.suite_kpc >= 95, ...
    'T5 fail: stressed MC KPC suite-mean = %.3f %% < 95 %%', S.suite_kpc);
fprintf('T5 ok: stressed MC KPC suite-mean = %.3f %% >= 95 %% bar\n', S.suite_kpc);

%% T6 paired t-test p < 1e-3
assert(S.p_t < 1e-3, ...
    'T6 fail: paired t-test p = %.4g >= 1e-3', S.p_t);
fprintf('T6 ok: paired t-test KPC < best baseline on unmet, p = %.4g\n', S.p_t);

%% T7 Cohen's d >= 0.5 with bootstrap CI excluding 0
assert(S.cohens_d >= 0.5, ...
    'T7 fail: Cohen''s d = %.3f < 0.5', S.cohens_d);
assert(S.d_ci(1) > 0, ...
    'T7 fail: bootstrap CI on d = [%+.2f, %+.2f] includes 0', S.d_ci(1), S.d_ci(2));
fprintf('T7 ok: Cohen''s d on unmet = %.3f (>= 0.5), CI = [%+.2f, %+.2f] (excludes 0)\n', ...
        S.cohens_d, S.d_ci(1), S.d_ci(2));

%% T8 scalability: plant-in-the-loop at n_user = 50 still clears the 95 % bar
assert(S.n50_worst >= 95, ...
    'T8 fail: n_user = 50 worst-consumer met %% = %.2f < 95 %%', S.n50_worst);
fprintf('T8 ok: n_user = 50 worst-consumer met %% = %.2f %% (>= 95 %%), median solve %.0f ms\n', ...
        S.n50_worst, S.n50_med_ms);

%% T9 longevity: 90-day closed loop shows no worst-consumer drift
assert(abs(S.drift_90) < 1e-6, ...
    'T9 fail: %g-day drift = %.3g pp (not zero)', S.n90_days, S.drift_90);
fprintf('T9 ok: %g-day worst-consumer drift = %.3g pp (zero-drift longevity)\n', ...
        S.n90_days, S.drift_90);

fprintf('\n=== Robustness validation: all 9 invariants passed ===\n');
