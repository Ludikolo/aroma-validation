# Robustness, scalability, sustainability: validation reference

## What this part does

The locked KPC v2 controller from `controller_closed_loop/` is
validated as deployment-ready along three axes:

  - **robustness** under perturbations (sensor noise, plant-parameter
    mismatch, demand forecast noise, capacity scaling, ...)
  - **scalability** to bigger districts (compute scales sub-quadratic
    in `n_user`, projects to 200-consumer deployments)
  - **sustainability** measured by 11 district-heating KPIs from the
    Buffa / Wirtz 5GDHC literature, with four acceptance bars

This part is a meta-analysis: every individual axis is its own
multi-hour Monte Carlo, parameter sweep, or long-horizon CL run.
The 30-second `demo_robustness.m` loads the committed `.mat`
artefacts, computes headline statistics, plots two summary figures,
and saves `demo.mat` for the validation step.

## 15 robustness axes (T9..T16 + Bucket-B extensions)

| # | axis                                  | phase test                             |
|---|---------------------------------------|----------------------------------------|
| 1 | capacity scaling 0.5--1.5x            | `phase_t9_robustness.m`                |
| 2 | sensor noise sigma in [0, 0.4]        | `phase_t9_robustness.m`                |
| 3 | warmup duration 15--120 min            | `phase_t9_robustness.m`                |
| 4 | plant-model mismatch (alpha, eps)     | `phase_t9d_mismatch.m`                 |
| 5 | demand forecast noise sigma            | `phase_t9e_demand_forecast.m`          |
| 6 | per-consumer demand heterogeneity     | `phase_t9f_heterogeneous_demand.m`     |
| 7 | design-capacity Monte Carlo (30 seeds) | `phase_t10_mc.m`                       |
| 8 | stressed-capacity Monte Carlo          | `phase_t10b_mc_stressed.m`             |
| 9 | sampling-time off-design               | `phase_t11_ts_sweep.m`                 |
| 10| multi-rate train-deploy                | `phase_t11b_multirate.m` / t11d        |
| 11| slack-bound theorem (post-hoc safety) | `phase_t11c_slack_bound.m`             |
| 12| out-of-distribution profile-swap      | `phase_t13_ood_profile_swap.m`         |
| 13| compute scaling in n_user             | `phase_t14_compute_scaling.m`          |
| 14| climate / ambient temperature         | `phase_t15_climate_shift.m`            |
| 15| long-horizon drift (7 days + 30 days) | `phase_t16_long_horizon.m` + B1        |

Bucket-B (this final pass) adds five non-overlapping extensions on
top of the 15 axes above: 30-day horizon, 20 random profile
permutations, extreme `T_ext` (-20..+25 C), source-state z_0 init
perturbation, and QP failure-mode recovery. Combined acceptance is
asserted by `phase_t18_bucket_b.m`.

## Sustainability headline (KPC v2 vs hold-nominal, 6-scenario suite)

| metric                              | KPC     | hold    | delta / bar                |
|-------------------------------------|---------|---------|----------------------------|
| suite source energy `E_plant`        | 3508 MJ | 3756 MJ | **-6.61 %**, bar: KPC<hold |
| suite delivered `E_delivered`        | 1698 MJ | 1816 MJ | -6.50 %                    |
| delivery ratio (KPC / hold)         | 0.935   | 1.000   | bar: >= 0.90               |
| per-scenario `T_0r` in [10, 25] band | 100 %   | 100 %   | bar: >= 90 %               |
| suite-mean `dT_source`              | 3.13 K  | 3.26 K  | bar: <= 15 K               |

The four bars are documented in `phase_t17_sustainability.m` in the
broader thesis codebase; the `validate_robustness.m` in this repo
asserts them against the precomputed `sustainability_headline.mat`.

## Headline numbers under stress (30-seed paired Monte Carlo)

`mdot_scale = 0.6 x design`, 30 seeds, randomised warmup + initial-
condition perturbation:

| controller        | suite-mean met % | mean solve [ms] |
|-------------------|------------------|-----------------|
| **KPC v2 locked** | **100.00**       | **95**          |
| RBC               | 97.95            | n/a             |
| hold-nominal      | 90.71            | n/a             |

Statistics on `unmet_total` (lower = better, KPC vs `min(hold, rbc)`):
  - paired t-test, one-sided: p ~ 1e-45
  - Wilcoxon signed-rank   :   all 30 seeds favour KPC (p well below 1e-3)
  - Cohen's d              :   31.8   95% bootstrap CI ~ [25.5, 37.7]

## Checks (validate_robustness.m)

| # | invariant                                                                       | result                          |
|---|---------------------------------------------------------------------------------|---------------------------------|
| 1 | source-energy saving > 0 (KPC < hold)                                            | +6.61 %                        |
| 2 | delivered ratio kpc/hold >= 0.90                                                | 0.935                          |
| 3 | suite-mean `dT_source` <= 15 K (5GDHC low-exergy)                                | 3.13 K                         |
| 4 | per-scenario `T_0r` band compliance >= 90 %                                      | min 100 %                      |
| 5 | stressed MC KPC suite-mean >= 95 %                                              | 100.0 %                        |
| 6 | paired t-test KPC < best(hold, RBC) on unmet, p < 1e-3                          | p ~ 1e-45                      |
| 7 | Cohen's d on unmet >= 0.5 with 95 % bootstrap CI > 0                            | d = 31.8, CI = [25.5, 37.7]    |
| 8 | scalability: n_user = 50 worst-consumer met % >= 95 %                           | 99.89 % (median solve 30.5 s)  |
| 9 | longevity: 90-day worst-consumer drift ~ 0                                       | 0 pp                           |

## Files in this folder

  - `ROBUSTNESS_NOTES.md`              this reference
  - `demo_robustness.m`                meta-analysis demo, ~5 s
  - `results/sustainability_headline.mat`  4-bar KPI baseline
  - `results/monte_carlo_stressed.mat`  30-seed stressed MC
  - `results/n50_closed_loop.mat`      50-consumer plant-in-the-loop result
  - `results/longevity_90day.mat`      90-day closed-loop drift result
  - `results/demo.mat`                 written by demo, read by validate
  - `figures/sustainability.{pdf,png}`  4-bar summary chart
  - `figures/mc_boxplot.{pdf,png}`     paired-MC box plot

The locked tune artefacts (`best_tune.mat`, `best_horizons.mat`,
`best_alpha.mat`) live in `controller_closed_loop/results/` and the
predictor fits (`vseq_fits_full.mat`) live in
`predictor_open_loop/results/`; this folder cross-references them
rather than duplicating.
