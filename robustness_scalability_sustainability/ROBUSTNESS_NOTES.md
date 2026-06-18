# Part 4: controller comparison + deployment checks

## What this part does

Two things, both shown as time-domain figures (no bar charts):

1. a head-to-head **comparison** of the KPC controller against three MPC baselines
   on a hard day, to show it is the best practical choice;
2. **deployment** checks: a 90-day longevity loop (no drift) and the compute case
   for scaling to a large district.

## The comparison (`visualizations/4_comparison/run_comparison.m` -> `viz_compare.m`)

Four controllers run the SAME plant, scenario, warmup, horizon and actuator limits,
so the only thing that differs is the controller (a fair comparison, not a strawman):

  - **KPC**            : direct multi-step Koopman prediction + one convex QP (this work)
  - **Koopman-LMPC**   : same Koopman lift, but an iterated one-step model + one QP
  - **Jacobian-LMPC**  : re-linearise the nonlinear plant every step + one QP
  - **NMPC**           : full nonlinear plant model + fmincon to convergence

Scenario: a severe stressed day (`mdot_scale = 0.35`, harder than the controller part),
in the capacity-binding regime, over a full day. Fairness: identical warmup
and demand forecast, the SAME demand-adequacy flow floor for all four, the full physical
actuator box (supply temperature and the network flows) for the baselines, a comparable,
large-enough finite-difference step for the gradient-based baselines, and a generous NMPC budget;
every controller's per-step solve and convergence is recorded. Worst-consumer demand-met
is a daily metric, so the controller is free to lead and recover around a peak; that is
what rewards the full-horizon planning the controllers differ on.

Result (live numbers, severe stressed day at mdot = 0.35):

| controller      | worst-cons met % | median solve | steps converged |
|-----------------|------------------|--------------|-----------------|
| **KPC**         | **97.6** (best)  | **~1.0 s**   | **96 / 96**     |
| NMPC            | 95.5             | ~60 s        | 86 / 96         |
| Koopman-LMPC    | 91.5             | ~0.6 s       | 96 / 96         |
| Jacobian-LMPC   | 87.4             | ~15 s        | 96 / 96         |

Under this severe bind KPC is the only controller that is BOTH comfortable and practical.
It keeps the worst consumer at 97.6% (the only one above the 95% comfort bar and the best
of all four), solves one convex QP in ~1 s, and stays feasible on every step. The
exact-model NMPC is close on comfort (95.5%) but pays heavily: it re-integrates the
nonlinear ODE inside every solve, so it is ~60x slower (~60 s/step, up to ~250 s) and, on
this hard problem, fails to converge in 10 of 96 steps. The iterated Koopman-LMPC (91.5%)
and the Jacobian-LMPC (87.4%) fall below the comfort bar; the gap to the iterated model is
the direct multi-step prediction. So KPC beats even the idealised exact-model NMPC here,
with no plant model and a convex QP.

## Stress sweep (`run_stress_sweep.m` -> `viz_stress_sweep.m`)

To show why the head-to-head is run at a hard point, the two fast Koopman controllers (KPC
direct multi-step vs the iterated one-step) are run across a range of flow capacities. They
track together at easy capacity and fan apart as it gets harder: KPC stays above the 95%
comfort bar down to about half-ish capacity while the iterated model drops below it sooner.
The direct multi-step prediction is what keeps KPC ahead, and the gap grows the more the
network is squeezed.

## Deployment (`run_longevity.m` -> `demo_robustness.m`)

  - **Longevity**: one continuous 90-day closed loop at the stressed half-capacity
    (`mdot_scale = 0.50`, ~8640 control steps). The per-day worst-consumer met % stays
    flat (drift +0.000 pp) and the per-day solve time stays flat (no compute creep), so
    nothing in the Koopman prediction or the QP degrades over a long run. A negligible
    0.012 % of steps the QP does not converge under the bind; the controller then holds
    its previous command and recovers, so comfort is unaffected.
  - **Scalability**: the KPC convex QP solves in about a second, roughly 0.1% of the
    `Ts = 900 s` sample budget, at the production size; its decision dimension grows
    linearly with the number of consumers, so there is large headroom to scale to a
    bigger district. The plant-rolling baselines re-integrate the full network ODE
    inside every solve (15-60 s at 5 consumers, and they grow with the network), so they
    do not scale.

## Checks (`tests/validate_robustness.m`)

| #  | invariant                                              | result                |
|----|--------------------------------------------------------|-----------------------|
| C1 | KPC QP feasible every step; baseline rates reported   | KPC 96/96; NMPC 86/96 |
| C2 | KPC worst-consumer met % >= 95 % bar AND best of all  | 97.6 % (next 95.5)    |
| C3 | KPC at least 10x faster than NMPC                      | ~60x                  |
| C4 | KPC median solve < 5 % of the Ts budget (real-time)   | 0.11 %                |
| L1 | 90-day worst-consumer drift <= 0.5 pp                 | +0.000 pp             |
| L2 | 90-day infeasible-step rate <= 0.1 %                  | 0.012 %               |

## Files

  - `../visualizations/4_comparison/run_comparison.m`  runs the four controllers, saves `comparison.mat`
  - `../visualizations/4_comparison/viz_compare.m`     renders the comparison figure
  - `../visualizations/4_comparison/baseline_cl_run.m` the shared fair harness for the plant-based baselines
  - `../visualizations/4_comparison/run_stress_sweep.m` runs the capacity sweep, saves `stress_sweep.mat`
  - `../visualizations/4_comparison/viz_stress_sweep.m` renders the stress-sweep figure
  - `../controller_closed_loop/lmpc/`                  the iterated Koopman-LMPC controller (build / solve / step loop)
  - `run_longevity.m`                                  the 90-day longevity run
  - `demo_robustness.m`                                the deployment summary + longevity figure
  - `tests/validate_robustness.m`                      the six acceptance invariants

The locked tune (`best_tune.mat`, `best_horizons.mat`, `best_alpha.mat`) lives in
`controller_closed_loop/results/` and the predictor fits (`vseq_fits_full.mat`) in
`predictor_open_loop/results/`; this part cross-references them rather than duplicating.
