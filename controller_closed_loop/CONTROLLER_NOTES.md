# KPC controller: closed-loop validation reference

## What this part does

```
min_{U, S_hi, S_lo, S_Ts, S_mix}
    -c_p * 1^T (F z_k + G U)
    + alpha * (T_0s - T_0r) * Q_net
    + rho_slack * ||(S_hi, S_lo, S_Ts, S_mix)||_1
s.t.
    c_p (F z_k + G U) <= d + S_hi          (demand contract, soft)
    c_p (F z_k + G U) >= -S_lo             (non-negative heat, soft)
    T_r <= T_s_pred + S_Ts                 (return below supply, soft)
    |g_N(z, U)| <= S_mix                   (mixing residual, soft)
    |Delta u_k| <= Delta u_max             (rate caps, hard)
    1^T r_q_supply = 1^T r_q_return        (Kirchhoff on r_q, hard)
```

One sparse convex QP per Ts = 15 min. The state-prediction matrices
`F`, `G` come from the V_seq predictor in `predictor_open_loop/`.
Apply only the first block of the optimal `U^star`, then repeat at
the next sample (receding horizon).

The cost has three terms: a tracking term that maximises delivered
heat (linear in `U` via the V_seq predictor), an energy-penalty
term `alpha * (T_0s - T_0r) * Q_net` that pushes the source supply
temperature down when load allows (linearised around the current
operating point), and three soft slack penalties that let the QP
breathe under temporary infeasibility without rejecting the step.

## Locked tune

| Symbol             | Value  | Selected by                            |
|--------------------|--------|----------------------------------------|
| `Np`               | 5      | cross-scenario sweep (S1..S6)          |
| `Nc`               | 1      | cross-scenario sweep (Occam pick)      |
| `alpha`            | 1      | two-pass Pareto over alpha             |
| `rho_slack`        | 1.195  | Bayesian optimisation                  |
| `adequacy_safety`  | 1.06   | demand-adequacy flow floor             |
| `dT_0s_max`        | 7.5    | rate cap (verified non-binding)        |
| `dr_q_max`         | 4.5    | rate cap                               |
| `dT_ir_max`        | 15.0   | rate cap                               |
| QP algorithm       | `interior-point-convex` | sensitivity sweep confirms active-set p99 is 4x worse |

Saved tune artefacts:
- `controller_closed_loop/results/best_horizons.mat` -> `Np_best`, `Nc_best`
- `controller_closed_loop/results/best_alpha.mat`    -> `alpha_star`
- `controller_closed_loop/results/best_tune.mat`     -> `best_tune` struct

## Baselines compared

| name            | structure                                                          |
|-----------------|--------------------------------------------------------------------|
| hold-nominal    | apply `r_q = mdotEdges`, `T_0s = T_in_nom`; no controller          |
| demand-following RBC | per-edge `r_q^(i)` proportional to demand `d^(i)`             |
| Koopman-LMPC    | same predictor lift, single one-step `(A, B)` instead of V_seq     |
| classical LMPC  | Jacobian-linearised plant, 765-state                               |
| NMPC multi-start| nonlinear plant model, `fmincon` with 3 random restarts            |
| **KPC v2 (this work)** | V_seq direct multi-step + convex QP, locked tune above      |

## Checks (validate_controller.m)

| #  | Invariant                                                                   | Result on the demo set         |
|----|------------------------------------------------------------------------------|--------------------------------|
| 1  | Every QP returns `exitflag > 0` (feasible / optimal)                          | 96/96 steps                    |
| 2  | KPC met % = 100 on every (consumer, scenario) cell at design capacity         | 30/30 cells                    |
| 3  | KPC suite-mean met % at stressed (mdot = 0.6) >= 95 % bar                     | 100 %                          |
| 4  | KPC worst-consumer met % at stressed > hold worst-consumer met %               | +17.8 pp                       |
| 5  | Substation contract `c <= d` on the realised plant trace                      | exactly                        |
| 6  | Non-negative delivered heat `c >= 0` on the realised plant trace              | exactly                        |
| 7  | Mean QP solve time <= 5 % of `Ts` budget                                      | ~95 ms / 900 s budget          |

## Spec equations -> code

| Equation                                                                                       | File                                                                |
|------------------------------------------------------------------------------------------------|---------------------------------------------------------------------|
| Sparse convex QP assembly (cost, slacks, hard + soft constraints)                              | `controller_closed_loop/kpc/kpc_v2_solve.m`                          |
| Per-step receding-horizon outer loop (build z, call solve, apply, simulate, advance)            | `controller_closed_loop/kpc/kpc_step_loop.m`                         |
| Predictor matrix preparation: F, G, G_d for the optimisation                                    | `controller_closed_loop/kpc/build_kpc_v2_matrices.m`                 |
| Horizon truncation of the predictor block (Np_best from sweep)                                  | `controller_closed_loop/kpc/truncate_pred_to_Np.m`                   |
| Supply / return incidence used for the Kirchhoff equality on `r_q`                              | `controller_closed_loop/kpc/build_incidence_v2.m`                    |
| Six-scenario suite at design capacity (T8 artefact)                                             | `controller_closed_loop/results/scenario_suite.mat`                  |

## Headline numbers

At the locked tune on stressed capacity (`mdot_scale = 0.6`), the 24 h
demo has KPC v2 at 100 % suite-mean and 100 % worst-consumer met %,
against 90.7 % worst-consumer for hold-nominal, at a mean QP solve of
about 95 ms (0.01 % of the 15 min budget).

The broader thesis compares KPC against five alternatives (hold,
demand-following RBC, iterated Koopman-LMPC, classical Jacobian-
linearised LMPC, multi-start NMPC) under a fair per-controller retune.
Under that retune KPC has the highest suite-mean met % on both the
stressed and the out-of-distribution scenarios (margins +0.4 to +1.9 pp
stressed, +1.3 to +6.7 pp out-of-distribution). The iterated Koopman-LMPC
ties KPC at the 100 % comfort ceiling under the matched tune, so KPC's
edge there is its closed-loop prediction quality (about 34 % lower
in-loop prediction RMSE and far less constraint slack), not a met %
knockout. KPC also uses about a third less source energy than the
classical and reduced-order baselines and solves 7-19 x faster, which is
the multi-criterion case the thesis actually makes.

This repository ships only the KPC + hold-nominal pair used by the demo
and validation; the alternative baselines live with the broader thesis
codebase.

Slack-bound theorem: the QP's optimal slack `S_hi^*` upper-bounds the
actual demand overshoot up to a bounded one-step predictor error
(`epsilon_pred`); this post-hoc safety bound is proved and checked in
the thesis.

## Data

The committed results in `controller_closed_loop/results/` are:

- `best_tune.mat`, `best_horizons.mat`, `best_alpha.mat`: the locked
  tune coordinates from cross-scenario Pareto + BO.
- `scenario_suite.mat`: the six 4 h closed loops at design capacity
  (S1..S6); demo_controller reads the KPC and hold-nominal traces to
  compute the headline 30-cell met % at design.
- `demo.mat`: produced by `demo_controller.m`; one deterministic 24 h
  closed loop at stressed capacity used by `validate_controller.m`.

All `.mat` files are committed so the demo runs in ~30 s on a clean
clone; no compute step is required to reproduce the headline numbers.
