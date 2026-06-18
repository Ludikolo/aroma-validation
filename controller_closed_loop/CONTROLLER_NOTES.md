# KPC controller: closed-loop validation reference

## What this part does

```
min_{U, S_hi, S_lo, S_Ts, S_mix}
    -c_p * 1^T (F z_k + G U + G_d d)
    + alpha * (T_0s - T_0r) * Q_net
    + (Delta u)^T R_du (Delta u)
    + rho_hi ||S_hi||_1 + rho_lo ||S_lo||_1 + rho_Ts ||S_Ts||_1 + rho_mix ||S_mix||_1
s.t.
    c_p (F z_k + G U + G_d d) <= d + S_hi   (demand contract, soft)
    c_p (F z_k + G U + G_d d) >= -S_lo      (non-negative heat, soft)
    T_r <= T_s_pred + S_Ts                  (return below supply, soft)
    |g_N(z, U)| <= S_mix                    (mixing residual, soft)
    |Delta u_k| <= Delta u_max              (rate caps, hard)
    r_q^(i) >= min(s_adeq d^(i)/(c_p dT_floor), 0.95 r_q_max^(i))   (demand-adequacy floor, hard)
    M_s r_q = 0,  M_r r_q = 0               (per-node Kirchhoff on r_q, hard)
```

One sparse convex QP per Ts = 15 min. The state-prediction matrices
`F`, `G` come from the V_seq predictor in `predictor_open_loop/`.
Apply only the first block of the optimal `U^star`, then repeat at
the next sample (receding horizon).

The cost has four terms: a tracking term that maximises delivered
heat (linear in `U` via the V_seq predictor); an energy-penalty
term `alpha * (T_0s - T_0r) * Q_net` that pushes the source supply
temperature down when load allows (linearised around the current
operating point); a quadratic move-suppression penalty
`(Delta u)^T R_du (Delta u)` that smooths the actuator moves (this is
the only quadratic term, so the QP is strictly convex); and four soft
slack penalties that let the QP breathe under temporary infeasibility
without rejecting the step. The slack prices are not equal:
`rho_hi = rho_slack` on the demand-cap slack, while the non-negativity
and return-below-supply slacks cost ten times more
(`rho_lo = rho_Ts = 10 * rho_slack`, since a breach there is
unphysical), and `rho_mix` is set per tune.

Two specifics of the realised QP. The predictor folds the demand forecast in as
a known input, `c = c_p (F z_k + G U + G_d d)` (demand-in-V); since `d` is known
at solve time, `G_d d` is a constant offset and
the QP stays convex in `U`. And two of the constraints above are hard and always
on: the demand-adequacy flow floor, with `dT_floor = T_in_nom - T_r_min` and
`s_adeq = adequacy_safety` (locked 1.06, capped at `0.95 r_q_max`), which keeps
the planned flow able to physically carry the demand; and the Kirchhoff balance,
realised as the full per-node incidence `M_s r_q = 0`, `M_r r_q = 0` (the
rank-reduced supply / return incidence from `build_incidence_v2.m`, 13 + 13 = 26
node-balance rows per horizon step), not a single scalar sum.

## Locked tune

| Symbol             | Value  | Selected by                            |
|--------------------|--------|----------------------------------------|
| `Np`               | 12     | horizon sweep (robustness plateaus; nominal stays 100 %) |
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
| classical LMPC  | re-linearise the plant each step (Jacobian) + one QP               |
| NMPC            | nonlinear plant model, `fmincon` to convergence                    |
| **KPC (this work)** | V_seq direct multi-step + convex QP, locked tune above      |

## Checks (validate_controller.m)

| #  | Invariant                                                                   | Result on the demo set         |
|----|------------------------------------------------------------------------------|--------------------------------|
| 1  | Every QP returns `exitflag > 0` (feasible / optimal)                          | 96/96 steps                    |
| 2  | KPC met % = 100 on every (consumer, scenario) cell at design capacity         | 30/30 cells                    |
| 3  | KPC suite-mean met % at stressed (mdot = 0.6) >= 95 % bar                     | 100.000 %                      |
| 4  | KPC worst-consumer met % at stressed (100.000) > hold worst (82.206)          | +17.79 pp                      |
| 5  | Substation contract `c <= d` on the realised plant trace                      | exactly                        |
| 6  | Non-negative delivered heat `c >= 0` on the realised plant trace              | exactly                        |
| 7  | Mean QP solve time <= 5 % of `Ts` budget                                      | 851 ms mean, 1205 ms max       |

## Spec equations -> code

| Equation                                                                                       | File                                                                |
|------------------------------------------------------------------------------------------------|---------------------------------------------------------------------|
| Sparse convex QP assembly (cost, slacks, hard + soft constraints)                              | `controller_closed_loop/kpc/kpc_v2_solve.m`                          |
| Per-step receding-horizon outer loop (build z, call solve, apply, simulate, advance)            | `controller_closed_loop/kpc/kpc_step_loop.m`                         |
| Predictor matrix preparation: F, G, G_d for the optimisation                                    | `controller_closed_loop/kpc/build_kpc_v2_matrices.m`                 |
| Horizon truncation of the predictor block (Np_best from sweep)                                  | `controller_closed_loop/kpc/truncate_pred_to_Np.m`                   |
| Supply / return incidence used for the Kirchhoff equality on `r_q`                              | `controller_closed_loop/kpc/build_incidence_v2.m`                    |
| Six-scenario suite at design capacity                                             | `controller_closed_loop/results/scenario_suite.mat`                  |

## Headline numbers

At the locked tune on stressed capacity (`mdot_scale = 0.6`), the 24 h
demo has KPC at 100 % on every design cell and 100 %
worst-consumer met % at stressed (100 % suite-mean), against 90.7 %
suite-mean for hold-nominal, at a mean QP solve of about 851 ms
(0.09 % of the 15 min budget).

Part 4 (`visualizations/4_comparison/` plus
`robustness_scalability_sustainability/`) ships a head-to-head comparison of
KPC against three MPC baselines on a severe stressed day (`mdot_scale = 0.35`,
in the capacity-binding regime that pushes the controllers hardest): the iterated Koopman-LMPC (same
lift, a one-step model), the Jacobian-linearised classical LMPC, and a full
nonlinear NMPC. On a fair comparison (same plant, warmup, forecast, horizon and
demand-adequacy floor) KPC is the only controller both comfortable and practical:
97.6 % worst-consumer met (the only one above the 95 % bar and the best of all),
~1 s per convex QP, feasible every step. The exact-model NMPC reaches 95.5 % but
is ~60x slower and fails to converge in 10 of 96 steps; the iterated Koopman-LMPC
(91.5 %) and the Jacobian-LMPC (87.4 %) fall below the comfort bar. KPC is
data-driven and convex; the plant-rolling baselines need the exact or linearised
plant. See `robustness_scalability_sustainability/ROBUSTNESS_NOTES.md`.

Slack-bound theorem: the QP's optimal slack `S_hi^*` upper-bounds the
actual demand overshoot up to a bounded one-step predictor error
(`epsilon_pred`); this post-hoc safety bound is derived and verified
separately.

## Data

The committed results in `controller_closed_loop/results/` are:

- `best_tune.mat`, `best_horizons.mat`, `best_alpha.mat`: the locked,
  validated controller configuration (prediction and control horizons,
  cost weights, slack penalty) used throughout this code.
- `scenario_suite.mat`: the six 4 h closed loops at design capacity
  (S1..S6); demo_controller reads the KPC and hold-nominal traces to
  compute the headline 30-cell met % at design.
- `demo.mat`: produced by `demo_controller.m`; one deterministic 24 h
  closed loop at stressed capacity used by `validate_controller.m`.

All `.mat` files are committed so the demo runs in ~30 s on a clean
clone; no compute step is required to reproduce the headline numbers.
