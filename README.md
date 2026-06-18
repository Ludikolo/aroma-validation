# aroma validation

This is the code behind my MSc thesis on direct multi-step Koopman
predictive control for fifth-generation district heating. It pairs the
AROMA 5GDHC plant simulator (heating mode, 23 nodes, 29 edges, 5
prosumers) with a four-part validation:

- `ground_truth/` checks the plant against its design operating points;
- `predictor_open_loop/` checks the Koopman lift and the direct
  multi-step V_seq predictor;
- `controller_closed_loop/` checks the convex-QP MPC controller (KPC)
  that runs on top of the predictor;
- `robustness_scalability_sustainability/` runs the deployment checks: a
  long-horizon longevity loop; the head-to-head comparison against the
  LMPC and NMPC baselines lives in `visualizations/4_comparison/`.

Run everything from the repo root:

```matlab
startup
demo_plant ;       validate_plant
demo_predictor ;   validate_predictor
demo_controller ;  validate_controller
demo_robustness ;  validate_robustness
viz_compare        % the four-controller head-to-head figure
```

Each `demo_*` recomputes its result, saves `results/demo.mat`, and plots
a short summary; each `validate_*` reloads that file and asserts the
invariants (it does not re-simulate). The plant demo is the slow one (a
few minutes); the rest run in seconds. The `*_NOTES.md` in each folder
hold the formulas and the list of checks.

The Part 4 figures are `visualizations/4_comparison/viz_compare.pdf` (the
four-controller comparison on a severe day) and `stress_sweep.pdf` (the
direct multi-step advantage growing as the network is squeezed), plus
`robustness_scalability_sustainability/figures/longevity.pdf` (the 90-day
deployment). Regenerating the comparison and the sweep is expensive
(`run_comparison`, `run_longevity`, `run_stress_sweep` take minutes to
hours); their saved `.mat` ship with the repo, so `viz_compare` and
`demo_robustness` render from them instantly.

The production lift is the 49-feature exergy dictionary and the locked
controller tune is `(Np, Nc, alpha, rho_slack) = (12, 1, 1, 1.195)` with
an adequacy-safety margin of `1.06`. On top of this controller the
thesis also develops a few extensions (a practical ISS certificate, an
online RLS predictor update, an ADMM split for larger networks, and a
disturbance-observer offset filter); those are covered in the thesis
rather than exercised here.
