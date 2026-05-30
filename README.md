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
- `robustness_scalability_sustainability/` runs the deployment checks:
  robustness sweeps, compute scaling, and the 5GDHC sustainability KPIs.

Run everything from the repo root:

```matlab
startup
demo_plant ;       validate_plant
demo_predictor ;   validate_predictor
demo_controller ;  validate_controller
demo_robustness ;  validate_robustness
```

Each `demo_*` recomputes its result, saves `results/demo.mat`, and plots
a short summary; each `validate_*` reloads that file and asserts the
invariants (it does not re-simulate). The plant demo is the slow one (a
few minutes); the rest run in seconds. The `*_NOTES.md` in each folder
hold the formulas and the list of checks.

The production lift is the 49-feature exergy dictionary and the locked
controller tune is `(Np, Nc, alpha, rho_slack) = (5, 1, 1, 1.195)` with
an adequacy-safety margin of `1.06`. On top of this controller the
thesis also develops a few extensions (a practical ISS certificate, an
online RLS predictor update, an ADMM split for larger networks, and a
disturbance-observer offset filter); those are covered in the thesis
rather than exercised here.
