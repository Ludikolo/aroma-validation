# Visualizations

Short scripts, grouped by part, that each make one figure showing that a part of
the system works: the plant physics, the predictor, the controller, and how the
controller compares against standard MPC baselines.

The part 1-3 scripts run live: they load only the trained models (the predictor
maps and the locked controller tune) and compute every trajectory, forecast and
metric on the spot. The part 4 comparison is the exception, because the baselines
roll the nonlinear plant and take minutes to hours; its expensive runs are saved
and the viz scripts render from them in seconds.

## How to run

From the repo root in MATLAB:

```matlab
startup            % puts all the code on the path
run('visualizations/1_plant/viz_plant_step.m')
```

Each script clears the workspace, adds the path itself, makes its figure, saves a
vector PDF next to itself, and prints a short summary to the console. Run them in
any order.

## 1. Plant (`1_plant/`)

- `viz_plant_step.m` steps the source supply temperature up by 6 C and shows the
  warm front reaching each consumer after its transport delay. The point is the
  1D pipe-transport physics: nearer consumers respond first, none exceed the
  source.
- `viz_plant_operation.m` runs a full day at nominal flow and checks two
  conservation laws on the trajectory (substation energy balance and Kirchhoff
  mass conservation). Both residuals come out at machine precision.

## 2. Predictor (`2_predictor/`)

- `viz_predictor_future.m` makes a 4 h-ahead forecast at every step over a fresh
  normal operating day and lays it on top of what the plant actually does. The
  direct multi-step Koopman map tracks the truth; the iterated one-step map falls
  behind.
- `viz_predictor_error.m` quantifies this: forecast error vs horizon over a normal
  day and all five consumers, for the direct multi-step map, the iterated one-step
  map, and a zero-order-hold null. The direct map is lowest at every horizon and the
  gap to the iterated map grows with the horizon.

## 3. Controller (`3_controller/`)

Both run the KPC controller in closed loop for 24 h at nominal flow (normal
operation; the stressed comparison is part 4).

- `viz_controller_tracking.m` shows, for every consumer, the delivered heat
  (solid) under the demand (dotted) over the whole day, so you can see the
  controller meet each consumer's demand through the diurnal cycle (100 % met).
- `viz_controller_signals.m` shows the control commands (source supply
  temperature and per-consumer flows), the resulting supply temperatures, and the
  solve time per step. It also prints a set of constraint checks (actuator
  rate limits, return-below-supply, QP feasibility) that all pass.

## 4. Comparison (`4_comparison/`)

- `run_comparison.m` runs all four controllers (KPC, the iterated Koopman-LMPC, a
  Jacobian-LMPC, and a full nonlinear NMPC) on the same severe stressed day through
  one shared harness, over a full 24 h, so the only difference is the controller. It
  saves `comparison.mat`. `viz_compare.m` renders it: how each controller tracks the
  demand for the hardest consumer, plus the worst-consumer demand met and the median
  solve time per step.
- `run_stress_sweep.m` runs KPC and the iterated Koopman-LMPC across a range of flow
  capacities and saves `stress_sweep.mat`; `viz_stress_sweep.m` renders the gap
  between them growing as the network is squeezed harder.
- `baseline_cl_run.m` is the shared closed-loop harness for the two plant-based
  baselines; it mirrors the KPC closed loop exactly so the comparison is fair.

The baselines roll the nonlinear plant inside their optimiser, so they are slow
(`run_comparison` takes minutes to hours, dominated by NMPC). That cost gap is itself
part of the result, which is why the saved `.mat` ship and the viz scripts render
from them in seconds.
