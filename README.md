# AROMA validation

Code behind my MSc thesis on direct multi-step Koopman predictive control for
fifth-generation district heating. It pairs the AROMA 5GDHC plant simulator
(heating mode, 23 nodes, 29 edges, 5 prosumers) with a five-part validation.
`CHANGELOG.md` lists what changed since the June package.

Requirements: MATLAB R2024b or later (validated on R2025b) with the Optimization
Toolbox (`quadprog`, `fmincon`). Runs are deterministic; the convex-controller
comfort and feasibility numbers reproduce exactly, the NMPC leg (`fmincon`) to
about half a percentage point.

## Layout

    aroma-build/
      README.md   CHANGELOG.md   run_all_tests.m   startup.m
      shared/                                 plant model, network, lift helpers
      ground_truth/                           plant checks; plant_validation.pdf
      predictor_open_loop/                    the 47-feature lift and fits; predictor_validation.pdf, dictionary_study.pdf
      controller_closed_loop/                 the KPC convex QP; controller_validation.pdf
      robustness_scalability_sustainability/  stress, sweep, longevity; robustness_validation.pdf
      visualizations/                         figures and the benchmark; fair_comparison/comparison_validation.pdf
      tests/                                  the five validators, validate_*.m

## Fifteen-minute route

From the repo root:

    run_all_tests     % the five validators, one PASS/FAIL line per component, ~1 min
    demo_pulse        % a +5 K source pulse travelling through the network, ~2 min
    demo_predictor    % evaluates the predictor on held-out data, ~1 min

`run_all_tests` reloads the shipped result files and asserts each component's
invariants; nothing is re-simulated.

## What is where

Each component has a folder (shown in the tree above), a short explanation PDF, a
demo that recomputes it, and a validator that checks it.

- **plant**: `ground_truth/`, `plant_validation.pdf`, demos `demo_plant` and
  `demo_pulse`, validator `validate_plant`.
- **predictor**: `predictor_open_loop/`, `predictor_validation.pdf` and
  `dictionary_study.pdf`, demo `demo_predictor`, validator `validate_predictor`.
- **controller**: `controller_closed_loop/`, `controller_validation.pdf`, demo
  `demo_controller`, validators `validate_controller` and `validate_qp`.
- **comparison**: `visualizations/fair_comparison/` holds the write-up
  `comparison_validation.pdf` and the benchmark, figure `viz_airtight`, validator
  `validate_robustness`. (`viz_compare`, in `visualizations/4_comparison/`, draws
  the wider-box variant of the same four-controller run.)
- **robustness**: `robustness_scalability_sustainability/`,
  `robustness_validation.pdf`, demo `demo_robustness`, validator
  `validate_robustness`.

One validator (`validate_robustness`) covers both the comparison and the
robustness rows.

## Running each component

    startup
    demo_plant ;       validate_plant
    demo_predictor ;   validate_predictor
    demo_controller ;  validate_controller
    demo_robustness ;  validate_robustness
    viz_airtight       % the four-controller benchmark figures

Each demo recomputes its result and plots a short summary; `demo_plant` is the
slow one (about half an hour), the rest run in seconds to minutes. `demo_robustness`
is the exception: it reloads `longevity.mat` and `comparison.mat` and plots them
rather than re-simulating, because those two runs take hours. Each validator
reloads the saved artefacts and asserts the invariants (no re-simulation). The
`*_validation.pdf` in each folder holds the formulas, discusses every figure, and
lists the checks, each with the coded assert threshold next to the observed value.
