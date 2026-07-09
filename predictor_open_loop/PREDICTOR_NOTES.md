# Koopman predictor: validation reference

## What this part does

```
c_k^(i) = c_p * theta_k^(i)
theta_{k+h}^(i) = c_{h,i}' z_k + d_{h,i}' V_k^(h)
V_k^(h) = [u_k; ...; u_{k+h-1};  d_{k+1}; ...; d_{k+h}]
```

The substation heat extraction `c^(i)` is exactly proportional to the
bilinear `theta^(i) = (T_s - T_r) * q`, so a predictor for `theta`
predicts `c`. The lifted state `z_k = psi(xi_k)` is a fixed vector of
47 functions of the augmented plant state (the base 36-feature lift
plus a 11-feature exergy block). `V_k^(h)` stacks the applied control
inputs `u = [T^(0,s); r_q; T^(i,r)]` and the demand forecast `d` over
the prediction horizon h. The demand is measured/forecast at each
substation, so it enters as a known input rather than something to
predict.

One linear fit per (consumer, horizon) pair gives V_seq the direct
multi-step structure. No recursion, no one-step error accumulation.

## Dictionary (47 features)

The base lift is 36 features in five blocks; the production lift adds a
11-feature exergy block on top, giving `n_z = 47`.

| block | symbol             | count | what it captures                                        |
|-------|--------------------|-------|---------------------------------------------------------|
| A     | thermal            | 12    | T_{F0}, T^(i,s), T^(i,r), T^(0,r)                       |
| B     | hydraulic          | 6     | q^(i), Q_net                                            |
| C     | **theta**          | 5     | theta^(i) = (T^(i,s) - T^(i,r)) * q^(i)                 |
| const | const              | 1     | intercept                                               |
| E     | source lag         | 2     | T_{F0}(k-1), T_{F0}(k-2)                                |
| F     | delayed bilinear   | 10    | T_{F0}(k-d) * q^(i),  d in {1, 2}                       |
| G     | **exergy**         | 11    | q^(i)*(T^(i,s)-T_amb), q^(i)*(T^(i,r)-T_amb), Q_net*(T^(0,r)-T_amb) |

Dropping the theta block from the base lift costs about 0.21 in R^2 at
h = 1 (over 20 x more than dropping any other single block): the
substation contract `c = c_p * theta` is exact, so theta in the lift
makes the heat-extraction cost linear in z. The exergy block adds
energy-transport-capacity features that help most when flow is scarce;
they make theta partly redundant for open-loop fit but improve the
closed-loop delivery the controller cares about.

## Predictors compared

| name              | structure                                                              |
|-------------------|------------------------------------------------------------------------|
| ZOH               | theta_{k+h}^(i) = theta_k^(i)                                           |
| iterated (A, B)   | z_{k+1} = A z_k + B [u_k; d_{k+1}], applied h times, theta read from z  |
| **V_seq (this work)** | theta_{k+h}^(i) = c_{h,i}' z_k + d_{h,i}' V_k^(h), per (h, i)       |

The iterated baseline gets the same lift, the same demand input, and the
same operational weighting as V_seq, and its ridge lambda is picked among
the stable fits (rho(A) < 1) so the rollout never blows up. That keeps it
an honest opponent: the only thing that differs is direct-multi-step
versus rolled-out one-step.

V_seq uses split-ridge regularisation: separate lambda for the
state-coefficient block c and the input-coefficient block d, picked
on the validation set per (h, i). Inputs scale 40 * h coefficients
(35 controls + 5 demand, per step) versus 47 state coefficients, so one
lambda is the wrong knob.

## Checks (validate_predictor.m)

| #  | Invariant                                                              | Result on the demo set    |
|----|-------------------------------------------------------------------------|---------------------------|
| 1  | V_seq RMSE <= ZOH RMSE at every horizon                                 | true, h = 1..16           |
| 2  | V_seq RMSE < iterated A,B RMSE at h = 16                                | 12.1 % lower, suite mean  |
| 3  | V_seq 1-step suite R^2 >= 0.7                                           | 0.712                     |
| 4  | Worst-consumer NRMSE at h = 1 <= 30 % of peak demand                    | 19.2 %                    |
| 5  | Theta-block contribution to R^2 at h = 1 >= 0.10                        | 0.207                     |
| 6  | q_{k+1} = a q_k + (1-a) r_q,k on every test trajectory                  | machine precision         |

These run on the held-out PRBS test set. The direct-multi-step gain over
the iterated baseline is modest there (12 % at h = 16) but large on a
normal operating day, the regime the controller runs in: there the
iterated rollout compounds and V_seq is about 3x lower RMSE from h = 5
onward (see the operational tracking and parity figures). At h = 1 the
two are equal, as expected for a one-step forecast.

## Spec equations -> code

| Equation                                                                       | File                                                       |
|--------------------------------------------------------------------------------|------------------------------------------------------------|
| psi(xi) = 47 features (base 36 + 11 exergy bilinears)                          | predictor_open_loop/koopman/candidate_library.m            |
| theta_{k+h}^(i) = c_{h,i}' z_k + d_{h,i}' V_k^(h), per (h, i), split-ridge     | predictor_open_loop/koopman/fit_vseq.m                     |
| z_{k+1} = A z_k + B [u_k; d_{k+1}] (demand-aware, stable iterated baseline)    | predictor_open_loop/koopman/fit_iterated.m                 |
| Per-junction mixing residual g_N(z) used as a soft constraint downstream       | predictor_open_loop/koopman/build_mixing.m                 |
| q_{k+1} = a q_k + (1 - a) r_q,k,  a = exp(-Ts / eps)                            | shared/plant/flow_dynamics.m                               |

## Data

PRBS-excited trajectories at Ts = 900 s, 24 h each, on the AROMA
plant. 20 train + 2 val + 3 test. The 3 test trajectories ship with
the repository under `predictor_open_loop/results/test_*.mat`. The
full 25-file set can be regenerated with
`predictor_open_loop/data/generate_data.m`.

The PRBS set gives broadband excitation for the dynamics but spends a lot
of its time capacity-limited (`c = c_max != d`), which dilutes the demand
coefficient. To calibrate the predictor for the regime it is deployed in,
the fit also uses a handful of smooth operating days
(`predictor_open_loop/results/op_*.mat`, made by
`predictor_open_loop/data/generate_operational.m`) and weights those
samples more (`p.o1.op_weight`, a sqrt weight in least squares). This is
what puts the multi-step parity on the diagonal on a normal day; the PRBS
test R^2 (T3) is the guard that it does not over-specialise.

PRBS dwell times use pairwise-coprime prime multipliers
(`[2, 3, 5, 7, 11] * Ts`) so the cycled excitation across edges does
not phase-lock at long horizons. Coprime dwell spacing is a practitioner
heuristic for decorrelating multi-channel PRBS; for PRBS / input
experiment design in general see Ljung 1999 ch. 14.
