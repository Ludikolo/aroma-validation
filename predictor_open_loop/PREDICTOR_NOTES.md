# Koopman predictor: validation reference

## What this part does

```
c_k^(i) = c_p * theta_k^(i)
theta_{k+h}^(i) = c_{h,i}' z_k + d_{h,i}' V_k^(h)
```

The substation heat extraction `c^(i)` is exactly proportional to the
bilinear `theta^(i) = (T_s - T_r) * q`, so a predictor for `theta`
predicts `c`. The lifted state `z_k = psi(xi_k)` is a fixed vector of
49 functions of the augmented plant state (the base 37-feature lift
plus a 12-feature exergy block); `V_k^(h)` is the stacked input +
demand sequence over the prediction horizon h.

One linear fit per (consumer, horizon) pair gives V_seq the direct
multi-step structure. No recursion, no one-step error accumulation.

## Dictionary (49 features)

The base lift is 37 features in five blocks; the production lift adds a
12-feature exergy block on top, giving `n_z = 49`.

| block | symbol             | count | what it captures                                        |
|-------|--------------------|-------|---------------------------------------------------------|
| A     | thermal            | 13    | T^(0,s), T_{F0}, T^(i,s), T^(i,r), T^(0,r)              |
| B     | hydraulic          | 6     | q^(i), Q_net                                            |
| C     | **theta**          | 5     | theta^(i) = (T^(i,s) - T^(i,r)) * q^(i)                 |
| const | const              | 1     | intercept                                               |
| E     | source lag         | 2     | T_{F0}(k-1), T_{F0}(k-2)                                |
| F     | delayed bilinear   | 10    | T_{F0}(k-d) * q^(i),  d in {1, 2}                       |
| G     | **exergy**         | 12    | q^(i)*(T^(i,s)-T_amb), q^(i)*(T^(i,r)-T_amb), Q_net*(T^(0,s)-T_amb), Q_net*(T^(0,r)-T_amb) |

Dropping the theta block from the base lift costs about 0.21 in R^2 at
h = 1 (over 20 x more than dropping any other single block): the
substation contract `c = c_p * theta` is exact, so theta in the lift
makes the heat-extraction cost linear in z. The exergy block adds
energy-transport-capacity features that help most when flow is scarce;
they make theta partly redundant for open-loop fit but improve the
closed-loop delivery the controller cares about.

## Predictors compared

| name              | structure                                                          |
|-------------------|--------------------------------------------------------------------|
| ZOH               | theta_{k+h}^(i) = theta_k^(i)                                       |
| iterated (A, B)   | z_{k+1} = A z_k + B u_k, applied h times, theta read from z         |
| **V_seq (this work)** | theta_{k+h}^(i) = c_{h,i}' z_k + d_{h,i}' V_k^(h), per (h, i)   |

V_seq uses split-ridge regularisation: separate lambda for the
state-coefficient block c and the input-coefficient block d, picked
on the validation set per (h, i). Inputs scale 35 * h coefficients
versus 49 state coefficients, so one lambda is the wrong knob.

## Checks (validate_predictor.m)

| #  | Invariant                                                              | Result on the demo set    |
|----|-------------------------------------------------------------------------|---------------------------|
| 1  | V_seq RMSE <= ZOH RMSE at every horizon                                 | true, h = 1..16           |
| 2  | V_seq RMSE < iterated A,B RMSE at h = 16                                | -13.1 % suite mean        |
| 3  | V_seq 1-step suite R^2 >= 0.7                                           | 0.733                     |
| 4  | Worst-consumer NRMSE at h = 1 <= 30 % of peak demand                    | 17.9 %                    |
| 5  | Theta-block contribution to R^2 at h = 1 >= 0.10                        | 0.207                     |
| 6  | q_{k+1} = a q_k + (1-a) r_q,k on every test trajectory                  | machine precision         |

## Spec equations -> code

| Equation                                                                       | File                                                       |
|--------------------------------------------------------------------------------|------------------------------------------------------------|
| psi(xi) = 49 features (base 37 + 12 exergy bilinears)                          | predictor_open_loop/koopman/candidate_library.m            |
| theta_{k+h}^(i) = c_{h,i}' z_k + d_{h,i}' V_k^(h), per (h, i), split-ridge     | predictor_open_loop/koopman/fit_vseq.m                     |
| z_{k+1} = A z_k + B u_k (iterated baseline)                                    | predictor_open_loop/koopman/fit_iterated.m                 |
| Per-junction mixing residual g_N(z) used as a soft constraint downstream       | predictor_open_loop/koopman/build_mixing.m                 |
| q_{k+1} = a q_k + (1 - a) r_q,k,  a = exp(-Ts / eps)                            | shared/plant/flow_dynamics.m                               |

## Data

PRBS-excited trajectories at Ts = 900 s, 24 h each, on the AROMA
plant. 20 train + 2 val + 3 test. The 3 test trajectories ship with
the repository under `predictor_open_loop/results/test_*.mat`. The
full 25-file set can be regenerated with
`predictor_open_loop/data/generate_data.m`.

PRBS dwell times use pairwise-coprime prime multipliers
(`[2, 3, 5, 7, 11] * Ts`) so the cycled excitation across edges does
not phase-lock at long horizons (Soderstrom and Stoica 1989 sec. 5;
Ljung 1999 ch. 13).
