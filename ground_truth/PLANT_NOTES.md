# AROMA plant: validation reference

## Scenarios

| name           | flow                  | T^(0,s)              | start | met % |
|---|---|---|---|---|
| nominal_midday | 1.0x                  | 20 C                 | 12:30 | 100  |
| morning_ramp   | 1.0x                  | 20 C                 | 06:00 | 100  |
| low_flow       | 0.5x                  | 20 C                 | 12:30 | 81.2 |
| high_flow      | 1.5x                  | 20 C                 | 12:30 | 100  |
| supply_step    | 1.0x                  | 20 -> 24 -> 18 C, 8 h phases | 12:30 | 90.1 |
| variable_flow  | 0.5 -> 1.5 -> 0.5x ramp | 20 C                | 12:30 | 81.5 |

90 min warmup before each. Same demand profile across all six.
low_flow / supply_step / variable_flow partial on purpose: substation
clips at c_max and pins T_r at T_r_min, leaves demand unmet, never
breaks the floor.

## Checks (validate_plant.m)

| #  | Invariant                                                                              | Result                          |
|----|----------------------------------------------------------------------------------------|---------------------------------|
| 1  | c^(i) <= d^(i)                                                                          | exact, every step               |
| 2  | T_r^(i) >= T_r_min                                                                      | exact, every step               |
| 3  | trunk Kirchhoff: q_F0->F1 = sum_i q^(i)                                                 | < 1e-5 kg/s                     |
| 4  | capacity-binding regime per scenario                                                    | each scenario                   |
| 5  | c = c_p (T_s - T_r) q                                                                   | < 4e-12 W                       |
| 6  | E_extracted / E_source                                                                  | 0.79 - 0.97 (within [0.5, 1])   |
| 7  | c^(i) >= 0 AND T_r^(i) <= T_s^(i)                                                       | exact, every step               |
| 8  | per-junction Kirchhoff (full incidence)                                                 | < 1e-12 kg/s                    |
| 9  | source-side return mixing: sum_i q^(i) T_r^(i) ~= Q_net T_0r                            | < 5 % residual                  |

## Spec equations -> code

| Equation                                                                | File                                                                                                  |
|-------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------|
| q_{k+1} = a q_k + (1-a) r_q,k,  a = exp(-Ts / eps)                       | shared/plant/flow_dynamics.m   ;  shared/plant/simulate_plant.m line 170                              |
| c = c_p (T_s - T_r) q   (mass-flow form)                                 | shared/network/physics/components/prosumer_house.m lines 21-26                                        |
| substation contracts: c <= d, T_r >= T_r_min, c >= 0, T_r <= T_s         | same file, by construction of clipping                                                                 |
| A_s q^(s,e) = [Q_net; -q^(i)],  A_r q^(r,e) = [-Q_net; q^(i)]            | shared/plant/simulate_plant.m (flow propagation),  shared/network/build_incidence.m (matrices)        |
| sum_in T_in q_in = T_node sum_out q_out  per junction                    | shared/network/physics/rhs_network.m lines 92-96                                                       |
| 1D pipe transport with first-order ambient loss                          | shared/network/physics/rhs_network.m  (f_DHC ODE called by ode15s in simulate_plant.m)                 |
| source CSTR: T_int = (1 - theta_1) T_in + theta_1 T_out, theta_1 = 0.7   | shared/network/physics/rhs_network.m lines 92-96 ;  shared/network/physics/measure_node_outlets.m 67-69 |

## Substation floor: how it stays >= T_r_min under demand spikes

prosumer_house.m:
```
c_max = c_p * q * (T_s - T_r_min)        % ceiling
c     = min(d, c_max)                     % clip
T_r   = T_s - c / (c_p * q)               % derived
```

If d <= c_max:  c = d,  T_r above floor.
If d >  c_max:  c = c_max,  T_r = T_r_min exactly,  shortfall in unmet column.

C2 example: T_s = 20, q = 1 kg/s, demand 30 kW.
c_max = 4180 * 1 * (20 - 15) = 20.9 kW.
c = 20.9,  T_r = 20 - 20.9 / 4180 = 15.0 C exactly.
Unmet: 9.1 kW.
