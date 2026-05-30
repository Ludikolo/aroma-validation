function p = params()
% PARAMS  Plant parameters for the AROMA 5GDHC simulator.
%   Every script reads parameters through this function so there are
%   no magic numbers scattered across the codebase.
%
%   Network: AROMA benchmark network. 23 nodes, 29
%   edges, 5 prosumers. Heating mode at the upper end of the 5GDHC
%   band (Tin_nom = 20 C, T_ext = 5 C).

%% Operating point (5GDHC heating mode)
p.mode     = 'heating';
p.Text     = 5;        % Ambient temperature [C]
p.Tin_nom  = 20;       % Nominal source supply [C]
p.mdot_nom = 3.0;      % Nominal mass flow [kg/s]
p.cp       = 4180;     % Specific heat of water [J/(kg*K)]
p.rho      = 1000;     % Water density [kg/m^3]

%% Manipulated input dT_c (offset on the source driver temperature)
% Effective supply range: Tin_nom + u in [14, 26] C, inside the 5GDHC band.
p.u_min = -6.0;
p.u_max = +6.0;
p.u_nom =  0.0;

%% Exogenous flow scaling alpha_m
p.w_min = 0.5;
p.w_max = 1.5;
p.w_nom = 1.0;

%% Output node (used by single-output scripts; harmless for the demo)
p.y_node     = 'F1';
p.y_node_idx = [];

%% Source CSTR
% rho cp V * dT/dt = rho cp F_in (T_ret - T) + kQ (T_c - T).
%   V      small volume so the source response settles in seconds and
%          does not mask network transport dynamics.
%   kQ     large enough that the source tracks its setpoint closely.
%   theta1 outlet weighting in the AROMA mixing rule.
p.source.V      = 0.05;
p.source.kQ     = 40000;
p.source.theta1 = 0.70;

%% Pipe heat loss to ambient
% Small coefficient (insulated 5GDHC pipes; Buffa et al. 2019).
p.pipe.alpha = 1e-6;

%% Prosumer (substation)
% dT_ref = 5 K nominal extraction.
%   Q_total = mdot_nom * cp * dT_ref = 62.7 kW network capacity.
% T_r_min = 15 C: the substation pins the return at the floor when
% capacity binds (rather than violating it). Yields ~21 kW c_max per
% consumer at design flow, right at the commercial-peak demand so
% the capacity envelope is genuinely meaningful.
p.consumer.dT_ref  = 5;
p.consumer.T_r_min = 15;

%% Flow mode
% 'exogenous'     : w(t) provided by the simulation script (default).
% 'demand_driven' : w(t) = min(w_max, sum demand / (cp * mdot_nom * dT_L)).
p.flow_mode = 'exogenous';

% Demand-driven configuration (ignored in exogenous mode).
% Per-consumer demand profile shape used by build_plant.
p.demand.profile_map = {'residential','commercial','shop','commercial','residential'};
p.demand.delta_T_L   = 5;
p.demand.w_max_flow  = 1.5;

%% Sample period
% Discrete cadence at which the trajectory is read out and the
% first-order flow update q_{k+1} = a q_k + (1-a) r_q,k runs. The
% pipe transport ODE itself integrates in continuous time (ode15s)
% and is independent of this choice.
p.Ts = 900;

%% Flow dynamics
% q_{k+1} = a * q_k + (1 - a) * r_q,k  with  a = exp(-Ts / eps).
% Time constants from the transport-delay scale (van der Heijde 2017):
% ring edges and stubs ~600 s, user-side stubs ~300 s.
p.flow_dyn.epsilon_supply = 600;
p.flow_dyn.epsilon_return = 600;
p.flow_dyn.epsilon_user   = 300;

%% Predictor knobs (used by predictor_open_loop/)
% Open-loop prediction horizon. The lift is fitted at every h in 1..H_max
% so the predictor RMSE-vs-horizon curve is directly readable.
p.o1.H_max = 16;

% Dictionary configuration.
%   delay_lags         lag indices for the source-temperature embedding
%                      (Takens 1981): keep short at Ts = 900 since one
%                      pipe transit fits in 1-3 samples.
%   use_delayed_bilinear adds T_F0(k-d) * q^(i) features (best v1 block).
%   use_substation_extras gates the saturation / regime / diurnal trig
%                      block. Off by default: it improves h = 1 R^2 but
%                      worsens long-horizon V_seq under the same data.
p.o2.use_delays            = true;
p.o2.use_flow_disc         = false;
p.o2.use_delayed_bilinear  = true;
p.o2.delay_lags            = [1, 2];
p.o2.use_substation_extras = false;
% Exergy block: 12 bilinears q*(T_s-Tamb), q*(T_r-Tamb), Qnet*(T0s-Tamb),
% Qnet*(T0r-Tamb) appended to the base lift, giving n_z = 49. This is the
% production lift; the controller tune below is fitted on it.
p.o2.use_exergy            = true;

%% Data set (PRBS-excited trajectories used to fit the predictor)
% 20 train + 2 val + 3 test trajectories at 24 h each. The held-out test
% set is the one the validation invariants read.
p.data.n_train    = 20;
p.data.n_val      = 2;
p.data.n_test     = 3;
p.data.traj_dur_h = 24;
p.data.traj_dur_s = 24 * 3600;
p.data.rng_seed   = 42;

%% Excitation (PRBS on r_q + dT_c)
% Dwell times are pairwise coprime ([2 3 5 7 11] * Ts at Ts = 900) so
% cycled PRBS across edges does not phase-lock at long horizons
% (Soderstrom and Stoica 1989, sec. 5; Ljung 1999, ch. 13).
p.excite.r_q_lo_factor = 0.3;
p.excite.r_q_hi_factor = 1.5;
p.excite.r_q_dwells_s  = [1800 2700 4500 6300 9900];
p.excite.T0s_amp       = 6;
p.excite.T0s_dwell_s   = 3600;

end
