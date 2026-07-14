% DEMO_PLANT  Open-loop demonstration and validation of the AROMA
%   5GDHC plant simulator. One file: builds the plant, runs six 24 h
%   scenarios that span the operating envelope, saves trajectories
%   plus per-scenario figures, and finally checks nine physical
%   invariants on the saved data.
%
%   See plant_validation.pdf, Section 3 (the figures) and Section 4 (the checks).
%
%   The plant simulator is the AROMA reference implementation; this
%   file is the demonstration and validation
%   layer on top of it; the plant code under shared/plant and
%   shared/network is not modified here.
%
%   Inputs to the plant (chosen per scenario):
%     - demand d^(i)(t) per consumer (compute_prosumer_Q)
%     - flow reference r_q^(*)(t) per edge / user
%     - source supply temperature T^(0,s)(t)
%     - initial network state (cold steady state at T_in_nom)
%
%   Outputs (read off the trajectory):
%     - per consumer T_s^(i), T_r^(i), c^(i)
%     - per edge q^(s,e), q^(r,e), q^(i) for users
%     - source return T^(0,r) and the full pipe-cell temperature field
%
%   The simulator integrates 1D pipe transport on dx = 10 m cells in
%   continuous time via ode45, runs first-order flow dynamics on
%   each edge, models the source as a CSTR (perfectly mixed control
%   volume), and lets the substation enforce c <= d and T_r >= T_r_min
%   by clipping.

clear; clc;
startup;
p = params();

% Run the demo at the fine sample rate (Ts = 60 s) regardless of what
% params.m carries on disk. The plant ODE itself runs in continuous
% time via ode45; p.Ts only sets the cadence at which q is updated and
% the trajectory is read out. Sixty samples per hour gives the figures
% enough resolution to show transient behaviour clearly.
p.Ts = 60;

fprintf('\n=== AROMA plant demo ===\n');
fprintf('Plant sampled at Ts = %g s with eps = (%g, %g, %g) s\n', ...
        p.Ts, p.flow_dyn.epsilon_supply, p.flow_dyn.epsilon_return, p.flow_dyn.epsilon_user);

fprintf('Building plant...\n');
[net, z0_cold] = build_plant(p);
ei = edge_user_index(net);
n_user = ei.n_user;

% Trunk supply edge F0 -> F1 (flow into the network). Used at the end
% to verify mass conservation: at steady state the trunk inflow must
% equal the sum of consumer flows.
F0_F1_edge = ei.Es(1);
R0_idx     = find(strcmp({net.Nodes.name}, 'R0'));

% Scenarios. Each scenario is a 24 h window run from a 90 min warmup
% so the pipe temperature field and the first-order flow state both
% reach their scenario-relative steady state before we record. The
% same diurnal demand profile is used across every scenario so that
% all observed differences come from the control inputs (flow factor
% and T_0s schedule), not from a different demand realisation.
T_warm = 90 * 60;
mdotE  = net.mdotEdges(:);

% A time-varying flow ramp for the variable_flow scenario is defined
% as a local function at the bottom of this script (variable_flow_profile).
% Holds at 0.5x for the first 4 h, ramps linearly to 1.5x by hour 8,
% holds at 1.5x for 8 h, ramps back to 0.5x by hour 20, holds. This
% exercises the simulator under a non-trivial dynamic flow input,
% the kind of dynamic flow input a downstream controller produces.

scenarios = struct( ...
    'name', {'nominal_midday', 'morning_ramp', 'low_flow', 'high_flow', 'supply_step', 'variable_flow'}, ...
    'desc', { ...
        'nominal flow + full 24 h diurnal demand starting 12:30', ...
        'nominal flow + full 24 h diurnal demand starting 06:00 (time-shifted check)', ...
        '0.5x nominal flow + diurnal demand (capacity-binding at peaks)', ...
        '1.5x nominal flow + diurnal demand (over-capacity)', ...
        'T_0s steps 20, 24, 18 in three 8 h phases, nominal flow, 12:30 start', ...
        'time-varying flow 0.5x to 1.5x ramp, T_0s = 20, 12:30 start (dynamic-flow check)' ...
    }, ...
    't_offset_h', {12.5, 6.0, 12.5, 12.5, 12.5, 12.5}, ...
    'T_sim_min', {1440, 1440, 1440, 1440, 1440, 1440}, ...
    'flow_factor', {1.0, 1.0, 0.5, 1.5, 1.0, NaN});  % NaN signals "use the function below"

% T_0s schedule per scenario. Default: hold at T_in_nom (= 20 C).
% Total length is 24 h for every scenario so each run covers a full
% diurnal demand cycle (residential morning peak around 07:30,
% commercial midday peak around 13:00, shop pm peak around 17:00,
% evening residential around 19:00, night, then back to morning).
% F0 -> F1 trunk transit alone is 500 m / 0.334 m/s about 25 min at
% nominal flow. The furthest consumer C5 is much slower than the trunk
% velocity suggests: its terminal pipes carry only 0.12-0.24 kg/s, so
% the first arrival is around 8.7 h (measured in demo_pulse.m), not the
% 94 min a trunk-velocity estimate over the 1883 m path would give.
% The supply_step scenario uses three 8 h phases so the main consumers
% have hours of settled response per phase after the wave arrives.
T_0s_funs = { ...
    @(t) p.Tin_nom, ...
    @(t) p.Tin_nom, ...
    @(t) p.Tin_nom, ...
    @(t) p.Tin_nom, ...
    @(t) p.Tin_nom + 4*(t >= 8*3600 & t < 16*3600) - 2*(t >= 16*3600), ...
    @(t) p.Tin_nom ...
};

%% Run scenarios
out = cell(1, numel(scenarios));
for k = 1:numel(scenarios)
    sc = scenarios(k);
    fprintf('\n[%d/%d] %s\n', k, numel(scenarios), sc.name);
    fprintf('  %s\n', sc.desc);

    t_offset_s = sc.t_offset_h * 3600;
    T_sim      = sc.T_sim_min * 60;

    % Choose the flow input. Numeric flow_factor means hold-constant;
    % NaN means use the variable_flow ramp defined above.
    if isnan(sc.flow_factor)
        r_q_fun = @(t) variable_flow_profile(t, mdotE);
    else
        r_q_const = sc.flow_factor * mdotE;
        r_q_fun = @(t) r_q_const;
    end

    % Warmup at scenario start with hold-nominal so each scenario sees
    % a representative thermal field, not a cold network.
    p_wu = p;
    p_wu.t_offset = t_offset_s;
    p_wu.r_q_fun  = r_q_fun;
    res_wu = simulate_plant(net, z0_cold, p_wu, @(t) 0, @(t) 1.0, T_warm);

    % Scenario run. simulate_plant initialises q from net.q0 (= cold
    % steady state at nominal flow) on every call, so we seed net.q0
    % with the warmup's final q to preserve the flow state across the
    % warmup -> scenario transition.
    net_run = net;
    net_run.q0 = res_wu.q_edges(:, end);

    p_run = p;
    p_run.t_offset = t_offset_s + T_warm;
    p_run.r_q_fun  = r_q_fun;
    T_0s_fun  = T_0s_funs{k};
    u_fun     = @(t) T_0s_fun(t) - p.Tin_nom;
    res = simulate_plant(net_run, res_wu.z_final, p_run, u_fun, @(t) 1.0, T_sim);

    % Pack results
    sc.t       = res.t;
    sc.T_0s    = p.Tin_nom + res.u;
    sc.T_0r    = res.Tout(R0_idx, :);
    sc.T_s_i   = res.T_s_i;
    sc.T_r_i   = res.T_r_i;
    sc.q_users = res.q_users;
    sc.q_F0_F1 = res.q_edges(F0_F1_edge, :);
    sc.q_edges = res.q_edges;       % full per-edge flow [n_edges x N], used by the per-junction Kirchhoff check
    sc.r_q     = res.r_q;
    sc.d_i     = res.d_i;
    sc.c_i     = res.c_i;
    out{k} = sc;

    % the |q_F0F1 - sum q_i| below is a max over the FULL run, so it picks
    % up the flow transients: each edge lags its reference with its own
    % first-order filter, which breaks trunk Kirchhoff for a few samples
    % after a flow move (about 3 percent, decaying within a few samples).
    % Conservation is asserted at steady state (T3 and T8 below).
    fprintf('  done (%d samples, max c-d %.2f W, min T_r %.2f C, max |q_F0F1 - sum q_i| %.2e)\n', ...
        numel(sc.t), max(sc.c_i(:) - sc.d_i(:)), min(sc.T_r_i(:)), ...
        max(abs(sc.q_F0_F1 - sum(sc.q_users, 1))));
end
scenarios = [out{:}];

%% Save
outdir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end
save(fullfile(outdir, 'demo.mat'), 'scenarios');
fprintf('\nSaved %s/demo.mat\n', outdir);

%% Plot
plot_plant_demo(scenarios, p);

%% Inline validation: nine physical invariants
% Each block checks one invariant on the saved scenarios. The same
% checks also live in tests/validate_plant.m as a standalone test
% (operates on the saved demo.mat without re-simulating).

fprintf('\n=== Physical invariants ===\n');

% T1: substation never delivers more than was asked for. The
% prosumer_house clips delivery at c_max; this invariant says no
% scenario sneaks past the clip.
for k = 1:numel(scenarios)
    sc = scenarios(k);
    viol = max(sc.c_i(:) - sc.d_i(:));
    assert(viol <= 1, 'T1 fail: %s c_i exceeds d_i by %.2f W', sc.name, viol);
end
fprintf('T1 pass: c^(i) <= d^(i) on every step in every scenario\n');

% T2: substation never returns water below the 5GDHC floor of 15 C.
% The clipping construction guarantees this whenever T_s >= T_r_min.
for k = 1:numel(scenarios)
    sc = scenarios(k);
    viol = p.consumer.T_r_min - min(sc.T_r_i(:));
    assert(viol <= 0.1, 'T2 fail: %s T_r_i below T_r_min by %.2f K', sc.name, viol);
end
fprintf('T2 pass: T_r^(i) >= T_r_min on every step in every scenario\n');

% T3: at steady state, total flow into the network (the trunk
% F0 -> F1) equals the sum of consumer flows. Mass cannot be created
% or destroyed inside the pipes.
for k = 1:numel(scenarios)
    sc = scenarios(k);
    n_last = max(1, numel(sc.t)-2):numel(sc.t);
    Q_net  = sum(sc.q_users(:, n_last), 1);
    res_K  = max(abs(Q_net - sc.q_F0_F1(n_last)));
    fprintf('  %-15s  Kirchhoff residual = %.2e kg/s\n', sc.name, res_K);
    assert(res_K < 1e-5, 'T3 fail: %s Kirchhoff residual %.2e kg/s', sc.name, res_K);
end
fprintf('T3 pass: trunk inflow = sum of consumer flows at steady state (residual < 1e-5)\n');

% T4: scenario-specific demand-met behaviour. Capacity-binding
% scenarios should drop below 95 percent met; adequate-capacity
% scenarios should reach at least 99 percent.
fprintf('\n  Demand-met summary\n');
for k = 1:numel(scenarios)
    sc = scenarios(k);
    pct = 100 * sum(sc.c_i(:)) / max(sum(sc.d_i(:)), 1e-9);
    fprintf('  %-15s  total met %% = %6.2f\n', sc.name, pct);
    switch sc.name
        case 'low_flow'
            assert(pct < 95, 'T4 fail: low_flow should be capacity-bound (got %.2f)', pct);
        case 'supply_step'
            assert(pct > 80 && pct < 99, ...
                'T4 fail: supply_step should be partially capacity-bound (got %.2f)', pct);
        case 'variable_flow'
            % capacity binds during the 0.5x phases that overlap demand
            % peaks; same 70-95 band as tests/validate_plant.m
            assert(pct > 70 && pct < 95, ...
                'T4 fail: variable_flow should bind at some peak (got %.2f)', pct);
        otherwise
            assert(pct >= 99, 'T4 fail: %s should reach >= 99 %% (got %.2f)', sc.name, pct);
    end
end
fprintf('T4 pass: capacity-bound vs adequate-capacity scenarios behave as expected\n');

% T5: storage-pipeline consistency. The saved c^(i) must equal
% c_p * (T_s - T_r) * q^(i) recomputed from the saved arrays (with q
% in kg/s, the AROMA convention, the factor is c_p not rho c_p).
% simulate_plant computes c_i from these same arrays, so this checks
% the save pipeline, not substation energy conservation.
fprintf('\n  Storage-pipeline consistency (c vs cp dT q)\n');
for k = 1:numel(scenarios)
    sc = scenarios(k);
    c_check = p.cp * (sc.T_s_i - sc.T_r_i) .* sc.q_users;
    res_sub = max(abs(sc.c_i - c_check), [], 'all');
    fprintf('  %-15s  max |c - cp dT q| = %.2e W\n', sc.name, res_sub);
    assert(res_sub < 1e-9, 'T5 fail: %s storage-pipeline residual %.2e W', sc.name, res_sub);
end
fprintf('T5 pass: saved c^(i) matches c_p (T_s - T_r) q^(i) (residual < 1e-9 W)\n');

% T6: network energy ratio. The thermal energy actually extracted by
% all consumers, divided by the water-side source power integrated
% over the run, should sit between 0.5 and 1.0. A perfectly mixed
% source with no losses would give 1.0; the AROMA CSTR mixing weight
% theta_1 = 0.7 puts the lower bound around
% 0.7; thermal-storage contributions push it higher.
fprintf('\n  Network energy ratio\n');
for k = 1:numel(scenarios)
    sc = scenarios(k);
    Ts = mean(diff(sc.t));
    Q_net = sum(sc.q_users, 1);
    E_in  = sum(p.cp * (sc.T_0s - sc.T_0r) .* Q_net) * Ts;
    E_ex  = sum(sc.c_i(:)) * Ts;
    ratio = E_ex / max(E_in, 1);
    fprintf('  %-15s  E_extr = %6.1f MJ  E_src = %6.1f MJ  ratio = %.3f\n', ...
        sc.name, E_ex/1e6, E_in/1e6, ratio);
    assert(ratio >= 0.5 && ratio <= 1.0, ...
        'T6 fail: %s energy ratio %.3f outside [0.5, 1.0]', sc.name, ratio);
end
fprintf('T6 pass: E_extracted / E_source_water in [0.5, 1.0] in every scenario\n');

% T7: substation upper contracts. The saved c_i is clipped at zero on
% save, so the falsifiable quantity is the raw product cp (T_s - T_r) q
% before the clip; plus T_r <= T_s (return below supply). Together with
% T1 + T2 these cover all four bounds of the substation contract.
fprintf('\n  Substation upper contracts\n');
for k = 1:numel(scenarios)
    sc = scenarios(k);
    c_raw = p.cp * (sc.T_s_i - sc.T_r_i) .* sc.q_users;
    viol_Tr_Ts = max(sc.T_r_i(:) - sc.T_s_i(:));
    fprintf('  %-15s  min raw cp dT q = %6.2f W  max(T_r - T_s) = %+6.3f K\n', ...
        sc.name, min(c_raw(:)), max(sc.T_r_i(:) - sc.T_s_i(:)));
    assert(min(c_raw(:)) > -1, 'T7 fail: %s raw cp dT q below -1 W (%.2f W)', sc.name, min(c_raw(:)));
    assert(viol_Tr_Ts <= 0.1, 'T7 fail: %s T_r above T_s by %.3f K', sc.name, viol_Tr_Ts);
end
fprintf('T7 pass: raw cp dT q > -1 W AND T_r^(i) <= T_s^(i) on every step in every scenario\n');

% T8: per-junction mass conservation on the FULL incidence matrices
% (not just the trunk that T3 covered). At steady state every
% internal supply and return node must obey Kirchhoff exactly.
fprintf('\n  Per-junction mass conservation\n');
K = build_incidence_v2(net);
for k = 1:numel(scenarios)
    sc = scenarios(k);
    n_last = max(1, numel(sc.t)-2):numel(sc.t);
    q_last = mean(sc.q_edges(:, n_last), 2);
    res_supply = max(abs(K.M_supply * q_last));
    res_return = max(abs(K.M_return * q_last));
    fprintf('  %-15s  supply residual = %.2e   return residual = %.2e   kg/s\n', ...
        sc.name, res_supply, res_return);
    assert(res_supply < 1e-9, 'T8 fail: %s supply per-junction residual %.2e', sc.name, res_supply);
    assert(res_return < 1e-9, 'T8 fail: %s return per-junction residual %.2e', sc.name, res_return);
end
fprintf('T8 pass: M_supply * q = 0 AND M_return * q = 0 at every internal node\n');

% T9: source-side return mixing law. Per-junction energy balance
% Sum_in T_in q_in = T_node Sum_out q_out evaluated at the source's
% return node. At steady state the mass-weighted average of the
% consumer return temperatures equals Q_net * T_0r, modulo small
% pipe-transport losses.
fprintf('\n  Source-side return mixing law\n');
for k = 1:numel(scenarios)
    sc = scenarios(k);
    n_last = max(1, numel(sc.t)-2):numel(sc.t);
    q_last  = sc.q_users(:, n_last);
    Tr_last = sc.T_r_i (:, n_last);
    T0r_last = sc.T_0r(n_last);
    Qnet_last = sum(q_last, 1);
    lhs = sum(q_last .* Tr_last, 1);
    rhs = Qnet_last .* T0r_last;
    rel_resid = mean(abs(lhs - rhs) ./ max(abs(rhs), 1e-9));
    fprintf('  %-15s  mean rel residual = %.3e   (tol 5e-2)\n', sc.name, rel_resid);
    assert(rel_resid < 0.05, ...
        'T9 fail: %s return mixing residual %.2e exceeds 5 %% tolerance', sc.name, rel_resid);
end
fprintf('T9 pass: source-side return mixing law holds at steady state (within pipe-loss tolerance)\n');

fprintf('\nAll nine physical invariants hold on these six scenarios.\n');


%% Local function for the variable-flow scenario
function rq = variable_flow_profile(t_rel_s, mdotE_kg)
% Hold 0.5x for 0..4 h, linear ramp to 1.5x by hour 8, hold at 1.5x
% to hour 16, ramp back to 0.5x by hour 20, hold.
t_h = t_rel_s / 3600;
if     t_h < 4
    f = 0.5;
elseif t_h < 8
    f = 0.5 + (t_h - 4) * (1.0 / 4);
elseif t_h < 16
    f = 1.5;
elseif t_h < 20
    f = 1.5 - (t_h - 16) * (1.0 / 4);
else
    f = 0.5;
end
rq = f * mdotE_kg;
end
