function result = simulate_plant(net, z0, p, u_fun, w_fun, T_sim)
% SIMULATE_PLANT  Simulate the AROMA network with time-varying inputs.
%   result = simulate_plant(net, z0, p, u_fun, w_fun, T_sim)
%
%   Inputs:
%     net    - network struct from build_plant
%     z0     - initial state vector
%     p      - params struct
%     u_fun  - @(t) -> dTc [K]
%     w_fun  - @(t) -> alpha_m [-]  (used when p.flow_mode = 'exogenous')
%     T_sim  - simulation duration [s]
%
%   Outputs:
%     result.t       - [1 x N] time vector
%     result.y       - [1 x N] output (T at p.y_node)
%     result.u       - [1 x N] applied dTc
%     result.w       - [1 x N] applied alpha_m (either from w_fun or derived
%                       from aggregate demand when flow_mode = 'demand_driven')
%     result.Tout    - [n_nodes x N] all node outlet temperatures
%     result.z_final - final ODE state vector (for warm-starting next sim)
%
%   Flow-mode extensions (p.flow_mode):
%     'exogenous' (default) - w comes from w_fun as in the paper. Back-compat.
%     'demand_driven'        - w(t) = min(w_max, Q_total(t) / (cp*mdot_nom*dT_L))
%                               Adds Q_demand / Q_delivered / Q_unmet time series
%                               plus violation summary (viol_time_pct, E_unmet_kWh).
%
%   Dynamic flow path: set p.r_q_fun = @(t) r_q_vec to drive q as a state
%     via q_{k+1} = a.*q_k + b.*r_q,k. When set, w_fun is ignored for
%     flow scaling. Legacy callers without p.r_q_fun keep the old
%     scalar-w behaviour.
%
% THEORY index:
%   line 215: sampled data

Ts = p.Ts;
t_sample = 0 : Ts : T_sim;
N = numel(t_sample);
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', Ts);
% optional solver-tolerance override for numerical verification; default unchanged
if isfield(p, 'ode_reltol') && ~isempty(p.ode_reltol)
    ode_opts = odeset(ode_opts, 'RelTol', p.ode_reltol, 'AbsTol', p.ode_reltol * 1e-2);
end

% Time offset so compute_prosumer_Q sees absolute time-of-day, not
% the simulator-internal 0..Ts. The CL runner sets this every step.
t_offset = 0;
if isfield(p, 't_offset') && ~isempty(p.t_offset)
    t_offset = p.t_offset;
end

% --- flow mode --------------------------------------------------------
flow_mode = 'exogenous';
if isfield(p, 'flow_mode') && ~isempty(p.flow_mode)
    flow_mode = lower(p.flow_mode);
end
is_demand_driven = strcmp(flow_mode, 'demand_driven');

if is_demand_driven %flow is calculated by driven (other test)--> now is off
    dT_L     = p.demand.delta_T_L;
    Q_design = p.cp * p.mdot_nom * dT_L;         % [W] denominator for w
    % w_max_flow may be a scalar (constant capacity) or @(t) for a
    % time-varying capacity (renewable supply, grid curtailment, etc.).
    if isa(p.demand.w_max_flow, 'function_handle')
        w_max_fun = p.demand.w_max_flow;
    else
        w_max_fun = @(tt) p.demand.w_max_flow;
    end
end

% --- node indices -----------------------------------------------------
y_idx   = find_node(net, p.y_node);
src_idx = find_node(net, 'F0');

% --- per-edge / per-user flow indexing (the spec) --
ei      = edge_user_index(net);
n_edges = numel(net.Edges);
all_e   = 1:n_edges;

% --- per-consumer indexing -------------------------------------------
% T_s^i = outlet of the F*->C_i supply stub (= what the substation sees).
% T_r^i = post-substation node temperature at C_i.
% Reading the stub outlet from z directly avoids the small lag artefact
% you otherwise get during fast q transients (post-tap T_s drifts away
% from the temperature actually arriving at the user).
con_idx = zeros(ei.n_user, 1);
stub_off = zeros(ei.n_user, 1);   % first cell offset of the F*->C_i stub
stub_n   = zeros(ei.n_user, 1);   % number of cells in that stub
for c = 1:ei.n_user
    con_idx(c)  = find_node(net, ei.consumers{c});
    e           = ei.user(c);
    stub_off(c) = net.Edges(e).offset;
    stub_n(c)   = net.Edges(e).n;
    % inject t_offset into consumer params so the substation's
    % compute_prosumer_Q call inside rhs_network and measure_node_outlets
    % sees absolute wall-clock time even when this simulator instance
    % only spans one Ts (closed-loop runner case).
    net.Nodes(con_idx(c)).params.t_offset = t_offset;
end

% --- Dynamic flow path: discrete-time first-order q_{k+1} = a.*q_k + b.*r_q --
is_v2_flow = isfield(p, 'r_q_fun') && ~isempty(p.r_q_fun);
if is_v2_flow
    r_q_fun_v = p.r_q_fun;
    Fq_a      = net.flow_dyn.a;
    Fq_b      = net.flow_dyn.b;
    q_state   = net.q0;
end

% --- allocation -------------------------------------------------------
result.t    = t_sample;
result.y    = zeros(1, N);
result.u    = zeros(1, N);
result.w    = zeros(1, N);
result.Tout = zeros(numel(net.Nodes), N);

% per-edge mass flow [kg/s] and per-user mass flow [kg/s]
result.q_edges = zeros(n_edges, N);
result.q_users = zeros(ei.n_user, N);
result.edge_idx = ei;

% per-consumer traces: T_s, T_r, demanded d_i, realised c_i [W]
result.T_s_i = zeros(ei.n_user, N);
result.T_r_i = zeros(ei.n_user, N);
result.d_i   = zeros(ei.n_user, N);
result.c_i   = zeros(ei.n_user, N);

% when a flow reference is provided, also record what flow reference was applied
if is_v2_flow
    result.r_q = zeros(n_edges, N);
end

if is_demand_driven
    result.Q_demand    = zeros(1, N);
    result.Q_delivered = zeros(1, N);
    result.Q_unmet     = zeros(1, N);
end

% --- main loop --------------------------------------------------------
z_cur = z0;
for k = 1:N
    t_now = t_sample(k);

    % controller output u
    dTc = max(p.u_min, min(p.u_max, u_fun(t_now)));

    % flow multiplier w
    % Per-consumer profile evaluations go through prosumer_house ->
    % compute_prosumer_Q which applies prm.t_offset internally, so the
    % simulator-relative time is enough there. Absolute time is only
    % needed for the demand-driven aggregate flow path.
    if is_demand_driven
        t_abs      = t_now + t_offset;
        Q_tot      = compute_total_demand(net, t_abs);
        w_max_now  = w_max_fun(t_abs);
        Q_cap_now  = w_max_now * Q_design;
        am_raw     = Q_tot / Q_design;
        am         = min(w_max_now, max(p.w_min, am_raw));
        result.Q_demand(k)    = Q_tot;
        result.Q_delivered(k) = min(Q_tot, Q_cap_now);
        result.Q_unmet(k)     = max(0, Q_tot - Q_cap_now);
    else
        am = max(p.w_min, min(p.w_max, w_fun(t_now)));
    end

    result.u(k) = dTc;
    result.w(k) = am;

    % apply to network
    net.Nodes(src_idx).params.Tc_fun = @(t) p.Tin_nom + dTc;
    if is_v2_flow
        % freeze q_state for the duration of this sample (ZOH on flow)
        net.flows = q_state;
        r_q_now   = r_q_fun_v(t_now);
        result.r_q(:, k) = r_q_now;
        % discrete-time first-order update for the next sample
        q_state_next = Fq_a .* q_state + Fq_b .* r_q_now;
    else
        net.flows = net.mdotEdges * am;
    end

    % record per-edge flows and per-user flows q^(i) for downstream analysis
    q_now = get_edge_flows(net, all_e, t_now);
    result.q_edges(:, k) = q_now;
    result.q_users(:, k) = q_now(ei.user);

    % measure node outlets
    Tin_fun = @(t) p.Tin_nom + dTc;
    Tout = measure_node_outlets(z_cur, net, p.cp, p.mode, p.Text, Tin_fun, t_now);
    result.Tout(:, k) = Tout;
    result.y(k) = Tout(y_idx);

    % per-consumer summary: T_s read at the F*->C_i stub outlet (this is
    % what the substation actually sees), T_r at the consumer node post-
    % substation. c_i computed from those two equals the substation's
    % own clipped value because Tout(C_i) = T_s - c_substation/(q*cp).
    qi = result.q_users(:, k);
    Ts_now = zeros(ei.n_user, 1);
    for c = 1:ei.n_user
        if qi(c) >= 0
            Ts_now(c) = z_cur(stub_off(c) + stub_n(c) - 1);
        else
            Ts_now(c) = z_cur(stub_off(c));
        end
    end
    Tr_now = Tout(con_idx);
    result.T_s_i(:, k) = Ts_now;
    result.T_r_i(:, k) = Tr_now;
    result.c_i(:, k)   = max(0, qi .* p.cp .* (Ts_now - Tr_now));
    for c = 1:ei.n_user
        % pass simulator-relative t; prm.t_offset (injected above) adds
        % the wall-clock offset inside compute_prosumer_Q
        result.d_i(c, k) = compute_prosumer_Q(t_now, net.Nodes(con_idx(c)).params);
    end

    % integrate one sample
    if k < N
        odefun = @(t, z) rhs_network(t, z, net, p.cp, p.mode, p.Text, Tin_fun);
        % THEORY (sampled data): plant runs continuously, controller acts every Ts = 900 s (zero order hold)
        [~, Z] = ode45(odefun, [t_now, t_now + Ts], z_cur, ode_opts);
        z_cur = Z(end, :)';
    end

    % advance the flow state to k+1 (no-op in legacy mode)
    if is_v2_flow
        q_state = q_state_next;
    end
end

result.z_final = z_cur;

% --- violation summary for demand-driven mode -------------------------
if is_demand_driven
    dt_h = Ts / 3600;                                  % hours per sample
    result.viol_time_pct = mean(result.Q_unmet > 0) * 100;
    result.E_unmet_kWh   = sum(result.Q_unmet) * dt_h / 1000;
end
end


function idx = find_node(net, name)
idx = [];
for n = 1:numel(net.Nodes)
    if strcmp(net.Nodes(n).name, name), idx = n; return; end
end
error('find_node: "%s" not found', name);
end
