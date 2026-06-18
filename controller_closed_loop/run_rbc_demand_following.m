function res = run_rbc_demand_following(net, z0, p, T_sim)
% RUN_RBC_DEMAND_FOLLOWING  Rule-based demand-following baseline (no predictor,
% no optimisation), for the controller comparison.
%
%   res = run_rbc_demand_following(net, z0, p, T_sim)
%
% At every sample it reads the current per-consumer demand and sets the
% consumer-stub flow to the value that would just carry that demand at the
% nominal temperature drop, r_q^(i) = d^(i) / (c_p * dT_design), clipped to the
% physical flow box. The rest of the edges scale proportionally so total flow
% tracks the sum of demands. The source supply temperature is held at nominal
% (no source modulation) - this is exactly what a simple feed-forward
% demand-driven controller does. Same plant, warmup and demand realisation as
% the KPC run, so the comparison is like-for-like.
%
% p.scenario_start (absolute wall-clock time of the loop start) is used to read
% the demand profiles at the right time of day.

Ts      = p.Ts;
N_steps = round(T_sim / Ts);
ei      = edge_user_index(net);
n_user  = ei.n_user;
n_edges = numel(net.Edges);

con_idx = zeros(n_user, 1);
for i = 1:n_user
    con_idx(i) = find(strcmp({net.Nodes.name}, ei.consumers{i}));
end
F0_idx = find(strcmp({net.Nodes.name}, 'F0'));
R0_idx = find(strcmp({net.Nodes.name}, 'R0'));

dT_design = p.consumer.dT_ref;                  % nominal user temperature drop [K]
r_q_hi    = p.excite.r_q_hi_factor * net.mdotEdges(:);   % physical flow box, per edge [kg/s]
r_q_lo    = p.excite.r_q_lo_factor * net.mdotEdges(:);

% absolute time of the loop start, so compute_prosumer_Q reads the right phase
abs_base = 0;
if isfield(p, 'scenario_start') && ~isempty(p.scenario_start)
    abs_base = p.scenario_start;
end

% storage, same field layout as a KPC closed-loop trajectory
res.t       = (0:N_steps) * Ts;
res.Tout    = zeros(numel(net.Nodes), N_steps + 1);
res.q_edges = zeros(n_edges, N_steps + 1);
res.q_users = zeros(n_user,  N_steps + 1);
res.r_q     = zeros(n_edges, N_steps + 1);
res.T_0s    = zeros(1, N_steps + 1);
res.T_0r    = zeros(1, N_steps + 1);
res.T_F0    = zeros(1, N_steps + 1);
res.T_is    = zeros(n_user, N_steps + 1);
res.T_ir    = zeros(n_user, N_steps + 1);
res.d_i     = zeros(n_user, N_steps + 1);
res.c_i     = zeros(n_user, N_steps + 1);

z_state = z0;
for k = 1:N_steps
    t_now = (k-1) * Ts;

    % per-user demand at the current absolute time, and the flow that carries it
    d_now = zeros(n_user, 1);
    for i = 1:n_user
        d_now(i) = compute_prosumer_Q(abs_base + t_now, net.Nodes(con_idx(i)).params);
    end
    r_q_user = d_now / (p.cp * dT_design);
    for i = 1:n_user
        e = ei.user(i);
        r_q_user(i) = min(r_q_hi(e), max(r_q_lo(e), r_q_user(i)));
    end

    % full reference: scale every edge by the total-user-flow ratio, then pin
    % the consumer stubs to their targets, then clip to the physical box
    scale    = sum(r_q_user) / max(sum(net.mdotEdges(ei.user)), 1e-9);
    r_q_full = net.mdotEdges(:) * scale;
    r_q_full(ei.user) = r_q_user;
    r_q_full = max(r_q_lo, min(r_q_hi, r_q_full));

    % one plant step at this flow reference, source held at nominal (u = 0)
    p_step = p;
    p_step.r_q_fun  = @(t) r_q_full;
    p_step.t_offset = abs_base + t_now;
    net_step = net;
    if k == 1, net_step.q0 = net.mdotEdges(:); end
    res_step = simulate_plant(net_step, z_state, p_step, @(t) 0, @(t) 1.0, Ts);

    res.Tout(:, k+1)    = res_step.Tout(:, end);
    res.q_edges(:, k+1) = res_step.q_edges(:, end);
    res.q_users(:, k+1) = res_step.q_users(:, end);
    res.r_q(:, k+1)     = res_step.r_q(:, end);
    res.T_0s(k+1)       = p.Tin_nom + res_step.u(end);
    res.T_0r(k+1)       = res_step.Tout(R0_idx, end);
    res.T_F0(k+1)       = res_step.Tout(F0_idx, end);
    res.T_is(:, k+1)    = res_step.T_s_i(:, end);
    res.T_ir(:, k+1)    = res_step.T_r_i(:, end);
    res.d_i(:, k+1)     = res_step.d_i(:, end);
    res.c_i(:, k+1)     = res_step.c_i(:, end);

    z_state = res_step.z_final;
end

% drop the k=1 column (initial state, before any move) so the window matches
% the KPC step loop's convention
fn = {'t','Tout','q_edges','q_users','r_q','T_0s','T_0r','T_F0','T_is','T_ir','d_i','c_i'};
for j = 1:numel(fn)
    res.(fn{j}) = res.(fn{j})(:, 2:end);
end
end
