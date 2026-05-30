function traj = generate_traj(net, z0, p, T_sim, seed)
% GENERATE_V2_TRAJ  One v2 plant trajectory under multi-channel
% PRBS excitation. Returns a struct with every signal the dictionary
% might want.

[r_q_fun, T0s_fun, ~] = excitation_multi(net, p, T_sim, seed);

% v2 mode: simulator takes the flow reference via p.r_q_fun
p.r_q_fun = r_q_fun;
u_fun = @(t) T0s_fun(t);
w_fun = @(t) 1.0;        % ignored in v2 mode but the interface needs it

res = simulate_plant(net, z0, p, u_fun, w_fun, T_sim);

F0_idx = find_node(net, 'F0');
R0_idx = find_node(net, 'R0');

% pack the trajectory
traj.t        = res.t;                                       % 1 x N
traj.Tout     = res.Tout;                                    % n_nodes x N (full nodal temps)
traj.q_edges  = res.q_edges;                                 % n_edges x N
traj.q_users  = res.q_users;                                 % |I| x N
traj.r_q      = res.r_q;                                     % n_edges x N
traj.T_0s     = p.Tin_nom + res.u;                           % 1 x N (commanded supply)
traj.T_0r     = res.Tout(R0_idx, :);                         % 1 x N (return at plant)
traj.T_F0     = res.Tout(F0_idx, :);                         % 1 x N (post-source supply)
traj.T_is     = res.T_s_i;                                   % |I| x N
traj.T_ir     = res.T_r_i;                                   % |I| x N
traj.d        = res.d_i;                                     % |I| x N
traj.c        = res.c_i;                                     % |I| x N
traj.theta    = (res.T_s_i - res.T_r_i) .* res.q_users;      % |I| x N (the bilinear obs.)

% bookkeeping
traj.seed     = seed;
traj.T_sim    = T_sim;
traj.Ts       = p.Ts;
traj.edge_idx = res.edge_idx;

end


function idx = find_node(net, name)
for i = 1:numel(net.Nodes)
    if strcmp(net.Nodes(i).name, name), idx = i; return; end
end
error('find_node: "%s" not found', name);
end
