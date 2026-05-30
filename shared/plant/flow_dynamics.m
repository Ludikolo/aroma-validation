function Fq = flow_dynamics(net, p)
% FLOW_DYNAMICS  Discrete first-order map per edge:
%   q_{k+1}^(*) = a^(*) q_k^(*) + (1 - a^(*)) r_q,k^(*)
%   a^(*) = exp(-Ts / eps^(*))
%
% Three groups: supply edges, return edges, user-side stubs.

ei = edge_user_index(net);
n_edges = numel(net.Edges);

a = zeros(n_edges, 1);
a(ei.Es) = exp(-p.Ts / p.flow_dyn.epsilon_supply);
a(ei.Er) = exp(-p.Ts / p.flow_dyn.epsilon_return);

% per-user stubs override the supply default with the (faster) user tau
for c = 1:ei.n_user
    a(ei.user(c)) = exp(-p.Ts / p.flow_dyn.epsilon_user);
end

if any(a == 0)
    error('flow_dynamics: edge %d got no group assignment', find(a == 0, 1));
end

Fq.a  = a;
Fq.b  = 1 - a;
Fq.Ts = p.Ts;
Fq.eps_supply = p.flow_dyn.epsilon_supply;
Fq.eps_return = p.flow_dyn.epsilon_return;
Fq.eps_user   = p.flow_dyn.epsilon_user;

end
