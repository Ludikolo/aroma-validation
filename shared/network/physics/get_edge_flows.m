function q = get_edge_flows(net, idx, t)
% GET_EDGE_FLOWS  Return mass flow(s) for edge index/indices idx at time t.
%
% Handles both constant (vector) and time-varying (function handle) net.flows.
%
% Inputs:  net - network struct, idx - edge index or vector of indices, t - time [s]
% Output:  q   - mass flow(s) [kg/s] for the requested edges


if isa(net.flows, 'function_handle')
    q_all = net.flows(t);     % should return one flow value per edge
else
    q_all = net.flows;        % constant flows (vector)
end

q = q_all(idx);

end
