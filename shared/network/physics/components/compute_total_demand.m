function Q_total = compute_total_demand(net, t)
% COMPUTE_TOTAL_DEMAND   Sum of Q_demand over all prosumer nodes at time t.
%
% Used by simulate_plant in demand-driven flow mode to derive the aggregate
% mass-flow set-point:
%     w(t) = min(w_max, Q_total(t) / (cp * mdot_nom * dT_L))

Q_total = 0;
for n = 1:numel(net.Nodes)
    if strcmpi(net.Nodes(n).component, 'house')
        Q_total = Q_total + compute_prosumer_Q(t, net.Nodes(n).params);
    end
end
end
