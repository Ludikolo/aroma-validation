function net = apply_demand_profiles(net, profile_map)
% APPLY_DEMAND_PROFILES   Assign per-prosumer profile types on the network.
%
%   net = apply_demand_profiles(net, {'residential','commercial','shop','commercial','residential'})
%
% Sets net.Nodes(Ci).params.profile for i = 1..5 from the given cell array.
% Called by scripts that want demand-driven flow with heterogeneous profiles;
% does not affect scripts that leave profile = 'diurnal' (legacy default).

if numel(profile_map) ~= 5
    error('apply_demand_profiles: expected 5 profile strings (C1..C5)');
end

Cnames = {'C1','C2','C3','C4','C5'};
for k = 1:5
    found = false;
    for n = 1:numel(net.Nodes)
        if strcmp(net.Nodes(n).name, Cnames{k})
            net.Nodes(n).params.profile = profile_map{k};
            found = true;
            break;
        end
    end
    if ~found
        error('apply_demand_profiles: node %s not found in network', Cnames{k});
    end
end
end
