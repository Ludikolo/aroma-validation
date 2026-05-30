function [net, z0] = build_plant(p)
% BUILD_PLANT  Build the AROMA 5GDHC network for simulation.
%   [net, z0] = build_plant(p)
%
%   Uses shared/network/ physics code (read-only). All parameters from p = params().
%
%   Inputs:
%     p  - struct from params()
%
%   Outputs:
%     net - network struct ready for rhs_network / ODE integration
%     z0  - initial state vector
%
%   The network has 23 nodes (F0-F8, R0-R8, C1-C5) and ~736 ODE states
%   (pipe cell temperatures + dynamic source state).

dx_target = 10;  % target pipe cell length [m]

%% 1) Build AROMA topology
[nodeNames, A, nodeData, mdotEdges, pipeDiameters] = ...
    build_aroma_network_branch(p.cp, p.mdot_nom, p.mode, p.Text);

%% 2) Discretize pipes into net struct
net = network_from_digraph(A, nodeNames, nodeData, p.cp, p.Text, dx_target, pipeDiameters);

%% 3) Apply params to source
src_idx = node_index(net, 'F0');
net.Nodes(src_idx).params.kQ     = p.source.kQ;
net.Nodes(src_idx).params.V      = p.source.V;
net.Nodes(src_idx).params.theta1 = p.source.theta1;

%% 4) Scale consumer Q_peak from dT_ref, inject T_r_min, set profiles
% Profiles override the build_aroma 'diurnal' placeholder with the
% realistic residential / commercial / shop shapes from params.
Q_new = p.mdot_nom * p.cp * p.consumer.dT_ref;
Q_old = p.mdot_nom * p.cp * 10;  % build_aroma default
scale = Q_new / Q_old;
for c = 1:5
    ci = node_index(net, sprintf('C%d', c));
    net.Nodes(ci).params.Q_peak  = net.Nodes(ci).params.Q_peak * scale;
    net.Nodes(ci).params.T_r_min = p.consumer.T_r_min;
    if isfield(p, 'demand') && isfield(p.demand, 'profile_map') ...
            && numel(p.demand.profile_map) >= c
        net.Nodes(ci).params.profile = p.demand.profile_map{c};
    end
end

%% 5) Apply pipe insulation
for e = 1:numel(net.Edges)
    net.Edges(e).alpha = p.pipe.alpha;
end

%% 6) Set nominal flows and source driver
net.flows = mdotEdges;
net.Nodes(src_idx).params.Tc_fun = @(t) p.Tin_nom;

%% 7) Initialize state vector with warm-start
[z0, net] = initialize_network_states(net, p.Text);

T_fwd = p.Tin_nom;
T_ret = (p.Tin_nom + p.Text) / 2;
for e = 1:numel(net.Edges)
    off = net.Edges(e).offset;
    ns  = net.Edges(e).n;
    fn  = net.Nodes(net.Edges(e).from).name;
    tn  = net.Nodes(net.Edges(e).to).name;
    if fn(1) == 'F' || tn(1) == 'C'
        z0(off:off+ns-1) = T_fwd;
    elseif fn(1) == 'R' || fn(1) == 'C'
        z0(off:off+ns-1) = T_ret;
    end
end

%% 8) Allocate dynamic source state
prm = net.Nodes(src_idx).params;
is_dyn   = isfield(prm, 'model') && strcmpi(prm.model, 'dynamic');
has_slot = net.Nodes(src_idx).offset > 0;
if is_dyn && ~has_slot
    net.Nodes(src_idx).offset = numel(z0) + 1;
    z0(end+1, 1) = p.Tin_nom;
elseif is_dyn && has_slot
    z0(net.Nodes(src_idx).offset) = p.Tin_nom;
end
net.n_states = numel(z0);

%% 9) Store metadata
net.nodeNames  = nodeNames;
net.mdotEdges  = mdotEdges;
net.p          = p;

% Dynamic flow state: q_0 = nominal flows; Fq is the per-edge map.
net.flow_dyn = flow_dynamics(net, p);
net.q0       = mdotEdges(:);

fprintf('build_plant: %d nodes, %d edges, %d states, mode=%s\n', ...
    numel(net.Nodes), numel(net.Edges), numel(z0), p.mode);
end

function idx = node_index(net, name)
for i = 1:numel(net.Nodes)
    if strcmp(net.Nodes(i).name, name), idx = i; return; end
end
error('node_index: "%s" not found', name);
end
