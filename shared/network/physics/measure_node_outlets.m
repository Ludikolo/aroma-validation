function Tnode_out = measure_node_outlets(z, net, cp, mode, Text, Tin_source, tAbs)
% MEASURE_NODE_OUTLETS  Re-create node outlet temperatures at one instant.
%
%   Tnode_out = measure_node_outlets(z, net, cp, mode, Text, Tin_source, tAbs)
%
%   This mirrors the logic in rhs_network:
%     (1) edge outlets, (2) node mixing,
%     (3) apply component physics at the node (source/house/etc.).
%
%   Extracted from log_traj_plant.m local function so it can be reused
%   by audit scripts and other tools.

E = numel(net.Edges);
N = numel(net.Nodes);

% (1) edge outlets (sign-aware)
Tout_edge = nan(E,1);
for e = 1:E
    off = net.Edges(e).offset; ns = net.Edges(e).n;
    Tcells = z(off:off+ns-1);
    q = get_edge_flows(net, e, tAbs);
    Tout_edge(e) = (q >= 0) * Tcells(end) + (q < 0) * Tcells(1);
end

% (2) node mixing (only incoming flows contribute)
Tnode_in = nan(N,1);
flow_in  = zeros(N,1);
for n = 1:N
    Qin = 0; num = 0;

    idx_to = net.Edges_to_by_node{n};
    if ~isempty(idx_to)
        q_to = get_edge_flows(net, idx_to, tAbs);
        mask = q_to > 0;
        Qin  = Qin + sum(q_to(mask));
        num  = num  + sum(q_to(mask) .* Tout_edge(idx_to(mask)));
    end

    idx_fr = net.Edges_from_by_node{n};
    if ~isempty(idx_fr)
        q_fr = get_edge_flows(net, idx_fr, tAbs);
        mask = q_fr < 0;
        Qin  = Qin + sum(-q_fr(mask));
        num  = num  + sum((-q_fr(mask)) .* Tout_edge(idx_fr(mask)));
    end

    flow_in(n)  = Qin;
    Tnode_in(n) = (Qin > 0) * (num/max(Qin,eps)) + (Qin == 0) * Text;
end

% (3) apply node components (source/house/dry_cooler/storage)
Tnode_out = Tnode_in;
for n = 1:N
    comp = lower(net.Nodes(n).component);
    prm  = net.Nodes(n).params;
    switch comp
        case 'source'
            if ~isfield(prm,'model') || strcmpi(prm.model,'dirichlet')
                Tnode_out(n) = Tin_source(tAbs);
            else
                off = net.Nodes(n).offset;
                if off <= 0
                    Tnode_out(n) = Tin_source(tAbs);
                else
                    Tmix = Tnode_in(n);
                    Tint = z(off);
                    theta1 = 0.7; if isfield(prm,'theta1'), theta1 = prm.theta1; end
                    theta1 = min(max(theta1,1e-6), 0.999999);
                    Tnode_out(n) = (Tint - (1-theta1)*Tmix)/theta1;
                end
            end

        case 'house'
            if exist('prosumer_house','file') == 2
                qloc = max(flow_in(n), 1e-12);
                Tin  = Tnode_in(n);
                [Tout, ~] = prosumer_house(tAbs, Tin, qloc, cp, prm, mode);
                Tnode_out(n) = Tout;
            else
                Tnode_out(n) = Tnode_in(n);
            end

        case 'dry_cooler'
            qloc = max(flow_in(n), 1e-12);
            if exist('dry_cooler','file') == 2
                Tnode_out(n) = dry_cooler(Tnode_in(n), qloc, cp, prm.UA, Text, prm.approach_min);
            else
                Tnode_out(n) = min(Tnode_in(n), Text);
            end

        case {'tes','ates','btes'}
            off = net.Nodes(n).offset;
            if off > 0
                if strcmp(comp,'tes'), Tnode_out(n) = z(off);
                else
                    C = 1; if isfield(prm,'C'), C = prm.C; end
                    Tnode_out(n) = z(off)/C;
                end
            end

        otherwise
            % pass-through (no component here)
    end
end
end
