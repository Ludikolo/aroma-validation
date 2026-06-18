function M = build_mixing(net)
% BUILD_MIXING  Find the mixing junctions for tilde g_N(z) = 0.
%
% A junction is any non-leaf node with >= 2 incoming edges, where the inflow
% balance sum_in (T_u - T_v) q_e must equal zero. T_u is the upstream-node outlet
% on each incoming edge, used as a proxy for the edge-outlet temperature (exact
% under negligible pipe loss, p.pipe.alpha = 1e-6; the residual target is built
% this way in extend_vseq_extras). In our topology this is F7 on the supply side
% and five return-side junctions R1, R2, R3, R4, R6.

n_nodes = numel(net.Nodes);
n_edges = numel(net.Edges);

% map nodes to incoming edges by name pattern
junctions = struct('name', {}, 'node_idx', {}, ...
                   'in_edges', {}, 'upstream', {}, 'is_supply', {});

for v = 1:n_nodes
    name = net.Nodes(v).name;
    if name(1) == 'C', continue; end           % consumer nodes are leaves
    if name(1) == 'F' && strcmp(name, 'F0'), continue; end  % source

    % find edges ending at v
    in_edges = [];
    upstream = [];
    for e = 1:n_edges
        if net.Edges(e).to == v
            in_edges(end+1) = e;                    %#ok<AGROW>
            upstream(end+1) = net.Edges(e).from;    %#ok<AGROW>
        end
    end

    if numel(in_edges) >= 2
        is_supply = name(1) == 'F';
        junctions(end+1).name      = name;          %#ok<AGROW>
        junctions(end).node_idx    = v;
        junctions(end).in_edges    = in_edges;
        junctions(end).upstream    = upstream;
        junctions(end).is_supply   = is_supply;
    end
end

M.junction     = junctions;
M.n_junctions  = numel(junctions);
M.n_edges      = n_edges;

end
