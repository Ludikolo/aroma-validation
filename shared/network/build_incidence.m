function K = build_incidence(net)
% BUILD_INCIDENCE  Supply and return incidence matrices for the
%   network's mass-conservation contracts. Set up so that
%       K.M_supply * q_edges = 0
%       K.M_return * q_edges = 0
%   at every internal node (one redundant source row dropped per side
%   so the matrices are full-rank). The same matrices encode
%   spec-style A_s, A_r incidence applied to r_q (= the
%   first-order flow reference).

ei      = edge_user_index(net);
n_edges = numel(net.Edges);

% Supply nodes: every F* node + every C* node.
% Return nodes: every R* node + every C* node.
N = numel(net.Nodes);
sup_node_mask = false(N, 1);
ret_node_mask = false(N, 1);
for n = 1:N
    nm = net.Nodes(n).name;
    if nm(1) == 'F' || nm(1) == 'C', sup_node_mask(n) = true; end
    if nm(1) == 'R' || nm(1) == 'C', ret_node_mask(n) = true; end
end
sup_nodes = find(sup_node_mask);
ret_nodes = find(ret_node_mask);
n_sn = numel(sup_nodes);
n_rn = numel(ret_nodes);

sn_row = zeros(N, 1);  sn_row(sup_nodes) = 1:n_sn;
rn_row = zeros(N, 1);  rn_row(ret_nodes) = 1:n_rn;

% Pure incidence: column e is +1 at the from-node row, -1 at the
% to-node row. Built across all 29 edge columns so the supply and
% return entries simply mask out the irrelevant edges.
A_s_full = zeros(n_sn, n_edges);
A_r_full = zeros(n_rn, n_edges);
for e = 1:n_edges
    fr = net.Edges(e).from;
    to = net.Edges(e).to;
    if any(ei.Es == e)
        if sn_row(fr) > 0, A_s_full(sn_row(fr), e) = +1; end
        if sn_row(to) > 0, A_s_full(sn_row(to), e) = -1; end
    end
    if any(ei.Er == e)
        if rn_row(fr) > 0, A_r_full(rn_row(fr), e) = +1; end
        if rn_row(to) > 0, A_r_full(rn_row(to), e) = -1; end
    end
end

% RHS contributions from the user-stub flows (these are the flows
% that close the supply and return circuits).
b_s_sel = zeros(n_sn, n_edges);
b_r_sel = zeros(n_rn, n_edges);

F0_row = sn_row(find_node(net, 'F0'));
R0_row = rn_row(find_node(net, 'R0'));

for c = 1:ei.n_user
    Ci_idx = find_node(net, ei.consumers{c});
    e_user = ei.user(c);
    b_s_sel(F0_row,        e_user) = +1;
    b_s_sel(sn_row(Ci_idx), e_user) = -1;
    b_r_sel(rn_row(Ci_idx), e_user) = +1;
end

% K.M = A - b_sel  =>  K.M * q = 0 expresses the per-junction balance.
M_kirch_supply = A_s_full - b_s_sel;
M_kirch_return = A_r_full - b_r_sel;

% Drop one redundant row per side (F0 / R0) so the matrices are
% full-rank.
keep_s = setdiff(1:n_sn, F0_row);
keep_r = setdiff(1:n_rn, R0_row);

K.M_supply = M_kirch_supply(keep_s, :);
K.M_return = M_kirch_return(keep_r, :);
K.n_sn     = numel(keep_s);
K.n_rn     = numel(keep_r);
K.n_edges  = n_edges;
K.ei       = ei;
end


function idx = find_node(net, name)
for i = 1:numel(net.Nodes)
    if strcmp(net.Nodes(i).name, name), idx = i; return; end
end
error('find_node: "%s" not found', name);
end
