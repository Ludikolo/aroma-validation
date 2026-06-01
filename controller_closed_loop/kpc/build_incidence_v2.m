function K = build_incidence_v2(net)
% BUILD_INCIDENCE_V2  Supply / return incidence matrices for Kirchhoff
% conservation, set up to act on the flow-reference r_q.
%
% Imposing A r_q = b at every horizon together with q_0 conservative is
% sufficient to keep q conservative for all k under q_{k+1} = a q_k +
% (1-a) r_q,k. This avoids fitting V_seq for every edge flow.

ei      = edge_user_index(net);
n_edges = numel(net.Edges);

% supply graph: nodes F0..F8, C1..C5 (14 nodes)
% return graph: nodes R0..R8, C1..C5 (14 nodes)
N      = numel(net.Nodes);
sup_node_mask = false(N, 1);
ret_node_mask = false(N, 1);
for n = 1:N
    nm = net.Nodes(n).name;
    if nm(1) == 'F' || nm(1) == 'C', sup_node_mask(n) = true; end
    if nm(1) == 'R' || nm(1) == 'C', ret_node_mask(n) = true; end
end
sup_nodes = find(sup_node_mask);
ret_nodes = find(ret_node_mask);
n_sn      = numel(sup_nodes);
n_rn      = numel(ret_nodes);

sn_row = zeros(N, 1);  sn_row(sup_nodes) = 1:n_sn;
rn_row = zeros(N, 1);  rn_row(ret_nodes) = 1:n_rn;

% pure incidence: A(:, e) = +1 at from-node row, -1 at to-node row.
% Build for the FULL r_q index space (n_edges columns) so the slicing
% below is just a regular linear matrix multiply.
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

% RHS of Kirchhoff depends on q^(i), which is a slice of r_q (the user
% stubs). Build a selector b_sel so that the RHS is b_sel * r_q, then
% the constraint becomes (A - b_sel) r_q = 0.
b_s_sel = zeros(n_sn, n_edges);
b_r_sel = zeros(n_rn, n_edges);

% find the F0 and R0 rows
F0_row = sn_row(find_node(net, 'F0'));
R0_row = rn_row(find_node(net, 'R0'));

for c = 1:ei.n_user
    Ci      = ei.consumers{c};
    e_user  = ei.user(c);                % the F* -> C_i supply stub edge
    Ci_idx  = find_node(net, Ci);

    % supply: F0 row gets +1 for this user's r_q (Q_net contribution).
    %         C_i row gets -1 for this user's r_q (consumer outflow).
    b_s_sel(F0_row,        e_user) = +1;
    b_s_sel(sn_row(Ci_idx), e_user) = -1;

    % return: R0 row stays zero because the closure absorbs Q_net here
    %         (R0 row of A_r_full sums to 0 already). C_i row gets +1.
    b_r_sel(rn_row(Ci_idx), e_user) = +1;
end

% per-horizon Kirchhoff matrix: M r_q = 0 with M = A - b_sel
M_kirch_supply = A_s_full - b_s_sel;       % 14 x n_edges
M_kirch_return = A_r_full - b_r_sel;       % 14 x n_edges

% Drop one redundant row per network (rank n_nodes - 1): F0 in supply,
% R0 in return. Otherwise Aeq is rank-deficient.
keep_s = setdiff(1:n_sn, F0_row);
keep_r = setdiff(1:n_rn, R0_row);

K.M_supply = M_kirch_supply(keep_s, :);    % (n_sn - 1) x n_edges
K.M_return = M_kirch_return(keep_r, :);    % (n_rn - 1) x n_edges
K.n_sn     = numel(keep_s);                % rank-1 reduced
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
