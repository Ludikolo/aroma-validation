function idx = edge_user_index(net)
% EDGE_USER_INDEX  Supply / return / user edge sets for v2.
%
% Three edge groups in the spec: E_s (supply), E_r (return),
% one user-flow q^(i) per consumer carried by the F*->C_i stub.
% I = {2,3,5,6,8} maps to our C1..C5 at taps
% F2/F3/F5/F6/F8.
%
% The order of idx.user(:) fixes the column order of q^(i) in the
% downstream pipeline; do not shuffle without re-running training.

n_edges = numel(net.Edges);

is_supply = false(n_edges, 1);
is_return = false(n_edges, 1);

% one supply-stub edge per consumer (Ci -> edge index)
n_users   = 5;
user_edge = zeros(n_users, 1);

for e = 1:n_edges
    fn = net.Nodes(net.Edges(e).from).name;
    tn = net.Nodes(net.Edges(e).to).name;

    f0 = fn(1);
    t0 = tn(1);

    if (f0 == 'F' && t0 == 'F') || (f0 == 'F' && t0 == 'C')
        % forward ring or supply stub into a consumer
        is_supply(e) = true;
        if t0 == 'C'
            ci = sscanf(tn, 'C%d');
            user_edge(ci) = e;
        end
    elseif (f0 == 'R' && t0 == 'R') || (f0 == 'C' && t0 == 'R') || (f0 == 'R' && t0 == 'F')
        % return ring, return stub from a consumer, or R0->F0 loop closure
        is_return(e) = true;
    else
        error('edge_user_index: edge %s->%s does not fit supply/return pattern', fn, tn);
    end
end

if any(user_edge == 0)
    error('edge_user_index: missing supply stub for at least one consumer');
end

% I = {2,3,5,6,8} attached at F2/F3/F5/F6/F8
idx.consumers     = {'C1','C2','C3','C4','C5'};
idx.tap_nodes     = {'F2','F3','F5','F6','F8'};

idx.Es   = find(is_supply);    % column of edge indices
idx.Er   = find(is_return);
idx.user = user_edge;          % length 5: idx.user(c) = edge index for q^(C_c)

idx.n_es   = numel(idx.Es);
idx.n_er   = numel(idx.Er);
idx.n_user = n_users;

end
