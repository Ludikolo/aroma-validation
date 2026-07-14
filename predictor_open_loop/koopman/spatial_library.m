% NOTE: reproduces the spatial dictionary comparison in dictionary_study.pdf; not on the run_all_tests path.
function [Z, names, meta] = spatial_library(traj, p, net, variant)
% SPATIAL_LIBRARY  Supply-ring spatial lift for the dictionary study.
%
% Wraps candidate_library (rows 1..47 unchanged) and appends the in-network
% supply-ring signals, following the masked-bilinear construction of the
% prior Koopman study, restricted to the forward ring:
%   S1 node temps   T_F1..T_F8                                     8
%   S2 ring flows   q_e for the F->F ring edges                    9
%   S3 masked products q_e * T_from(e) for every supply edge      14
% giving n_z = 47 + 31 = 78. The node temperatures are junction mixing
% temperatures (measure_node_outlets), so every added signal is one a real
% network could measure at a manifold; no in-pipe states are used.
%
% variant = 'full' (default) keeps the production lift as rows 1..47.
% variant = 'nolag' drops the source-delay block from the base lift before
% appending (n_z = 35 + 31 = 66); used only for the lag ablation, to test
% whether the spatial block makes the T_F0 lags redundant.
% variant = 'return' is 'full' plus the return-ring junction temperatures
% T_R1..T_R8 as linear features (n_z = 78 + 8 = 86); no return products,
% since the C->R stub products are spanned by the exergy block already.
% variant = 'none' appends nothing (n_z = 47): the production lift under the
% same fit protocol, used as the in-study reference of the dictionary ladder.
%
% net is optional (built from p when absent); pass it in fit loops so the
% plant is not rebuilt per trajectory.

if nargin < 3 || isempty(net)
    net = build_plant(p);
end
if nargin < 4 || isempty(variant)
    variant = 'full';
end

p_base = p;
if strcmp(variant, 'nolag')
    p_base.o2.use_delays = false;      % base lift without the delay block
end

[Z0, nm0, meta] = candidate_library(traj, p_base);

if strcmp(variant, 'none')
    Z = Z0; names = nm0;
    meta.dict = 'spatial_none';
    return;
end

% supply-side topology from the network itself (parametric, no hardcoding)
ei   = edge_user_index(net);
from = net.Edges_from;

% forward-ring nodes F1..F8 by name; F0 is already lift row idx.T_F0
sup_nodes = [];
for n = 1:numel(net.Nodes)
    nn = net.Nodes(n).name;
    if nn(1) == 'F' && ~strcmp(nn, 'F0')
        sup_nodes(end+1) = n; %#ok<AGROW>
    end
end

ring_edges = setdiff(ei.Es, ei.user);   % F->F ring edges (no consumer stubs)

feats = {};
nm    = {};

% S1: junction temperatures on the forward ring
idx_T = zeros(1, numel(sup_nodes));
for k = 1:numel(sup_nodes)
    n = sup_nodes(k);
    feats{end+1} = traj.Tout(n, :);
    nm{end+1}    = sprintf('T_%s', net.Nodes(n).name);
    idx_T(k) = meta.n_z + numel(feats);
end

% S2: realized ring flows (states, not the r_q input)
idx_q = zeros(1, numel(ring_edges));
for k = 1:numel(ring_edges)
    e = ring_edges(k);
    feats{end+1} = traj.q_edges(e, :);
    nm{end+1}    = sprintf('q_%s_%s', net.Nodes(net.Edges(e).from).name, ...
                                      net.Nodes(net.Edges(e).to).name);
    idx_q(k) = meta.n_z + numel(feats);
end

% S3: masked bilinears q_e * T_from(e), one per supply edge (outgoing mask:
% the product exists only where the edge physically leaves the node)
idx_P = zeros(1, numel(ei.Es));
for k = 1:numel(ei.Es)
    e = ei.Es(k);
    feats{end+1} = traj.q_edges(e, :) .* traj.Tout(from(e), :);
    nm{end+1}    = sprintf('qT_%s_%s', net.Nodes(net.Edges(e).from).name, ...
                                       net.Nodes(net.Edges(e).to).name);
    idx_P(k) = meta.n_z + numel(feats);
end

% optional return-ring junction temperatures R1..R8 (R0 is lift row idx.T_0r)
idx_R = [];
if strcmp(variant, 'return')
    for n = 1:numel(net.Nodes)
        nn = net.Nodes(n).name;
        if nn(1) == 'R' && ~strcmp(nn, 'R0')
            feats{end+1} = traj.Tout(n, :);
            nm{end+1}    = sprintf('T_%s', nn);
            idx_R(end+1) = meta.n_z + numel(feats); %#ok<AGROW>
        end
    end
end

Z     = [Z0; vertcat(feats{:})];
names = [nm0, nm];

meta.idx.spatial_T = idx_T;
meta.idx.spatial_q = idx_q;
meta.idx.spatial_P = idx_P;
if ~isempty(idx_R)
    meta.idx.spatial_R = idx_R;
end
meta.n_z  = size(Z, 1);
meta.dict = ['spatial_' variant];

end
