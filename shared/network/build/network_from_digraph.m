function net = network_from_digraph(A, nodeNames, nodeData, cp, Text, dx_target, pipeDiameters)
% NETWORK_FROM_DIGRAPH  Build a simulation-ready network struct from an adjacency matrix.
%
% Every nonzero A(i,j) = L defines a pipe from node i to j with length L [m].
% Each pipe is discretized into cells of approximately dx_target [m].
% The edge direction only sets cell index order; actual flow direction is
% determined at runtime by the sign of the mass flow.
%
% INPUTS:
%   A            - Sparse adjacency matrix, A(i,j) = pipe length [m]
%   nodeNames    - Cell array of node name strings
%   nodeData     - containers.Map with node component data
%   cp           - Specific heat capacity [J/(kg·K)]
%   Text         - Ambient temperature [°C]
%   dx_target    - Target cell length for discretization [m]
%   pipeDiameters- (optional) containers.Map, keys 'FromNode>ToNode', values [m]
%                  Falls back to uniform D = 0.1 m (DN100) if not provided
%
% OUTPUT: net struct with Nodes, Edges, Edges_to_by_node, Edges_from_by_node,
%   flows (placeholder), cp, Text, rho, dx_target.

%% 1) Input checks
if nargin < 6
    error('dx_target is required.');
end
if ~issparse(A)
    A = sparse(A);
end

% Count nodes and extract edges from the adjacency matrix
N = numel(nodeNames);
[I, J, W] = find(A);   % I = from-node idx, J = to-node idx, W = pipe length [m]
E = numel(I);

%% 2) Create node objects
net.Nodes = repmat(struct('name','', 'component','', 'params',struct(), 'offset',0), N, 1);
for n = 1:N
    net.Nodes(n).name = nodeNames{n};
    % If there is metadata for this node, attach it (component type + params)
    if isa(nodeData,'containers.Map') && isKey(nodeData, nodeNames{n})
        nd = nodeData(nodeNames{n});
        net.Nodes(n).component = nd.type;    % e.g., 'source', 'house'
        net.Nodes(n).params    = nd.params;  % whatever that component needs
    end
end

%% 3) Create edge objects (one per nonzero A(i,j))
% Check if per-pipe diameters were provided
usePipeDiameters = nargin >= 7 && isa(pipeDiameters, 'containers.Map');

net.Edges = repmat(struct('from',0,'to',0,'L',0,'n',0,'D',0,'Acs',0,'alpha',0,'offset',0), E, 1);
for e = 1:E
    % Who is connected?
    net.Edges(e).from = I(e);
    net.Edges(e).to   = J(e);
    net.Edges(e).L    = W(e);   % length [m]

    % Number of cells: aim for dx_target, minimum 1
    nCells = max(1, round(W(e)/dx_target));
    net.Edges(e).n = nCells;

    % Pipe diameter: use per-pipe map if provided, else default DN100
    if usePipeDiameters
        key = sprintf('%s>%s', nodeNames{I(e)}, nodeNames{J(e)});
        if isKey(pipeDiameters, key)
            D = pipeDiameters(key);
        else
            warning('No diameter for edge %s, using default 0.1m', key);
            D = 0.1;
        end
    else
        % Fallback: uniform diameter (DN100 = 0.1m)
        D = 0.1;
    end
    Acs = pi*(D/2)^2;       % [m^2] cross-sectional area
    net.Edges(e).D   = D;
    net.Edges(e).Acs = Acs;

    % Ambient exchange coefficient [1/s]
    net.Edges(e).alpha = 1e-5;
end

%% 4) Precompute adjacency lists for fast node mixing
net.Edges_from = arrayfun(@(e) e.from, net.Edges(:));
net.Edges_to   = arrayfun(@(e) e.to,   net.Edges(:));
net.Edges_to_by_node   = arrayfun(@(n) find(net.Edges_to   == n), 1:N, 'UniformOutput', false).';
net.Edges_from_by_node = arrayfun(@(n) find(net.Edges_from == n), 1:N, 'UniformOutput', false).';

%% 5) Store global constants and placeholders
% net.flows is a dummy; the caller overwrites it with real constant or time-varying flows.
net.cp        = cp;           % specific heat [J/kg/K]
net.Text      = Text;         % ambient temperature [°C]
net.dx_target = dx_target;    % desired cell length [m]
net.rho       = 1000;         % water density [kg/m^3]
net.flows     = ones(E,1);    % placeholder; to be set by caller (main script)
net.nodeData  = nodeData;     % keep original node metadata for reference

end
