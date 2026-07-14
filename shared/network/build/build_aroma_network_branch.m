function [nodeNames, A, nodeData, mdotEdges, pipeDiameters] = build_aroma_network_branch(cp, mdot, mode, Text, varargin)
% BUILD_AROMA_NETWORK_BRANCH  Build the AROMA 5GDHC ring network.
%
% One producer (F0) and five consumers (C1..C5) on forward/return rings
% with short 1m consumer stubs. Returns node names, adjacency matrix A
% (pipe lengths in m), node component data, edge mass flows (matching
% find(A) order), and per-pipe diameters from AROMA definition.pdf Table 2.
%
% Pipe diameters (AROMA PDF Table 2):
%   107 mm: F0-F1, F1-F2, F1-F6, R1-R0, R2-R1, R6-R1
%    83 mm: F2-F3, F3-F4, F4-F7, F6-F7, R3-R2, R4-R3, R7-R4, R7-R6
%    70 mm: F4-F5, F7-F8, R5-R4, R8-R7, stubs (assumed)

%% 1) Inputs and defaults
if nargin < 1 || isempty(cp),   cp   = 4180;      end
if nargin < 2 || isempty(mdot), mdot = 3.0;       end
if nargin < 3 || isempty(mode), mode = 'heating'; end
if nargin < 4 || isempty(Text), Text = 15;        end %#ok<NASGU> kept for the build_plant interface; ambient enters via network_from_digraph

% Optional: how to split C5 branch flow over two forward paths
% We have a split at node F1 and a split at node R7.
p = inputParser;
p.addParameter('alpha_split', 0.5, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
p.parse(varargin{:});
alpha_split = p.Results.alpha_split;

%% 2) Consumer weights, node list, and adjacency matrix A (meters)
% Relative heat demand weights, normalized to sum to 1.
w = struct('F2',0.11,'F3',0.34,'F5',0.08,'F6',0.38,'F8',0.08);
wSum = w.F2 + w.F3 + w.F5 + w.F6 + w.F8;            % ≈ 0.99 in table
w.F2 = w.F2 / wSum; w.F3 = w.F3 / wSum;
w.F5 = w.F5 / wSum; w.F6 = w.F6 / wSum; w.F8 = w.F8 / wSum;

% Reference loop power for a 10 K temperature drop at mdot
Q_total_W = mdot * cp * 10;

% Node names: forward F0..F8, return R0..R8, consumers C1..C5.
nodeNames = [ arrayfun(@(k) sprintf('F%d',k), 0:8, 'UniformOutput', false), ...
              arrayfun(@(k) sprintf('R%d',k), 0:8, 'UniformOutput', false), ...
              {'C1','C2','C3','C4','C5'} ];
N  = numel(nodeNames);
ix = containers.Map(nodeNames, 1:N);               % fast name -> index

% Build sparse adjacency A where A(i,j) = pipe length [m] from i to j.
A = sparse(N, N);

% Forward ring + branches
A(ix('F0'),ix('F1')) = 500.0;
A(ix('F1'),ix('F2')) = 282.8;
A(ix('F2'),ix('F3')) = 500.0;
A(ix('F3'),ix('F4')) = 282.8;
A(ix('F4'),ix('F5')) = 400.0;
A(ix('F4'),ix('F7')) = 282.8;
A(ix('F1'),ix('F6')) = 282.8;
A(ix('F6'),ix('F7')) = 500.0;
A(ix('F7'),ix('F8')) = 600.0;

% Return ring
A(ix('R1'),ix('R0')) = 500.0;
A(ix('R2'),ix('R1')) = 282.8;
A(ix('R3'),ix('R2')) = 500.0;
A(ix('R4'),ix('R3')) = 282.8;
A(ix('R5'),ix('R4')) = 400.0;
A(ix('R7'),ix('R4')) = 282.8;
A(ix('R6'),ix('R1')) = 282.8;
A(ix('R7'),ix('R6')) = 500.0;
A(ix('R8'),ix('R7')) = 600.0;

% Loop closure + consumer stubs (1 m so houses sit exactly on nodes)
L_stub = 1.0;
A(ix('R0'),ix('F0')) = L_stub;

A(ix('F2'),ix('C1')) = L_stub;  A(ix('C1'),ix('R2')) = L_stub;
A(ix('F3'),ix('C2')) = L_stub;  A(ix('C2'),ix('R3')) = L_stub;
A(ix('F5'),ix('C3')) = L_stub;  A(ix('C3'),ix('R5')) = L_stub;
A(ix('F6'),ix('C4')) = L_stub;  A(ix('C4'),ix('R6')) = L_stub;
A(ix('F8'),ix('C5')) = L_stub;  A(ix('C5'),ix('R8')) = L_stub;

%% 3) Node components: dynamic source at F0 + five houses at C1..C5
% Map that tells which component sits at each node and its parameters.
nodeData = containers.Map;

% Source at F0 (dynamic model by default; values can be overwritten later).
nodeData('F0') = struct('type','source','params', struct( ...
    'model','dynamic', ...
    'V',      0.05, ...   % m^3
    'rho',    1000, ...
    'cp',     cp, ...
    'theta1', 0.70, ...
    'kQ',     4000, ...   % W/K
    'Tc_fun', []));

% Houses C1..C5: size them by weight and give a simple diurnal profile.
Cnames = {'C1','C2','C3','C4','C5'};
wlist  = [w.F2, w.F3, w.F5, w.F6, w.F8];
for k = 1:5
    nodeData(Cnames{k}) = struct('type','house','params', struct( ...
        'Q_peak',  wlist(k) * Q_total_W, ...  % [W] house size by weight
        'profile', 'diurnal', ...              % simple placeholder
        'mode',    mode));                     % 'cooling' or 'heating'
end

%% 4) Constant edge flows satisfying mass balance
q2 = w.F2 * mdot;
q3 = w.F3 * mdot;
q5 = w.F5 * mdot;
q6 = w.F6 * mdot;
q8 = w.F8 * mdot;

% Split C5 branch flow q8 over two paths using alpha_split.
qa = alpha_split       * q8; % via F1->F6->F7
qb = (1 - alpha_split) * q8; % via F0->...->F4->F7

flows = containers.Map('KeyType','char','ValueType','double');

% Forward
flows('F0>F1') = mdot;
flows('F1>F2') = mdot - (q6 + qa);
flows('F2>F3') = mdot - (q6 + qa) - q2;
flows('F3>F4') = mdot - (q6 + qa) - q2 - q3;
flows('F4>F5') = q5;
flows('F4>F7') = qb;
flows('F1>F6') = q6 + qa;
flows('F6>F7') = qa;
flows('F7>F8') = q8;

% Return
flows('R8>R7') = q8;
flows('R7>R6') = qa;
flows('R7>R4') = qb;
flows('R6>R1') = q6 + qa;
flows('R5>R4') = q5;
flows('R4>R3') = q5 + qb;
flows('R3>R2') = q3 + q5 + qb;
flows('R2>R1') = q2 + q3 + q5 + qb;
flows('R1>R0') = mdot;

% Closure
flows('R0>F0') = mdot;

% Stubs (consumer legs carry their own branch flow)
flows('F2>C1') = q2; flows('C1>R2') = q2;
flows('F3>C2') = q3; flows('C2>R3') = q3;
flows('F5>C3') = q5; flows('C3>R5') = q5;
flows('F6>C4') = q6; flows('C4>R6') = q6;
flows('F8>C5') = q8; flows('C5>R8') = q8;

%% 5) Convert flows Map -> mdotEdges vector in find(A) order
% Edge order must match [I,J,~] = find(A).
[I, J, ~] = find(A); E = numel(I);
mdotEdges = zeros(E,1);
for e = 1:E
    key = sprintf('%s>%s', nodeNames{I(e)}, nodeNames{J(e)});
    mdotEdges(e) = flows(key);   % should exist for every nonzero A(i,j)
end

%% 6) Per-pipe diameter map (AROMA PDF Table 2)

D_107 = 0.107;  % DN100 equivalent
D_83  = 0.083;  % DN80 equivalent
D_70  = 0.070;  % DN65 equivalent

pipeDiameters = containers.Map('KeyType','char','ValueType','double');

% ---- Da = 107 mm (0.107 m) pipes ----
% Forward ring main trunk
pipeDiameters('F0>F1') = D_107;
pipeDiameters('F1>F2') = D_107;
pipeDiameters('F1>F6') = D_107;
% Return ring main trunk
pipeDiameters('R1>R0') = D_107;
pipeDiameters('R2>R1') = D_107;
pipeDiameters('R6>R1') = D_107;

% ---- Da = 83 mm (0.083 m) pipes ----
% Forward ring secondary
pipeDiameters('F2>F3') = D_83;
pipeDiameters('F3>F4') = D_83;
pipeDiameters('F4>F7') = D_83;
pipeDiameters('F6>F7') = D_83;
% Return ring secondary
pipeDiameters('R3>R2') = D_83;
pipeDiameters('R4>R3') = D_83;
pipeDiameters('R7>R4') = D_83;
pipeDiameters('R7>R6') = D_83;

% ---- Da = 70 mm (0.070 m) pipes ----
% Forward ring terminal
pipeDiameters('F4>F5') = D_70;
pipeDiameters('F7>F8') = D_70;
% Return ring terminal
pipeDiameters('R5>R4') = D_70;
pipeDiameters('R8>R7') = D_70;

% ---- Consumer stubs & loop closure: 70 mm (assumed) ----
pipeDiameters('R0>F0') = D_70;  % Loop closure

pipeDiameters('F2>C1') = D_70;  pipeDiameters('C1>R2') = D_70;
pipeDiameters('F3>C2') = D_70;  pipeDiameters('C2>R3') = D_70;
pipeDiameters('F5>C3') = D_70;  pipeDiameters('C3>R5') = D_70;
pipeDiameters('F6>C4') = D_70;  pipeDiameters('C4>R6') = D_70;
pipeDiameters('F8>C5') = D_70;  pipeDiameters('C5>R8') = D_70;

%% 7) Mass balance check
if abs(q2 + q3 + q5 + q6 + q8 - mdot) > 1e-9
    warning('Sum of consumer flows ~= mdot (%.6f vs %.6f).', q2+q3+q5+q6+q8, mdot);
end

% Sanity check: every edge should have a diameter defined
[Ie, Je, ~] = find(A);
for ee = 1:numel(Ie)
    key = sprintf('%s>%s', nodeNames{Ie(ee)}, nodeNames{Je(ee)});
    if ~isKey(pipeDiameters, key)
        warning('Missing diameter for edge %s. Using default 70mm.', key);
        pipeDiameters(key) = D_70;
    end
end

end
