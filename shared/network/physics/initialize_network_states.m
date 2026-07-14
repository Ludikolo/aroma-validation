function [z0, net] = initialize_network_states(net, Text)
% INITIALIZE_NETWORK_STATES  Allocate the state vector z0 for the network ODE.
%
% Assigns an offset (starting index in z0) to every pipe and to dynamic
% node components (TES/ATES/BTES). Pipe cells are initialized to ambient
% temperature Text; storage states are set to Text (or C*Text for energy-
% based storage). Also records net.n_states for downstream use.
%
% everything in net.xxx comes from network_from_digraph
%% 1) Reserve slots for all pipe cells and record their offsets
offset = 0;  % "where the next block will start" in the big state vector

% For each pipe, reserve 'n' entries and note where they start (offset)
for e = 1:numel(net.Edges)                  % net.edges = number of physical pipes
    ns = max(1, net.Edges(e).n);            % safety: at least 1 cell. net.edges.n is number of cells.
    net.Edges(e).offset = offset + 1;       % first index for this pipe
    offset = offset + ns;                   % move the cursor forward by n cells
end

%% 2) Reserve slots for nodes with their own state (TES/ATES/BTES)
% If a node has an internal state, it gets exactly 1 slot.
% Otherwise its offset stays 0 (meaning: "no state stored here").
for n = 1:numel(net.Nodes)
    comp = lower(net.Nodes(n).component);
    switch comp
        case {'tes','ates','btes'}
            net.Nodes(n).offset = offset + 1;  % this node's state index
            offset = offset + 1;               % advance by one slot
        otherwise
            net.Nodes(n).offset = 0;           % no dedicated state
    end
end

%% 3) Create z0 and fill pipe cells with ambient
z0 = zeros(offset,1);

% All pipe cells start at ambient temperature Text (safe neutral start)
for e = 1:numel(net.Edges)
    off = net.Edges(e).offset;
    ns  = net.Edges(e).n;
    z0(off:off+ns-1) = Text;
end

%% 4) Initialize dynamic node states (TES/ATES/BTES)
% TES stores temperature directly; ATES/BTES may store energy = C*T.
for n = 1:numel(net.Nodes)
    off = net.Nodes(n).offset;
    if off > 0
        comp = lower(net.Nodes(n).component);
        prm  = net.Nodes(n).params;
        switch comp
            case 'tes'
                z0(off) = Text;
            case {'ates','btes'}
                % ATES/BTES store C*T; if C is missing, default to 1
                C = 1;
                if isfield(prm,'C'), C = prm.C; end
                z0(off) = C * Text;
        end
    end
end

%% 5) Record total state count
net.n_states = numel(z0);

end
