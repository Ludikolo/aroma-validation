function dz = rhs_network(t, z, net, cp, mode, Text, Tin_source)
% RHS_NETWORK  ODE right-hand side for the AROMA 5GDHC network.
%
% This is the function that ode45/ode15s calls at every time step.
% It computes dz/dt for the full network state vector z, which contains
% all pipe cell temperatures plus a few node states (TES, source, etc.).
%
% Core physics for each pipe cell (first-order upwind advection + ambient loss):
%   dT/dt = -v * dT/dx - UA/(m*cp) * (T - T_amb)
% where v = mdot/(rho*A_cs) is the bulk velocity, and the spatial
% derivative dT/dx is discretised with an upwind finite-difference scheme.
%
% The rest of the function handles node mixing (flow-weighted averaging of
% incoming streams) and component models (source, house, dry cooler,
% TES/ATES/BTES) that provide boundary conditions to the pipe cells.
%
% Pipes move heat via advection; nodes mix incoming flows and apply device models.
% We return dz so the ODE solver can step forward.
%
% THEORY index:
%   line 49: node mixing
%   line 138: pipe advection
%   line 216: source


E = numel(net.Edges);
N = numel(net.Nodes);

%% 1) Edge outlet temperatures (sign-aware)
% for each pipe, we look at the last cell if flow is
% positive (going i -> j), or the first cell if flow is negative
% (going j -> i). That chosen cell is the pipe's outlet temperature.
Tout_edge = nan(E,1);
for e = 1:E
    off    = net.Edges(e).offset;
    ns     = net.Edges(e).n;
    Tcells = z(off:off+ns-1);

    q = get_edge_flows(net, e, t);                  % current flow for this edge
    Tout_edge(e) = (q >= 0) * Tcells(end) + ...     % outlet is last cell if q >= 0
                   (q <  0) * Tcells(1);            % outlet is first cell if q < 0
end

% THEORY (node mixing): outlet temp = flow-weighted average of incoming streams
%% 2) Node mixing (incoming streams only)
% Flow-weighted average of all streams entering the node.
% Falls back to ambient if no inflow.
Tnode_in = nan(N,1);
flow_in  = zeros(N,1);
for n = 1:N
    Qin = 0; num = 0;

    % Edges that end at node n (i -> n). Only positive flows are incoming.
    idx_to = net.Edges_to_by_node{n};
    if ~isempty(idx_to)
        q_to = get_edge_flows(net, idx_to, t);
        mask = q_to > 0;
        Qin  = Qin + sum(q_to(mask));
        num  = num  + sum(q_to(mask) .* Tout_edge(idx_to(mask)));
    end

    % Edges that start at node n (n -> j). Negative flows mean "into n".
    idx_from = net.Edges_from_by_node{n};
    if ~isempty(idx_from)
        q_from = get_edge_flows(net, idx_from, t);
        mask   = q_from < 0;
        Qin    = Qin + sum(-q_from(mask));
        num    = num  + sum((-q_from(mask)) .* Tout_edge(idx_from(mask)));
    end

    flow_in(n)  = Qin;                                        % total incoming mass flow
    Tnode_in(n) = (Qin > 0) * (num/max(Qin,eps)) + ...        % weighted average
                  (Qin == 0) * Text;                          % no inflow → ambient
end

%% 3) Node components → outlet temperatures
% Take the node inlet temperature and apply the component's physics to get
% the node outlet temperature (what the pipe connected downstream sees).
Tnode_out = Tnode_in;
for n = 1:N
    comp = lower(net.Nodes(n).component);
    prm  = net.Nodes(n).params;

    switch comp
        case 'source'
            % Default: if no model given or 'dirichlet' → fixed setpoint
            if ~isfield(prm,'model') || strcmpi(prm.model,'dirichlet')
                Tnode_out(n) = Tin_source(t);
            else
                % Dynamic source: uses its own state to set outlet temp
                off = net.Nodes(n).offset;
                if off <= 0
                    % Safety fallback: if no state allocated, behave like Dirichlet
                    Tnode_out(n) = Tin_source(t);
                else
                    Tmix = Tnode_in(n);    % water-side inlet into source
                    Tint = z(off);         % source internal average temperature (state)
                    % theta1 weighting; clamp to safe range
                    if isfield(prm,'theta1'), theta1 = prm.theta1; else, theta1 = 0.7; end
                    theta1 = min(max(theta1,1e-6), 0.999999);
                    % T = (1-theta1)*Tin + theta1*Tout  =>  Tout = (T - (1-theta1)*Tin)/theta1
                    Tnode_out(n) = (Tint - (1-theta1)*Tmix)/theta1;
                end
            end

        case 'house'
            % A house adds/removes heat depending on the mode and its profile
            qloc = max(flow_in(n), 1e-12); % avoid zero division
            Tnode_out(n) = prosumer_house(t, Tnode_in(n), qloc, cp, prm, mode);

        case 'dry_cooler'
            % Dry cooler pushes towards ambient with a minimum approach
            qloc = max(flow_in(n), 1e-12);
            Tnode_out(n) = dry_cooler(Tnode_in(n), qloc, cp, prm.UA, Text, prm.approach_min);

        case {'tes','ates','btes'}
            % Storage nodes output their stored temperature/mean temperature
            off = net.Nodes(n).offset;
            if off > 0
                if strcmp(comp,'tes')
                    Tnode_out(n) = z(off);      % TES stores temperature directly
                else
                    C = 1; if isfield(prm,'C'), C = prm.C; end
                    Tnode_out(n) = z(off)/C;   % ATES/BTES may store energy (C*T)
                end
            end

        otherwise
            % Pass-through node: outlet = inlet
    end
end

% THEORY (pipe advection): upwind transport along the flow plus heat loss to ambient
%% 4) Pipe dynamics (upwind advection + ambient exchange)
% For each pipe, compute the time derivative for all its cells:
% - advection: shift temperature along the flow direction
% - ambient loss: relax toward ambient Text with rate alpha
dz = zeros(size(z));
for e = 1:E
    off    = net.Edges(e).offset;
    ns     = net.Edges(e).n;
    Tcells = z(off:off+ns-1);

    q   = get_edge_flows(net, e, t);          % mass flow [kg/s]
    Acs = max(net.Edges(e).Acs, 1e-6);        % cross-sectional area [m^2]
    v   = (q / net.rho) / Acs;                % bulk velocity [m/s]

    % Package parameters for the pipe RHS function
    p.n     = ns;
    p.dx    = net.Edges(e).L / ns;
    p.v     = v;
    p.alpha = net.Edges(e).alpha;
    p.Text  = Text;

    % Upstream boundary condition depends on the sign of the flow
    Tin = Tnode_out( (q >= 0) * net.Edges(e).from + (q < 0) * net.Edges(e).to );

    % First-order upwind; direction set by sign(q)
    dT = pipe_rhs_bidir(t, Tcells, p, Tin);
    dz(off:off+ns-1) = dT;
end

%% 5) Dynamic node state ODEs (TES/ATES/BTES/Source)
% Some nodes have their own internal state (one number in the big vector).
% Here we write how that state changes over time (its ODE).
for n = 1:N
    off = net.Nodes(n).offset;
    if off > 0
        comp = lower(net.Nodes(n).component);
        prm  = net.Nodes(n).params;

        switch comp
            case 'tes'
                % Cold-side TES with a simple energy balance model
                qloc   = max(flow_in(n), 1e-12);
                inputs = struct('mdot', qloc, 'Tin', Tnode_in(n));
                ptes   = struct('cp',cp,'Text',Text,'C',prm.C,'UA',prm.UA);
                dz(off) = tes_cold(t, z(off), inputs, ptes);

            case {'ates','btes'}
                % Generic storage: d/dt(C*T) = mdot*cp*(Tin - Tmean) - UA*(Tmean - Text)
                qloc = max(flow_in(n), 1e-12);
                C  = 1; if isfield(prm,'C'),  C = prm.C;  end
                UA = 0; if isfield(prm,'UA'), UA = prm.UA; end
                Tmean = z(off)/C;
                dz(off) = qloc*cp*(Tnode_in(n) - Tmean) - UA*(Tmean - Text);

            case 'source'
                % Dynamic source internal energy balance
                if isfield(prm,'model') && strcmpi(prm.model,'dynamic')
                    % Safe defaults if some parameters are missing
                    if isfield(prm,'rho'),    rho    = prm.rho;   else, rho    = net.rho; end
                    if isfield(prm,'cp'),     cp_loc = prm.cp;    else, cp_loc = cp;      end
                    if isfield(prm,'V'),      V      = prm.V;     else, V      = 0.05;    end
                    if isfield(prm,'kQ'),     kQ     = prm.kQ;    else, kQ     = 4.0e3;   end

                    Tin_w = Tnode_in(n);     % water-side inlet
                    T_int = z(off);          % internal average water temperature (state)

                    % Driver temperature (refrigerant side)
                    if isfield(prm,'Tc_fun') && ~isempty(prm.Tc_fun)
                        Tc = prm.Tc_fun(t);
                    else
                        Tc = Tin_source(t);
                    end

                    qloc = max(flow_in(n), 1e-12);       % mass flow through source [kg/s]
                    Fin  = qloc / max(rho,1e-9);         % volumetric inflow [m^3/s]
                    denom = max(rho*cp_loc*V, 1e-9);     % avoid division by zero

                    % THEORY (source): first order energy balance, outlet follows the command T_0s with a lag (this is T_F0)
                    % Energy balance: rho*cp*V * dT/dt = rho*cp*Fin*(Tin - T) + kQ*(Tc - T)
                    dTdt = (rho*cp_loc*Fin*(Tin_w - T_int) + kQ*(Tc - T_int)) / denom;
                    dz(off) = dTdt;
                end
        end
    end
end
end
