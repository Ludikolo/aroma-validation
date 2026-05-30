function [r_q_fun, T0s_fun, info] = excitation_multi(net, p, T_sim, seed_base)
% EXCITATION_MULTI  Multi-channel excitation for the v2 plant.
%
% Per-edge PRBS on r_q in [lo*mdot, hi*mdot] with dwell times cycled
% across channels so the data matrix stays well conditioned. T_0s gets
% its own PRBS on the dTc offset.

if nargin < 4 || isempty(seed_base), seed_base = 0; end

n_edges = numel(net.Edges);
mdot_nom = net.mdotEdges(:);
t = (0:p.Ts:T_sim).';
N = numel(t);

% pull excitation config from params
lo_fac    = p.excite.r_q_lo_factor;
hi_fac    = p.excite.r_q_hi_factor;
dwells    = p.excite.r_q_dwells_s;
T0s_amp   = p.excite.T0s_amp;
T0s_dwell = p.excite.T0s_dwell_s;

% per-edge PRBS for r_q. Cycle through the dwell list so neighbouring
% channels do not switch on the same grid -> better conditioning.
R_q = zeros(n_edges, N);
for j = 1:n_edges
    dwell_j = dwells(mod(j-1, numel(dwells)) + 1);
    R_q(j, :) = excitation('prbs', t, struct( ...
        'lo',    lo_fac * mdot_nom(j), ...
        'hi',    hi_fac * mdot_nom(j), ...
        'dwell', dwell_j, ...
        'seed',  seed_base + 100 + j));
end

% T^(0,s) excitation as dTc offset, kept on its own seed line
dTc = excitation('prbs', t, struct( ...
    'amp',   T0s_amp, ...
    'dwell', T0s_dwell, ...
    'seed',  seed_base + 1));
dTc = dTc(:)';

% nearest-step lookup: simulate_plant calls these once per sample so the
% indexing cost is negligible
Ts = p.Ts;
r_q_fun = @(tt) R_q(:, sample_idx(tt, Ts, N));
T0s_fun = @(tt) dTc(sample_idx(tt, Ts, N));

info.t   = t;
info.R_q = R_q;
info.dTc = dTc;

end


function k = sample_idx(t, Ts, N)
% clamp to [1, N] so callers can query slightly past T_sim safely
k = min(N, max(1, floor(t / Ts) + 1));
end
