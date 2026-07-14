function [A, bcon, lb, ub] = air_build_cons(u_prev, cfg, d_seq)
% AIR_BUILD_CONS  Shared HARD constraints for the plant-based baselines, made
% identical to KPC: flow box [rq_lo_factor*md, 1.5*md], per-step rate caps on
% T0s and flow, T0s box [T0s_lo, T0s_hi], demand-adequacy floor, and Kirchhoff
% (satisfied structurally because r_q = md + N*x with N a basis of null[M]).
% Depends only on (u_prev, d_seq), not the plant state.
%
% Adequacy floor: KPC imposes the floor PER HORIZON on r_q^(i)_h, and with
% move-blocking Nc=1 the flow is constant over the horizon, so KPC's binding floor
% is the MAX over the horizon. The baselines also apply one constant flow over the
% horizon (roll), so to match KPC's effective floor exactly we take the max over h:
%   floor^(i) = max_h min( adeq*d^(i)_h/(cp*dT_floor), 0.95*r_q_hi^(i) ).
nxi = size(cfg.N, 2);
P = [zeros(numel(cfg.md),1), cfg.N];        % r_q = md + P*x
rqp = u_prev(2:end); md = cfg.md; e1 = [1, zeros(1,nxi)];
% P*x <= 0.5 md => r_q <= 1.5 md (1.5 = p.excite.r_q_hi_factor, KPC's upper box);
% -P*x <= (1-rq_lo_factor) md => r_q >= rq_lo_factor md
A = [ P; -P; P; -P; e1; -e1 ];
bcon = [ 0.5*md; (1-cfg.rq_lo_factor)*md; cfg.drq + (rqp - md); cfg.drq - (rqp - md); ...
         u_prev(1) + cfg.dT0s; -(u_prev(1) - cfg.dT0s) ];
floor_h   = min(cfg.adeq * d_seq ./ (cfg.p.cp * cfg.dT_floor), 0.95 * cfg.rqhi_usr);  % n_user x Np
floor_val = max(floor_h, [], 2);                                                      % n_user x 1
Pu = P(cfg.user, :);
A = [A; -Pu]; bcon = [bcon; md(cfg.user) - floor_val];
lb = [cfg.T0s_lo; -inf(nxi,1)]; ub = [cfg.T0s_hi; inf(nxi,1)];
end
