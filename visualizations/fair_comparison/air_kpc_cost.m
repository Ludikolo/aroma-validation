function [J, parts] = air_kpc_cost(x, C, T0r_seq, Qnet_seq, u_prev, cfg)
% AIR_KPC_COST  KPC's production objective evaluated on a rolled trajectory,
% scaled by 1/cost_scale (a constant, so the argmin and relative weights are
% unchanged; only the numeric magnitude is brought to O(1e3) for fmincon).
%
%   J = [ -sum_i,h c_ih                                  (delivery reward)
%         + alpha*sum_h (T0s - T0r_h)*Qnet_h             (source energy, cp folded into alpha)
%         + R_du_T0s*(T0s-T0s_prev)^2 + R_du_rq*||r_q-r_q_prev||^2  (move suppression)
%       ] / cost_scale
%
% This matches kpc_v2_solve's f/H exactly (with Nc=1 the move term reduces to the
% first move u_0-u_prev; T_ir is absent because the baselines do not actuate it).
T0s = x(1); rq = cfg.md + cfg.N * x(2:end);
delivery = -sum(C(:));
energy   = cfg.alpha_energy * sum((T0s - T0r_seq) .* Qnet_seq);
move     = cfg.R_du_T0s * (T0s - u_prev(1))^2 + cfg.R_du_rq * sum((rq - u_prev(2:end)).^2);
J = (delivery + energy + move) / cfg.cost_scale;
if nargout > 1
    parts = struct('delivery', delivery, 'energy', energy, 'move', move, 'scale', cfg.cost_scale);
end
end
