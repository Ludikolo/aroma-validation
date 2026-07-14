function [u_opt, sms, ef, it] = air_nmpc_solve(z_state, q0, d_seq, u_prev, t_base, cfg)
% AIR_NMPC_SOLVE  Nonconvex NMPC step: minimise the TRUE objective on the exact
% nonlinear-plant rollout (fmincon SQP, exact-plant predictor). cfg.cost_mode
% selects the neutral symmetric least-squares tracking cost ('ls', the primary
% run) or KPC's production objective ('kpc', the equalised-cost sensitivity).
x0 = [u_prev(1); cfg.N' * (u_prev(2:end) - cfg.md)];
[A, bcon, lb, ub] = air_build_cons(u_prev, cfg, d_seq);
if strcmp(cfg.cost_mode, 'ls')
    obj = @(x) guard(@() sum((roll_only(x, z_state, q0, t_base, cfg)/1e3 - d_seq/1e3).^2, 'all'));
else
    obj = @(x) guard(@() obj_kpc(x, z_state, q0, t_base, u_prev, cfg));
end
t0 = tic;
% guard the whole solve too: if fmincon still cannot proceed (e.g. the objective is
% undefined even at the start), report a failed step (hold u_prev) rather than crash
% the run. This happens only for the plant baselines under the aggressive delivery-max
% cost, where the exact-plant rollout can diverge at points fmincon probes, a
% nonconvex-NLP fragility counted as a non-converged step.
try
    [xo, ~, ef, output] = fmincon(obj, x0, A, bcon, [], [], lb, ub, [], cfg.fopts);
    it = output.iterations;
catch
    xo = x0; ef = -98; it = 0;
end
% Receding-horizon fallback, identical to KPC's QP (kpc_v2_solve: on exitflag<=0 it
% holds the warm-start = previous applied input). On a solver failure fmincon may
% return an INFEASIBLE last iterate; applying it would drive the plant to garbage.
% Holding u_prev makes a failed step under-deliver rather than apply an infeasible
% iterate, the same fallback as KPC.
if ef <= 0
    xo = x0;
end
sms = toc(t0) * 1000;
u_opt = [xo(1); cfg.md + cfg.N * xo(2:end)];
end

function J = guard(fn)
% return a large finite penalty if the rollout/objective is non-finite, so fmincon
% steers away from that region instead of erroring on an undefined objective.
try
    J = fn();
    if ~isfinite(J), J = 1e6; end
catch
    J = 1e6;
end
end

function J = obj_kpc(x, z_state, q0, t_base, u_prev, cfg)
[C, T0r_seq, Qnet_seq] = air_roll(x, z_state, q0, t_base, cfg);
J = air_kpc_cost(x, C, T0r_seq, Qnet_seq, u_prev, cfg);
end

function C = roll_only(x, z_state, q0, t_base, cfg)
[C, ~, ~] = air_roll(x, z_state, q0, t_base, cfg);
end
