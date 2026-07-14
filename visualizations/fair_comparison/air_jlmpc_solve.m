function [u_opt, sms, ef, it] = air_jlmpc_solve(z_state, q0, d_seq, u_prev, t_base, cfg)
% AIR_JLMPC_SOLVE  Linearise-the-exact-plant LMPC step: one convex QP.
% cfg.cost_mode 'ls' (the primary run) uses Gauss-Newton least-squares tracking;
% 'kpc' (the equalised-cost sensitivity) uses KPC's objective with delivery+energy
% LINEARISED via the finite-difference rollout and move-suppression kept as an
% EXACT convex quadratic in x.
t0 = tic;
x0 = [u_prev(1); cfg.N' * (u_prev(2:end) - cfg.md)]; nx = numel(x0);
step = [cfg.fd_T0s; cfg.fd_flow * ones(nx-1,1)];
[A, bcon, lb, ub] = air_build_cons(u_prev, cfg, d_seq);

if strcmp(cfg.cost_mode, 'ls')
    C0 = roll_only(x0, z_state, q0, t_base, cfg)/1e3; c0 = C0(:); dvec = d_seq(:)/1e3;
    Jac = zeros(numel(c0), nx);
    for j = 1:nx
        dx = zeros(nx,1); dx(j) = step(j);
        Cp = roll_only(x0 + dx, z_state, q0, t_base, cfg)/1e3;
        Jac(:,j) = (Cp(:) - c0) / step(j);
    end
    H = 2*(Jac'*Jac) + 1e-6*eye(nx); f = 2*Jac'*(c0 - Jac*x0 - dvec);
else
    % delivery + energy (the rollout part of J), linearised about x0 by FD
    [C0, T0r0, Q0] = air_roll(x0, z_state, q0, t_base, cfg);
    de0 = (-sum(C0(:)) + cfg.alpha_energy * sum((x0(1) - T0r0) .* Q0)) / cfg.cost_scale;
    g_de = zeros(nx, 1);
    for j = 1:nx
        dx = zeros(nx,1); dx(j) = step(j);
        [Cp, T0rp, Qp] = air_roll(x0 + dx, z_state, q0, t_base, cfg);
        dep = (-sum(Cp(:)) + cfg.alpha_energy * sum((x0(1) + dx(1) - T0rp) .* Qp)) / cfg.cost_scale;
        g_de(j) = (dep - de0) / step(j);
    end
    % move-suppression: exact convex quadratic in actual-x coordinates.
    % J_move/scale = [R_du_T0s (x1-up1)^2 + R_du_rq ||N v + (md - rqp)||^2]/scale
    Nmat = cfg.N; md = cfg.md; rqp = u_prev(2:end);
    Hm = zeros(nx); fm = zeros(nx,1);
    Hm(1,1) = 2*cfg.R_du_T0s;  fm(1) = -2*cfg.R_du_T0s*u_prev(1);
    Hm(2:end,2:end) = 2*cfg.R_du_rq*(Nmat'*Nmat);
    fm(2:end) = 2*cfg.R_du_rq*(Nmat'*(md - rqp));
    Hm = Hm / cfg.cost_scale; fm = fm / cfg.cost_scale;
    % total QP: 0.5 x'H x + f'x ; H = Hm (move), f = fm + g_de (linearised delivery+energy)
    H = Hm + 1e-9*eye(nx);
    f = fm + g_de;
end
[xo, ~, ef, output] = quadprog(H, f, A, bcon, [], [], lb, ub, x0, cfg.qopts);
sms = toc(t0) * 1000; it = output.iterations;
if ef <= 0      % hold u_prev on solver failure (same fallback as KPC); convex QP rarely needs it
    xo = x0;
end
u_opt = [xo(1); cfg.md + cfg.N * xo(2:end)];
end

function C = roll_only(x, z_state, q0, t_base, cfg)
[C, ~, ~] = air_roll(x, z_state, q0, t_base, cfg);
end
