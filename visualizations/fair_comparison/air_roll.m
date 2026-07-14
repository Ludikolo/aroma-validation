function [C, T0r_seq, Qnet_seq] = air_roll(x, z_state, q0, t_base, cfg)
% AIR_ROLL  Roll the EXACT plant forward Np steps with the constant input
% (T0s, r_q) implied by x = [T0s; null-space flow coeffs], r_q = md + N*x(2:end).
% Returns delivered heat C (n_user x Np) and, for the energy term, the per-horizon
% source return temperature T0r_seq (1 x Np) and total network flow Qnet_seq (1 x Np).
% This is the baselines' perfect-model predictor (oracle): exact plant, exact state.
T0s = x(1); rq = cfg.md + cfg.N * x(2:end);
C = zeros(cfg.n_user, cfg.Np); T0r_seq = zeros(1, cfg.Np); Qnet_seq = zeros(1, cfg.Np);
z = z_state; qc = q0;
for h = 1:cfg.Np
    ps = cfg.p; ps.r_q_fun = @(t) rq; ps.t_offset = t_base + (h-1)*cfg.Ts;
    nh = cfg.net; nh.q0 = qc;
    r = simulate_plant(nh, z, ps, @(t) T0s - cfg.p.Tin_nom, @(t) 1.0, cfg.Ts);
    C(:,h)      = r.c_i(:,end);
    T0r_seq(h)  = r.Tout(cfg.R0_idx, end);
    Qnet_seq(h) = sum(r.q_users(:,end));
    z = r.z_final; qc = r.q_edges(:,end);
end
end
