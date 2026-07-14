function res = air_baseline_run(net, res_wu, p, scenario_start, T_warm, T_sim, ei, F0_idx, R0_idx, con_idx, Np, solve_fn)
% AIR_BASELINE_RUN  Closed-loop runner for the plant-based baselines. Mirrors the
% production baseline_cl_run EXACTLY (same warmup, demand forecast at h*Ts, one-Ts
% plant apply of T0s + r_q, recorded outputs) so the only difference between KPC and
% a baseline is the optimiser inside solve_fn. solve_fn must return a 4th output
% (solver iterations), which we log alongside solve time and exitflag.
%   [u_opt, solve_ms, exitflag, iterations] = solve_fn(z_state, q0, d_seq, u_prev, t_base)
Ts = p.Ts; N_cl = round(T_sim/Ts); N_wu = numel(res_wu.t);
n_edges = numel(net.Edges); n_user = ei.n_user; N_tot = N_wu + N_cl;
traj.t = zeros(1,N_tot); traj.q_edges = zeros(n_edges,N_tot); traj.q_users = zeros(n_user,N_tot);
traj.r_q = zeros(n_edges,N_tot); traj.T_0s = zeros(1,N_tot); traj.T_0r = zeros(1,N_tot);
traj.T_F0 = zeros(1,N_tot); traj.T_is = zeros(n_user,N_tot); traj.T_ir = zeros(n_user,N_tot);
traj.d = zeros(n_user,N_tot); traj.c = zeros(n_user,N_tot);
iw = 1:N_wu;
traj.t(iw)=res_wu.t; traj.q_edges(:,iw)=res_wu.q_edges; traj.q_users(:,iw)=res_wu.q_users;
traj.r_q(:,iw)=res_wu.r_q; traj.T_0s(iw)=p.Tin_nom+res_wu.u; traj.T_0r(iw)=res_wu.Tout(R0_idx,:);
traj.T_F0(iw)=res_wu.Tout(F0_idx,:); traj.T_is(:,iw)=res_wu.T_s_i; traj.T_ir(:,iw)=res_wu.T_r_i;
traj.d(:,iw)=res_wu.d_i; traj.c(:,iw)=res_wu.c_i;
z_state = res_wu.z_final; u_prev = [traj.T_0s(N_wu); traj.r_q(:,N_wu)];
solve_ms = zeros(1,N_cl); exitflag = zeros(1,N_cl); iters_k = zeros(1,N_cl);
% per-step solve inputs, logged so a robustness probe can re-solve any step exactly
nz = numel(z_state);
z_steps  = zeros(nz, N_cl); q0_steps = zeros(n_edges, N_cl);
up_steps = zeros(1+n_edges, N_cl); tb_steps = zeros(1, N_cl);
for k = 1:N_cl
    cur = N_wu + k - 1; t_now = traj.t(cur);
    d_seq = zeros(n_user, Np);
    for h = 1:Np
        for i = 1:n_user
            d_seq(i,h) = compute_prosumer_Q(scenario_start + t_now + h*Ts, net.Nodes(con_idx(i)).params);
        end
    end
    z_steps(:,k)=z_state; q0_steps(:,k)=traj.q_edges(:,cur); up_steps(:,k)=u_prev; tb_steps(k)=scenario_start+t_now;
    [u_opt, sms, ef, it] = solve_fn(z_state, traj.q_edges(:,cur), d_seq, u_prev, scenario_start + t_now);
    solve_ms(k)=sms; exitflag(k)=ef; iters_k(k)=it;
    T_0s_apply = u_opt(1); r_q_apply = u_opt(2:1+n_edges);
    p_step = p; p_step.r_q_fun = @(t) r_q_apply; p_step.t_offset = scenario_start + t_now;
    net_step = net; net_step.q0 = traj.q_edges(:,cur);
    res_step = simulate_plant(net_step, z_state, p_step, @(t) T_0s_apply - p.Tin_nom, @(t) 1.0, Ts);
    nx = cur + 1;
    traj.t(nx)=t_now+Ts; traj.q_edges(:,nx)=res_step.q_edges(:,end); traj.q_users(:,nx)=res_step.q_users(:,end);
    traj.r_q(:,nx)=res_step.r_q(:,end); traj.T_0s(nx)=p.Tin_nom+res_step.u(end); traj.T_0r(nx)=res_step.Tout(R0_idx,end);
    traj.T_F0(nx)=res_step.Tout(F0_idx,end); traj.T_is(:,nx)=res_step.T_s_i(:,end); traj.T_ir(:,nx)=res_step.T_r_i(:,end);
    traj.d(:,nx)=res_step.d_i(:,end); traj.c(:,nx)=res_step.c_i(:,end);
    z_state = res_step.z_final; u_prev = [T_0s_apply; r_q_apply];
end
ci = N_wu + (1:N_cl);
res.d_i=traj.d(:,ci); res.c_i=traj.c(:,ci); res.r_q=traj.r_q(:,ci); res.T_0s=traj.T_0s(ci);
res.T_0r=traj.T_0r(ci); res.T_F0=traj.T_F0(ci); res.T_is=traj.T_is(:,ci); res.T_ir=traj.T_ir(:,ci);
res.q_users=traj.q_users(:,ci); res.Q_net=sum(traj.q_users(:,ci),1);
res.solve_ms=solve_ms; res.exitflag=exitflag; res.iters=iters_k;
% per-step solve inputs, exported for offline analysis
res.z_steps=z_steps; res.q0_steps=q0_steps; res.up_steps=up_steps; res.tb_steps=tb_steps;
end
