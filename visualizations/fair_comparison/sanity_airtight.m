% SANITY_AIRTIGHT  Fast gate before the full (~3h) comparison run.
% (1) KPC + Koopman-LMPC full 24h reproduce production comparison.mat (97.561/97.735)
%     -> the harness setup (plant/warmup/cfg/loops) matches production.
% (2) At two representative plant states (warmup + ~mid-day), the new baseline
%     solvers converge AND the optimiser DECREASES the true KPC objective vs the
%     warm-start x0 -> the cost wiring + Jacobian linearisation are sound.
% (3) NMPC iterations stay under the raised cap; build_cons polytope is feasible.
clearvars; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(root); addpath(here); addpath(fullfile(root, 'controller_closed_loop', 'lmpc')); startup();

p = params(); mdot_scale = 0.35; sc_start = 0; T_warm = 30*60; T_sim = 24*3600;
pred_dir = fullfile(root, 'predictor_open_loop', 'results');
ctrl_dir = fullfile(root, 'controller_closed_loop', 'results');
F = load(fullfile(pred_dir,'vseq_fits_full.mat'));
Hb = load(fullfile(ctrl_dir,'best_horizons.mat'));
Ab = load(fullfile(ctrl_dir,'best_alpha.mat'));
Bb = load(fullfile(ctrl_dir,'best_tune.mat'));
Np = Hb.Np_best;
pred = truncate_pred_to_Np(build_kpc_v2_matrices(F.fits), Np);
tune = Bb.best_tune;
tune.dT_0s_max=7.5; tune.dr_q_max=4.5; tune.dT_ir_max=15.0;
tune.use_kirchhoff=true; tune.use_Tr_le_Ts=true; tune.use_mixing=true; tune.rho_slack_mix=1;
tune.Nc=Hb.Nc_best; tune.alpha_energy=Ab.alpha_star;
S = load(fullfile(pred_dir,'iterated_AB.mat'));
if isfield(S,'fit_iter'), fit_iter=S.fit_iter; else, fit_iter=S; end
tune_lmpc = tune; tune_lmpc.use_mixing=false;

[net,z0_cold]=build_plant(p);
net.mdotEdges=mdot_scale*net.mdotEdges; net.flows=net.mdotEdges; net.q0=net.mdotEdges(:);
ei=edge_user_index(net); n_user=ei.n_user; n_edges=numel(net.Edges);
F0_idx=find(strcmp({net.Nodes.name},'F0')); R0_idx=find(strcmp({net.Nodes.name},'R0'));
con_idx=arrayfun(@(i) find(strcmp({net.Nodes.name},ei.consumers{i})),1:n_user);
p_wu=p; p_wu.t_offset=sc_start; p_wu.r_q_fun=@(t) net.mdotEdges(:);
res_wu=simulate_plant(net,z0_cold,p_wu,@(t)0,@(t)1.0,T_warm);

K=build_incidence_v2(net);
cfg.Np=Np; cfg.net=net; cfg.p=p; cfg.Ts=p.Ts; cfg.n_user=n_user;
cfg.md=net.mdotEdges(:); cfg.N=null([K.M_supply;K.M_return]); cfg.R0_idx=R0_idx;
cfg.T0s_lo=max(p.Tin_nom+p.u_min, tune.T_0s_min_floor); cfg.T0s_hi=p.Tin_nom+p.u_max;
cfg.dT0s=tune.dT_0s_max; cfg.drq=tune.dr_q_max; cfg.rq_lo_factor=tune.r_q_lo_factor;
cfg.user=ei.user; cfg.rqhi_usr=p.excite.r_q_hi_factor*cfg.md(ei.user);
cfg.adeq=tune.adequacy_safety; cfg.dT_floor=p.Tin_nom-p.consumer.T_r_min;
cfg.cost_mode='kpc'; cfg.alpha_energy=tune.alpha_energy;
cfg.R_du_T0s=tune.R_du_T0s; cfg.R_du_rq=tune.R_du_rq; cfg.cost_scale=1e3;
cfg.fd_T0s=2.0; cfg.fd_flow=0.05;
cfg.fopts=optimoptions('fmincon','Algorithm','sqp','Display','off',...
    'MaxFunctionEvaluations',10000,'MaxIterations',1000,...
    'FiniteDifferenceStepSize',cfg.fd_flow,'OptimalityTolerance',1e-6,'StepTolerance',1e-8);
cfg.qopts=optimoptions('quadprog','Display','off','OptimalityTolerance',1e-9,...
    'ConstraintTolerance',1e-9,'StepTolerance',1e-12,'MaxIterations',200);

fprintf('=== (1) setup gate: KPC + Koopman-LMPC full 24h ===\n');
rk = kpc_step_loop(net,res_wu,p,sc_start,T_warm,T_sim,pred,tune,ei,F0_idx,R0_idx,con_idx);
rl = lmpc_step_loop(net,res_wu,p,sc_start,T_warm,T_sim,fit_iter,Np,tune_lmpc,ei,F0_idx,R0_idx,con_idx);
wk=min(100*sum(rk.c_i,2)./max(sum(rk.d_i,2),1e-9));
wl=min(100*sum(rl.c_i,2)./max(sum(rl.d_i,2),1e-9));
fprintf('  KPC worst=%.3f (expect 97.561)   Koopman-LMPC worst=%.3f (expect 97.735)\n', wk, wl);
assert(abs(wk-97.561)<0.05, 'KPC setup mismatch'); assert(abs(wl-97.735)<0.05, 'Koopman-LMPC setup mismatch');
fprintf('  PASS setup reproduces production.\n');

fprintf('\n=== (2)/(3) baseline cost wiring at 2 states ===\n');
% state A: warmup; state B: plant stepped forward 35 steps at constant nominal input
q0A=res_wu.q_edges(:,end); zA=res_wu.z_final; tA=res_wu.t(end); uA=[p.Tin_nom; net.mdotEdges(:)];
pB=p; pB.t_offset=sc_start; pB.r_q_fun=@(t) net.mdotEdges(:);
resB=simulate_plant(net,z0_cold,pB,@(t)0,@(t)1.0,T_warm+35*p.Ts);
zB=resB.z_final; q0B=resB.q_edges(:,end); tB=resB.t(end);
% u_prev at B = last applied (nominal): [20; md]
uB=[p.Tin_nom; net.mdotEdges(:)];
states={struct('z',zA,'q0',q0A,'up',uA,'tb',tA,'name','warmup'), ...
        struct('z',zB,'q0',q0B,'up',uB,'tb',tB,'name','step35')};
for s=1:numel(states)
  st=states{s};
  d_seq=zeros(n_user,Np);
  for h=1:Np, for i=1:n_user, d_seq(i,h)=compute_prosumer_Q(sc_start+st.tb+h*p.Ts, net.Nodes(con_idx(i)).params); end; end
  x0=[st.up(1); cfg.N'*(st.up(2:end)-cfg.md)];
  % feasibility, on the same polytope the solvers use (max-over-horizon floor)
  [A,b,lb,ub]=air_build_cons(st.up,cfg,d_seq);
  [~,~,eflp]=linprog(zeros(numel(x0),1),A,b,[],[],lb,ub,optimoptions('linprog','Display','off'));
  % true objective at x0
  [C0,T0r0,Q0]=air_roll(x0,st.z,st.q0,st.tb,cfg); J0=air_kpc_cost(x0,C0,T0r0,Q0,st.up,cfg); dlv0=sum(C0(:));
  % jlmpc
  [uj,smsj,efj,itj]=air_jlmpc_solve(st.z,st.q0,d_seq,st.up,st.tb,cfg);
  xj=[uj(1); cfg.N'*(uj(2:end)-cfg.md)];
  [Cj,T0rj,Qj]=air_roll(xj,st.z,st.q0,st.tb,cfg); Jj=air_kpc_cost(xj,Cj,T0rj,Qj,st.up,cfg); dlvj=sum(Cj(:));
  % nmpc
  [un,smsn,efn,itn]=air_nmpc_solve(st.z,st.q0,d_seq,st.up,st.tb,cfg);
  xn=[un(1); cfg.N'*(un(2:end)-cfg.md)];
  [Cn,T0rn,Qn]=air_roll(xn,st.z,st.q0,st.tb,cfg); Jn=air_kpc_cost(xn,Cn,T0rn,Qn,st.up,cfg); dlvn=sum(Cn(:));
  fprintf('--- state %s: feas(linprog ef)=%d ---\n', st.name, eflp);
  fprintf('  J0=%.4f  Jjlmpc=%.4f (ef%d,%dit,%.0fms,dlv %.0f->%.0f)  Jnmpc=%.4f (ef%d,%dit,%.0fms,dlv %.0f->%.0f)\n', ...
      J0, Jj,efj,itj,smsj,dlv0,dlvj, Jn,efn,itn,smsn,dlv0,dlvn);
  assert(eflp==1, 'polytope infeasible');
  assert(efj>0, 'jlmpc did not converge'); assert(efn~=0, 'nmpc hit iteration cap');
  assert(Jj<=J0+1e-6, 'jlmpc did not reduce true J'); assert(Jn<=J0+1e-6, 'nmpc did not reduce true J');
  assert(itn<1000, 'nmpc hit MaxIterations');
  fprintf('  PASS (both solvers converge + reduce the true KPC objective; cap not hit).\n');
end
fprintf('\nALL SANITY CHECKS PASSED.\n');
