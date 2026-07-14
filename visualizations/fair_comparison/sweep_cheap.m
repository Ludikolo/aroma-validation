% SWEEP_CHEAP  Per-controller sensitivity for the CHEAP controllers (no NMPC):
%   (fairness check: horizon) Np sweep for KPC / Koopman-LMPC / Jacobian-LMPC -> are the baselines
%            disadvantaged by sharing KPC's tuned Np=12?
%   (fairness check: FD step) finite-difference step sweep for Jacobian-LMPC -> is worst-met
%            insensitive to the FD step?
% NMPC is excluded (a per-Np NMPC sweep is ~2.5h/value; done separately/coarsely).
% Uses the SAME equalised cost mode + shared box as the main run (cost_mode arg).
%
% RUN: matlab -batch "SW_COST='ls'; run('<abs>/sweep_cheap.m')"
clearvars -except SW_COST; clc;
here = fileparts(mfilename('fullpath')); root = fileparts(fileparts(here));
addpath(root); addpath(here); addpath(fullfile(root,'controller_closed_loop','lmpc')); startup();
if ~exist('SW_COST','var')||isempty(SW_COST), SW_COST='ls'; end
p=params(); mdot_scale=0.35; sc_start=0; T_warm=30*60; T_sim=24*3600;
pred_dir=fullfile(root,'predictor_open_loop','results'); ctrl_dir=fullfile(root,'controller_closed_loop','results');
F=load(fullfile(pred_dir,'vseq_fits_full.mat')); Hb=load(fullfile(ctrl_dir,'best_horizons.mat'));
Ab=load(fullfile(ctrl_dir,'best_alpha.mat')); Bb=load(fullfile(ctrl_dir,'best_tune.mat'));
S=load(fullfile(pred_dir,'iterated_AB.mat')); if isfield(S,'fit_iter'),fit_iter=S.fit_iter; else,fit_iter=S; end
tune0=Bb.best_tune; tune0.dT_0s_max=7.5; tune0.dr_q_max=4.5; tune0.dT_ir_max=15.0;
tune0.use_kirchhoff=true; tune0.use_Tr_le_Ts=true; tune0.use_mixing=true; tune0.rho_slack_mix=1;
tune0.Nc=Hb.Nc_best; tune0.alpha_energy=Ab.alpha_star;
[net,z0_cold]=build_plant(p);
net.mdotEdges=mdot_scale*net.mdotEdges; net.flows=net.mdotEdges; net.q0=net.mdotEdges(:);
ei=edge_user_index(net); n_user=ei.n_user; n_edges=numel(net.Edges);
F0_idx=find(strcmp({net.Nodes.name},'F0')); R0_idx=find(strcmp({net.Nodes.name},'R0'));
con_idx=arrayfun(@(i) find(strcmp({net.Nodes.name},ei.consumers{i})),1:n_user);
p_wu=p; p_wu.t_offset=sc_start; p_wu.r_q_fun=@(t) net.mdotEdges(:);
res_wu=simulate_plant(net,z0_cold,p_wu,@(t)0,@(t)1.0,T_warm);
K=build_incidence_v2(net);
basecfg.net=net; basecfg.p=p; basecfg.Ts=p.Ts; basecfg.n_user=n_user;
basecfg.md=net.mdotEdges(:); basecfg.N=null([K.M_supply;K.M_return]); basecfg.R0_idx=R0_idx;
basecfg.T0s_lo=max(p.Tin_nom+p.u_min,tune0.T_0s_min_floor); basecfg.T0s_hi=p.Tin_nom+p.u_max;
basecfg.dT0s=tune0.dT_0s_max; basecfg.drq=tune0.dr_q_max; basecfg.rq_lo_factor=tune0.r_q_lo_factor;
basecfg.user=ei.user; basecfg.rqhi_usr=p.excite.r_q_hi_factor*basecfg.md(ei.user);
basecfg.adeq=tune0.adequacy_safety; basecfg.dT_floor=p.Tin_nom-p.consumer.T_r_min;
basecfg.cost_mode=SW_COST; basecfg.alpha_energy=tune0.alpha_energy;
basecfg.R_du_T0s=tune0.R_du_T0s; basecfg.R_du_rq=tune0.R_du_rq; basecfg.cost_scale=1e3;
basecfg.fd_T0s=2.0; basecfg.fd_flow=0.05;
basecfg.qopts=optimoptions('quadprog','Display','off','OptimalityTolerance',1e-9,...
    'ConstraintTolerance',1e-9,'StepTolerance',1e-12,'MaxIterations',200);
worst=@(r) min(100*sum(r.c_i,2)./max(sum(r.d_i,2),1e-9));

fprintf('=== Np sweep, cost=%s (KPC / Koopman-LMPC / Jacobian-LMPC) ===\n', SW_COST);
NPs=[8 12 16];
swp=struct('Np',{},'kpc',{},'klmpc',{},'jlmpc',{});
for q=1:numel(NPs)
  Np=NPs(q);
  pred=truncate_pred_to_Np(build_kpc_v2_matrices(F.fits),Np);
  tune=tune0; tune_l=tune; tune_l.use_mixing=false;
  cfg=basecfg; cfg.Np=Np;
  rk=kpc_step_loop(net,res_wu,p,sc_start,T_warm,T_sim,pred,tune,ei,F0_idx,R0_idx,con_idx);
  rl=lmpc_step_loop(net,res_wu,p,sc_start,T_warm,T_sim,fit_iter,Np,tune_l,ei,F0_idx,R0_idx,con_idx);
  rj=air_baseline_run(net,res_wu,p,sc_start,T_warm,T_sim,ei,F0_idx,R0_idx,con_idx,Np,...
        @(z,qq,d,up,tb) air_jlmpc_solve(z,qq,d,up,tb,cfg));
  swp(q)=struct('Np',Np,'kpc',worst(rk),'klmpc',worst(rl),'jlmpc',worst(rj));
  fprintf('  Np=%2d : KPC=%.3f  Koopman-LMPC=%.3f  Jacobian-LMPC=%.3f\n',Np,worst(rk),worst(rl),worst(rj));
end

fprintf('\n=== Jacobian-LMPC FD-step sweep at Np=12, cost=%s (T0s step; flow fixed 0.05) ===\n', SW_COST);
Np=12; cfg=basecfg; cfg.Np=Np;
fds=[0.5 2.0]; fdsw=[];   % 2.0 = the main-run value, 0.5 checks a 4x smaller step
for a=1:numel(fds)
    c2=cfg; c2.fd_T0s=fds(a); c2.fd_flow=0.05;
    rj=air_baseline_run(net,res_wu,p,sc_start,T_warm,T_sim,ei,F0_idx,R0_idx,con_idx,Np,...
        @(z,qq,d,up,tb) air_jlmpc_solve(z,qq,d,up,tb,c2));
    fdsw(end+1,:)=[fds(a) 0.05 worst(rj) sum(rj.exitflag>0)]; %#ok<SAGROW>
    fprintf('  fd_T0s=%.2f fd_flow=0.05 : Jacobian-LMPC worst=%.3f  conv=%d/96\n',fds(a),worst(rj),sum(rj.exitflag>0));
end
save(fullfile(here,sprintf('sweep_cheap_%s.mat',SW_COST)),'swp','fdsw','NPs','SW_COST');
fprintf('\nSaved sweep_cheap_%s.mat\n',SW_COST);
