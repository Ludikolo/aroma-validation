% RUN_COMPARISON_AIRTIGHT  Airtight-fair 4-controller benchmark.
%
% Extends the part-4 run_comparison.m so the four controllers differ ONLY in the
% prediction model and the optimiser class, by removing every remaining cost/box
% asymmetry that harness still had:
%
%   (equalisation: cost)  In cost mode 'kpc' the two plant-based baselines
%                   (Jacobian-LMPC, NMPC) minimise KPC's exact production objective
%                       J = -sum_i,h c_ih  +  alpha*sum_h (T0s - T0r_h)*Qnet_h
%                           +  R_du_T0s*(dT0s)^2 + R_du_rq*||drq||^2
%                   (delivery reward + source energy + move suppression), instead
%                   of the old symmetric least-squares sum((c-d)/1e3)^2. NMPC
%                   minimises the TRUE nonlinear J on the exact-plant rollout;
%                   Jacobian-LMPC linearises delivery+energy via the FD rollout and
%                   keeps move suppression as an exact convex quadratic. The whole
%                   objective is scaled by 1/cost_scale (a constant, no effect on the
%                   argmin or relative weights) so fmincon is well conditioned. KPC
%                   and Koopman-LMPC are UNCHANGED (they already use this objective).
%                   The default mode 'ls' instead keeps the baselines on the neutral
%                   least-squares tracking cost (conservative for KPC, see the
%                   README); the table is reported under both conventions.
%   (equalisation: flow box)   Baseline flow lower bound 0.3*md -> r_q_lo_factor*md (=0.5694).
%   (equalisation: temperature floor)   Baseline T0s lower floor 14C -> max(14, T_0s_min_floor)=17.69C.
%   (equalisation: NMPC stop rule)       NMPC runs to a tight SQP stop (StepTolerance 1e-8) under a
%                   generous bounded budget (MaxIterations 1000, MaxFunctionEvaluations
%                   10000; the shipped ls run uses at most 4 iterations); cap hits are
%                   counted and reported below.
%   (equalisation: solver logging)        NMPC + Jacobian-LMPC log solver iterations + exitflag.
%
% Intentionally KEPT and documented (NOT "fixed" in KPC's favour):
%   - Perfect-model advantage: NMPC + Jacobian-LMPC roll the EXACT plant and read
%     the EXACT internal state (an oracle predictor); KPC + Koopman-LMPC use the
%     learned lift. This HANDICAPS KPC.
%   - KPC's soft slacks on c<=d / Tr<=Ts / mixing are a feasibility mechanism for
%     the learned predictor's model error; the baselines have no such constraint
%     (the exact-plant rollout enforces c<=d physically). The shared HARD
%     constraints (box, rate caps, adequacy floor, Kirchhoff) are identical for all.
%   - T_ir is an internal QP decision in KPC/Koopman-LMPC, never applied to the
%     plant; the actuated DOF (T0s + flow null-space) is identical for all four.
%
% RUN: matlab -batch "run('<abs>/run_comparison_airtight.m')"
% Optional workspace overrides set before running:
%   AIR_NP        prediction horizon (default best_horizons.Np_best = 12)
%   AIR_COSTMODE  'ls' (default, the headline run) or 'kpc' (equalised-cost sensitivity)
%   AIR_FDT0S     Jacobian-LMPC T0s finite-difference step (default 2.0; see the FD-step sweep in sweep_cheap.m)
%   AIR_NMPC_STEPTOL  NMPC StepTolerance (default 1e-8; see the cfg note below)
%   AIR_TAG       suffix for the saved .mat (default '')

clearvars('-except', 'AIR_NP', 'AIR_COSTMODE', 'AIR_FDT0S', 'AIR_TAG', 'AIR_NMPC_STEPTOL'); clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));          % visualizations/fair_comparison -> root
addpath(root); addpath(here); addpath(fullfile(root, 'controller_closed_loop', 'lmpc')); startup();

mdot_scale = 0.35; sc_start = 0; T_warm = 30*60; T_sim = 24*3600;
p = params();
pred_dir = fullfile(root, 'predictor_open_loop', 'results');
ctrl_dir = fullfile(root, 'controller_closed_loop', 'results');

if exist('AIR_COSTMODE','var') && ~isempty(AIR_COSTMODE), cost_mode = AIR_COSTMODE; else, cost_mode = 'ls'; end   % 'ls' is the headline run of the thesis table
if exist('AIR_FDT0S','var')    && ~isempty(AIR_FDT0S),    fd_T0s    = AIR_FDT0S;    else, fd_T0s = 2.0; end
if exist('AIR_TAG','var')      && ~isempty(AIR_TAG),      tag       = AIR_TAG;      else, tag = ''; end

% --- KPC (direct multi-step) ---
F  = load(fullfile(pred_dir, 'vseq_fits_full.mat'));
Hb = load(fullfile(ctrl_dir, 'best_horizons.mat'));
Ab = load(fullfile(ctrl_dir, 'best_alpha.mat'));
Bb = load(fullfile(ctrl_dir, 'best_tune.mat'));
if exist('AIR_NP','var') && ~isempty(AIR_NP), Np = AIR_NP; else, Np = Hb.Np_best; end
pred = truncate_pred_to_Np(build_kpc_v2_matrices(F.fits), Np);
tune = Bb.best_tune;
tune.dT_0s_max = 7.5; tune.dr_q_max = 4.5; tune.dT_ir_max = 15.0;
tune.use_kirchhoff = true; tune.use_Tr_le_Ts = true; tune.use_mixing = true; tune.rho_slack_mix = 1;
tune.Nc = Hb.Nc_best; tune.alpha_energy = Ab.alpha_star;
fprintf('AIRTIGHT comparison: n_z=%d, Np=%d, mdot=%.2f, cost=%s, fd_T0s=%.2f\n', ...
        F.fits.n_z, Np, mdot_scale, cost_mode, fd_T0s);

% --- iterated one-step Koopman model for Koopman-LMPC ---
S = load(fullfile(pred_dir, 'iterated_AB.mat'));
if isfield(S, 'fit_iter'), fit_iter = S.fit_iter; else, fit_iter = S; end
tune_lmpc = tune; tune_lmpc.use_mixing = false;

% --- plant + shared warmup (identical for all four) ---
[net, z0_cold] = build_plant(p);
net.mdotEdges = mdot_scale * net.mdotEdges; net.flows = net.mdotEdges; net.q0 = net.mdotEdges(:);
ei = edge_user_index(net); n_user = ei.n_user; n_edges = numel(net.Edges);
F0_idx = find(strcmp({net.Nodes.name}, 'F0'));
R0_idx = find(strcmp({net.Nodes.name}, 'R0'));
con_idx = arrayfun(@(i) find(strcmp({net.Nodes.name}, ei.consumers{i})), 1:n_user);
p_wu = p; p_wu.t_offset = sc_start; p_wu.r_q_fun = @(t) net.mdotEdges(:);
res_wu = simulate_plant(net, z0_cold, p_wu, @(t) 0, @(t) 1.0, T_warm);

% --- shared config for the plant baselines (equalised to KPC) ---
K = build_incidence_v2(net);
cfg.Np = Np; cfg.net = net; cfg.p = p; cfg.Ts = p.Ts; cfg.n_user = n_user;
cfg.md = net.mdotEdges(:); cfg.N = null([K.M_supply; K.M_return]);
cfg.R0_idx = R0_idx;
cfg.T0s_lo = max(p.Tin_nom + p.u_min, tune.T_0s_min_floor);   % equalisation: temperature floor
cfg.T0s_hi = p.Tin_nom + p.u_max;
cfg.dT0s = tune.dT_0s_max; cfg.drq = tune.dr_q_max;
cfg.rq_lo_factor = tune.r_q_lo_factor;                        % equalisation: flow box
cfg.user = ei.user; cfg.rqhi_usr = p.excite.r_q_hi_factor * cfg.md(ei.user);
cfg.adeq = tune.adequacy_safety; cfg.dT_floor = p.Tin_nom - p.consumer.T_r_min;
cfg.cost_mode    = cost_mode;        % equalisation: cost
cfg.alpha_energy = tune.alpha_energy;
cfg.R_du_T0s     = tune.R_du_T0s;
cfg.R_du_rq      = tune.R_du_rq;
cfg.cost_scale   = 1e3;
cfg.fd_T0s       = fd_T0s; cfg.fd_flow = 0.05;                % equalisation: FD step
% NMPC with the settings of the shipped headline run: StepTolerance 1e-8, tighter than
% fmincon's 1e-6 default (the QP baselines run at 1e-9, which fmincon cannot reliably
% reach), under a 1000-iteration budget that is never binding (at most 4 iterations
% used). A rerun at the looser natural 1e-6 stop with a 50/500 cap moves the NMPC row
% by 0.44 pp and the convex rows not at all (airtight_comparison_ls_rerun.mat).
% Iterations + exitflag are logged and the cap-hit count is reported (equalisation: solver logging and stop rule).
cfg.nmpc_steptol = 1e-8; cfg.nmpc_maxiter = 1000; cfg.nmpc_maxfun = 10000;
if exist('AIR_NMPC_STEPTOL','var') && ~isempty(AIR_NMPC_STEPTOL), cfg.nmpc_steptol = AIR_NMPC_STEPTOL; end
cfg.fopts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
                         'MaxFunctionEvaluations', cfg.nmpc_maxfun, 'MaxIterations', cfg.nmpc_maxiter, ...
                         'FiniteDifferenceStepSize', cfg.fd_flow, ...
                         'OptimalityTolerance', 1e-6, 'StepTolerance', cfg.nmpc_steptol);
% same algorithm (interior-point-convex, the quadprog default) and tolerances as kpc_v2_solve
cfg.qopts = optimoptions('quadprog', 'Display', 'off', ...
                         'OptimalityTolerance', 1e-9, 'ConstraintTolerance', 1e-9, ...
                         'StepTolerance', 1e-12, 'MaxIterations', 200);

% --- run all four ---
fprintf('Running KPC ...\n'); t0 = tic;
res_kpc = kpc_step_loop(net, res_wu, p, sc_start, T_warm, T_sim, pred, tune, ei, F0_idx, R0_idx, con_idx);
fprintf('  KPC %.1f s wall\n', toc(t0));
fprintf('Running Koopman-LMPC ...\n'); t0 = tic;
res_klmpc = lmpc_step_loop(net, res_wu, p, sc_start, T_warm, T_sim, fit_iter, Np, tune_lmpc, ei, F0_idx, R0_idx, con_idx);
fprintf('  Koopman-LMPC %.1f s wall\n', toc(t0));
fprintf('Running Jacobian-LMPC ...\n'); t0 = tic;
res_jlmpc = air_baseline_run(net, res_wu, p, sc_start, T_warm, T_sim, ei, F0_idx, R0_idx, con_idx, Np, ...
                             @(z,q,d,up,tb) air_jlmpc_solve(z, q, d, up, tb, cfg));
fprintf('  Jacobian-LMPC %.1f s wall\n', toc(t0));
fprintf('Running NMPC (raised budget, to convergence) ...\n'); t0 = tic;
res_nmpc = air_baseline_run(net, res_wu, p, sc_start, T_warm, T_sim, ei, F0_idx, R0_idx, con_idx, Np, ...
                            @(z,q,d,up,tb) air_nmpc_solve(z, q, d, up, tb, cfg));
fprintf('  NMPC %.1f s wall\n', toc(t0));

% --- metrics ---
names = {'KPC', 'Koopman-LMPC', 'Jacobian-LMPC', 'NMPC'};
R = {res_kpc, res_klmpc, res_jlmpc, res_nmpc};
if strcmp(cost_mode,'kpc')
    costnote = 'identical cost (delivery+energy+move) + box+floor+rate+Kirchhoff';
else
    costnote = 'baselines on neutral LS cost (conservative); identical box+floor+rate+Kirchhoff';
end
fprintf('\nAIRTIGHT (mdot=%.2f, %.1fh, Np=%d, cost=%s): %s\n', ...
        mdot_scale, T_sim/3600, Np, cost_mode, costnote);
fprintf('%-14s  worst%%   unmet[kWh]  Edlv[kWh]  Esrc[kWh]  intens  med-it  max-it  med-ms  max-ms  conv\n', 'controller');
M = struct();
for j = 1:numel(R)
    r = R{j};
    pc  = 100 * sum(r.c_i, 2) ./ max(sum(r.d_i, 2), 1e-9);
    unm = sum(max(r.d_i - r.c_i, 0), 'all') * p.Ts / 3.6e6;
    Edl = sum(r.c_i, 'all') * p.Ts / 3.6e6;                         % delivered heat [kWh]
    Esr = p.cp * sum(r.Q_net .* (r.T_F0 - r.T_0r)) * p.Ts / 3.6e6;  % source energy [kWh]
    intens = Esr / max(Edl, 1e-9);                                  % matched-service intensity
    conv = sum(r.exitflag > 0);
    if isfield(r, 'iters'), mit = median(r.iters); xit = max(r.iters); else, mit = NaN; xit = NaN; end
    M.(matlab.lang.makeValidName(names{j})) = struct('percons', pc, 'worst', min(pc), ...
        'unmet', unm, 'Edlv', Edl, 'Esrc', Esr, 'intensity', intens, 'med_iters', mit, 'max_iters', xit, ...
        'med_ms', median(r.solve_ms), 'max_ms', max(r.solve_ms), 'converged', conv, ...
        'n', numel(r.exitflag), 'exitflag', r.exitflag(:)');
    fprintf('%-14s  %6.3f  %9.3f  %8.2f  %8.2f  %6.4f  %6.1f  %6.0f  %6.0f  %6.0f  %d/%d\n', ...
        names{j}, min(pc), unm, Edl, Esr, intens, mit, xit, median(r.solve_ms), max(r.solve_ms), conv, numel(r.exitflag));
end
fprintf('\nper-consumer met%%:\n');
for j = 1:numel(R)
    pc = 100 * sum(R{j}.c_i, 2) ./ max(sum(R{j}.d_i, 2), 1e-9);
    fprintf('  %-14s [%s]\n', names{j}, sprintf('%.3f ', pc));
end
ef = res_nmpc.exitflag; u = unique(ef);
fprintf('NMPC exitflag: '); for k = 1:numel(u), fprintf('ef%d:%d  ', u(k), sum(ef==u(k))); end; fprintf('\n');
% report cap behaviour; a cap-hit (exitflag 0) is an honest "did not converge in
% budget" finding, not a reason to discard the run, so report it rather than throw.
ncap = sum(res_nmpc.iters >= cfg.nmpc_maxiter | res_nmpc.exitflag == 0);
if ncap == 0
    fprintf('NMPC max iters used = %d (cap %d) -> cap never binds\n', max(res_nmpc.iters), cfg.nmpc_maxiter);
else
    fprintf('NMPC HIT the iteration cap on %d/%d steps (max iters %d, cap %d) -> reported, not optimal there\n', ...
        ncap, numel(res_nmpc.iters), max(res_nmpc.iters), cfg.nmpc_maxiter);
end
fprintf('NMPC converged-to-optimality (ef1) = %d/%d; stalled (ef2) = %d; failed (ef<=0) = %d\n', ...
    sum(res_nmpc.exitflag==1), numel(res_nmpc.exitflag), sum(res_nmpc.exitflag==2), sum(res_nmpc.exitflag<=0));

fname = sprintf('airtight_comparison%s.mat', local_suffix(tag));
save(fullfile(here, fname), 'res_kpc', 'res_klmpc', 'res_jlmpc', 'res_nmpc', ...
     'names', 'M', 'mdot_scale', 'Np', 'cost_mode', 'fd_T0s', 'p', 'ei', ...
     'cfg', 'con_idx', 'sc_start', '-v7.3');
fprintf('\nSaved %s\n', fname);

function s = local_suffix(tag)
if isempty(tag), s = ''; else, s = ['_' tag]; end
end
