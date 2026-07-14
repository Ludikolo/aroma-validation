% RUN_COMPARISON  Part 4: run the four controllers head to head on a hard
% stressed scenario and save the traces + metrics. The figure (viz_compare.m)
% renders from the saved result, so the expensive run happens once.
%
% Controllers (all on the SAME plant, scenario, warmup and rate caps; the two
% plant-based baselines keep the wider default boxes, see the fair_comparison harness):
%   KPC            : direct multi-step Koopman prediction + one QP (this work)
%   Koopman-LMPC   : same Koopman lift but an iterated one-step model + one QP
%   Jacobian-LMPC  : re-linearise the nonlinear plant every step + one QP
%   NMPC           : full nonlinear plant model + fmincon to convergence
%
% The two Koopman controllers carry the demand-adequacy flow floor from the MPC
% formulation; for a clean apples-to-apples we give the SAME floor to the two
% plant-based baselines, so the only thing the comparison tests is the predictor
% and the optimiser, not who got a demand-aware constraint. Baseline solves apply
% their last iterate even on a failed exitflag (10 NMPC steps report ef = -2 in the
% saved run); the fair_comparison harness holds u_prev instead, the fallback KPC uses.

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(root); addpath(fullfile(root, 'controller_closed_loop', 'lmpc')); startup();

% --- scenario: a severe stressed day (harder than the controller part), in the capacity-binding regime ---
% A full 24 h is the right window: worst-consumer demand-met is a DAILY adequacy
% metric (the controller is free to lead/recover around a peak), so it rewards the
% full-horizon planning that distinguishes the controllers. A short peak-only
% window just measures the shared capacity bind at the peak and hides that.
mdot_scale = 0.35;
sc_start   = 0;
T_warm     = 30 * 60;
T_sim      = 24 * 3600;

p = params();
pred_dir = fullfile(root, 'predictor_open_loop', 'results');
ctrl_dir = fullfile(root, 'controller_closed_loop', 'results');

% --- KPC controller (predictor maps + locked tune) ---
F  = load(fullfile(pred_dir, 'vseq_fits_full.mat'));
Hb = load(fullfile(ctrl_dir, 'best_horizons.mat'));
Ab = load(fullfile(ctrl_dir, 'best_alpha.mat'));
Bb = load(fullfile(ctrl_dir, 'best_tune.mat'));
pred = truncate_pred_to_Np(build_kpc_v2_matrices(F.fits), Hb.Np_best);
tune = Bb.best_tune;
tune.dT_0s_max = 7.5; tune.dr_q_max = 4.5; tune.dT_ir_max = 15.0;
tune.use_kirchhoff = true; tune.use_Tr_le_Ts = true; tune.use_mixing = true; tune.rho_slack_mix = 1;
tune.Nc = Hb.Nc_best; tune.alpha_energy = Ab.alpha_star;
Np = Hb.Np_best;

% iterated one-step Koopman model for the Koopman-LMPC baseline
S = load(fullfile(pred_dir, 'iterated_AB.mat'));
if isfield(S, 'fit_iter'), fit_iter = S.fit_iter; else, fit_iter = S; end
tune_lmpc = tune; tune_lmpc.use_mixing = false;

% --- plant at stressed capacity + a shared warmup state ---
[net, z0_cold] = build_plant(p);
net.mdotEdges = mdot_scale * net.mdotEdges; net.flows = net.mdotEdges; net.q0 = net.mdotEdges(:);
ei = edge_user_index(net); n_user = ei.n_user;
F0_idx = find(strcmp({net.Nodes.name}, 'F0'));
R0_idx = find(strcmp({net.Nodes.name}, 'R0'));
con_idx = arrayfun(@(i) find(strcmp({net.Nodes.name}, ei.consumers{i})), 1:n_user);

p_wu = p; p_wu.t_offset = sc_start; p_wu.r_q_fun = @(t) net.mdotEdges(:);
res_wu = simulate_plant(net, z0_cold, p_wu, @(t) 0, @(t) 1.0, T_warm);

% --- shared config for the plant-based baselines ---
K = build_incidence_v2(net);
cfg.Np = Np; cfg.net = net; cfg.p = p; cfg.Ts = p.Ts; cfg.n_user = n_user;
cfg.md = net.mdotEdges(:);
cfg.N  = null([K.M_supply; K.M_return]);
cfg.T0s_lo = p.Tin_nom + p.u_min; cfg.T0s_hi = p.Tin_nom + p.u_max;
cfg.dT0s = 7.5; cfg.drq = 4.5;
% demand-adequacy floor parameters (same as the Koopman MPCs)
cfg.user     = ei.user;
cfg.rqhi_usr = p.excite.r_q_hi_factor * cfg.md(ei.user);
cfg.adeq     = tune.adequacy_safety;
cfg.dT_floor = p.Tin_nom - p.consumer.T_r_min;
% NMPC runs to convergence: the budget is large enough that fmincon stops on its
% optimality tolerance, not on the iteration cap, so the solve-time and feasibility
% it reports are genuine (no artificial cap distorting the comparison).
cfg.fopts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
                         'MaxFunctionEvaluations', 3000, 'MaxIterations', 300, ...
                         'FiniteDifferenceStepSize', 0.05, ...   % large enough to see the delayed delivery response, like the LMPC FD step (a tiny default step gives noisy gradients on this plant)
                         'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-8);
% same QP tolerance as KPC so the two QP-based controllers are solved alike
cfg.qopts = optimoptions('quadprog', 'Display', 'off', ...
                         'OptimalityTolerance', 1e-9, 'ConstraintTolerance', 1e-9, ...
                         'StepTolerance', 1e-12, 'MaxIterations', 200);

% --- run all four on the same window ---
fprintf('Running KPC ...\n'); t0 = tic;
res_kpc = kpc_step_loop(net, res_wu, p, sc_start, T_warm, T_sim, pred, tune, ei, F0_idx, R0_idx, con_idx);
fprintf('  KPC done (%.1f s wall)\n', toc(t0));

fprintf('Running Koopman-LMPC ...\n'); t0 = tic;
res_klmpc = lmpc_step_loop(net, res_wu, p, sc_start, T_warm, T_sim, fit_iter, Np, tune_lmpc, ei, F0_idx, R0_idx, con_idx);
fprintf('  Koopman-LMPC done (%.1f s wall)\n', toc(t0));

fprintf('Running Jacobian-LMPC ...\n'); t0 = tic;
res_jlmpc = baseline_cl_run(net, res_wu, p, sc_start, T_warm, T_sim, ei, F0_idx, R0_idx, con_idx, Np, ...
                            @(z, q, d, up, tb) jlmpc_solve(z, q, d, up, tb, cfg));
fprintf('  Jacobian-LMPC done (%.1f s wall)\n', toc(t0));

fprintf('Running NMPC ...\n'); t0 = tic;
res_nmpc = baseline_cl_run(net, res_wu, p, sc_start, T_warm, T_sim, ei, F0_idx, R0_idx, con_idx, Np, ...
                           @(z, q, d, up, tb) nmpc_solve(z, q, d, up, tb, cfg));
fprintf('  NMPC done (%.1f s wall)\n', toc(t0));

% --- metrics ---
names = {'KPC', 'Koopman-LMPC', 'Jacobian-LMPC', 'NMPC'};
R = {res_kpc, res_klmpc, res_jlmpc, res_nmpc};
fprintf('\nComparison (stressed mdot=%.2f, %.1f h window, Np=%d):\n', mdot_scale, T_sim/3600, Np);
fprintf('%-15s  worst-met%%  unmet[kWh]  med-solve  max-solve  feas\n', 'controller');
M = struct();
for j = 1:numel(R)
    r = R{j};
    pc  = 100 * sum(r.c_i, 2) ./ max(sum(r.d_i, 2), 1e-9);
    unm = sum(max(r.d_i - r.c_i, 0), 'all') * p.Ts / 3.6e6;
    Esr = p.cp * sum(r.Q_net .* (r.T_F0 - r.T_0r)) * p.Ts / 3.6e6;
    feas = sum(r.exitflag > 0);
    M.(matlab.lang.makeValidName(names{j})) = struct('percons', pc, 'worst', min(pc), ...
        'unmet', unm, 'Esrc', Esr, 'med', median(r.solve_ms), 'mx', max(r.solve_ms), 'feas', feas, 'n', numel(r.exitflag));
    fprintf('%-15s  %8.3f  %9.2f  %8.0f  %8.0f  %d/%d\n', ...
        names{j}, min(pc), unm, median(r.solve_ms), max(r.solve_ms), feas, numel(r.exitflag));
end

save(fullfile(here, 'comparison.mat'), 'res_kpc', 'res_klmpc', 'res_jlmpc', 'res_nmpc', ...
     'names', 'M', 'mdot_scale', 'sc_start', 'T_sim', 'Np', 'p', 'ei', '-v7.3');
fprintf('\nSaved comparison.mat\n');


% --- local baseline optimisers (Jacobian-LMPC + NMPC) ---
function C = roll(x, z_state, q0, t_base, cfg)
T0s = x(1); rq = cfg.md + cfg.N * x(2:end);
C = zeros(cfg.n_user, cfg.Np); z = z_state; qc = q0;
for h = 1:cfg.Np
    ps = cfg.p; ps.r_q_fun = @(t) rq; ps.t_offset = t_base + (h-1)*cfg.Ts;
    nh = cfg.net; nh.q0 = qc;
    r = simulate_plant(nh, z, ps, @(t) T0s - cfg.p.Tin_nom, @(t) 1.0, cfg.Ts);
    C(:, h) = r.c_i(:, end); z = r.z_final; qc = r.q_edges(:, end);
end
end

function [A, bcon, lb, ub] = build_cons(u_prev, cfg, d_now)
nxi = size(cfg.N, 2);
P   = [zeros(numel(cfg.md), 1), cfg.N];     % r_q = md + P*x
rqp = u_prev(2:end); md = cfg.md; e1 = [1, zeros(1, nxi)];
% rows: flow box (P*x <= 0.5 md -> r_q <= 1.5 md; -P*x <= 0.7 md ->
% r_q >= 0.3 md, the PRBS lower factor), flow rate caps, T0s rate caps
A = [ P; -P; P; -P; e1; -e1 ];
bcon = [ 0.5*md; 0.7*md; cfg.drq + (rqp - md); cfg.drq - (rqp - md); ...
         u_prev(1) + cfg.dT0s; -(u_prev(1) - cfg.dT0s) ];
% demand-adequacy floor on the user-stub flows (same as the Koopman MPCs)
floor_val = min(cfg.adeq * d_now ./ (cfg.p.cp * cfg.dT_floor), 0.95 * cfg.rqhi_usr);
Pu = P(cfg.user, :);
A    = [A; -Pu];                            % md_u + Pu*x >= floor
bcon = [bcon; md(cfg.user) - floor_val];
lb = [cfg.T0s_lo; -inf(nxi, 1)];
ub = [cfg.T0s_hi;  inf(nxi, 1)];
end

function [u_opt, sms, ef] = nmpc_solve(z_state, q0, d_seq, u_prev, t_base, cfg)
x0 = [u_prev(1); cfg.N' * (u_prev(2:end) - cfg.md)];
[A, bcon, lb, ub] = build_cons(u_prev, cfg, d_seq(:, 1));
obj = @(x) sum((roll(x, z_state, q0, t_base, cfg)/1e3 - d_seq/1e3).^2, 'all');
t0 = tic;
[xo, ~, ef] = fmincon(obj, x0, A, bcon, [], [], lb, ub, [], cfg.fopts);
sms = toc(t0) * 1000;
u_opt = [xo(1); cfg.md + cfg.N * xo(2:end)];
end

function [u_opt, sms, ef] = jlmpc_solve(z_state, q0, d_seq, u_prev, t_base, cfg)
t0 = tic;
x0 = [u_prev(1); cfg.N' * (u_prev(2:end) - cfg.md)]; nx = numel(x0);
C0 = roll(x0, z_state, q0, t_base, cfg)/1e3; c0 = C0(:); dvec = d_seq(:)/1e3;
Jac = zeros(numel(c0), nx);
% FD step: 2 C on the supply temp, large enough to see the delayed delivery
% response through the transport lag (a smaller step gives a noisy gradient and
% the controller stalls at the temperature floor); 0.05 kg/s on each flow DOF
step = [2.0; 0.05 * ones(nx-1, 1)];
for j = 1:nx
    dx = zeros(nx, 1); dx(j) = step(j);
    Cp = roll(x0 + dx, z_state, q0, t_base, cfg)/1e3;
    Jac(:, j) = (Cp(:) - c0) / step(j);
end
H = 2*(Jac'*Jac) + 1e-6*eye(nx); f = 2*Jac'*(c0 - Jac*x0 - dvec);
[A, bcon, lb, ub] = build_cons(u_prev, cfg, d_seq(:, 1));
[xo, ~, ef] = quadprog(H, f, A, bcon, [], [], lb, ub, x0, cfg.qopts);
sms = toc(t0) * 1000;
u_opt = [xo(1); cfg.md + cfg.N * xo(2:end)];
end
