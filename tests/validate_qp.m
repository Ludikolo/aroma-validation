% VALIDATE_QP  Independent verification of the KPC quadratic program on
%   one representative closed-loop step, rebuilt from the shipped demo
%   artefact. Nothing is re-simulated: the mid-day (12:00) state is
%   lifted from the logged demo trajectory exactly as kpc_step_loop
%   lifts it (the delay lags reach 2 samples back, so the truncated
%   history reproduces the demo's z_0), kpc_v2_solve builds and solves
%   the same QP it solves in closed loop, and the returned solution is
%   then checked against the exported QP data (info.H, info.f, ...)
%   rather than against anything the solver reports about itself.
%
%   Q1  Feasibility: the returned x satisfies every inequality,
%       equality and bound of the QP (max residuals reported).
%   Q2  KKT: the multipliers close the stationarity condition
%       H x + f + A' lam_in + Aeq' lam_eq + lam_up - lam_lo = 0,
%       have the right sign, and complementary slackness holds on the
%       inequalities and on the finite bounds.
%   Q3  Cross-solver: re-solving the identical (H, f, A, b, Aeq, beq,
%       LB, UB) with the active-set algorithm reproduces the objective
%       to 1e-6 relative and the same solution on the strictly convex
%       block (the inputs U, where the rate penalty makes H positive
%       definite; the slacks only enter the cost linearly).

clear; clc;
startup;
p = params();

root = fileparts(fileparts(mfilename('fullpath')));
res_dir  = fullfile(root, 'controller_closed_loop', 'results');
pred_res = fullfile(root, 'predictor_open_loop', 'results');

S = load(fullfile(res_dir, 'demo.mat'));
F = load(fullfile(pred_res, 'vseq_fits_full.mat'));
Hz = load(fullfile(res_dir, 'best_horizons.mat'));
pred = truncate_pred_to_Np(build_kpc_v2_matrices(F.fits), Hz.Np_best);
tune = S.tune;

fprintf('\n=== QP validation: one mid-day step, solved and cross-checked ===\n');

%% Rebuild the demo's plant and the mid-day solve inputs
% plant at stressed capacity, exactly as demo_controller builds it
net = build_plant(p);
net.mdotEdges = S.mdot_scale * net.mdotEdges;
net.flows     = net.mdotEdges;
net.q0        = net.mdotEdges(:);
ei = edge_user_index(net);
n_user = ei.n_user;
con_idx = zeros(n_user, 1);
for i = 1:n_user
    con_idx(i) = find(strcmp({net.Nodes.name}, ei.consumers{i}));
end

% closed-loop sample k sits at wall clock T_warm + k*Ts, so the QP the
% demo solved at noon is step k with T_warm + (k-1)*Ts = 12 h
Ts = p.Ts;
k  = round((12 * 3600 - S.T_warm) / Ts) + 1;
t_now = S.T_warm + (k - 1) * Ts;

% controller-side trajectory the lift saw at that step, from the logged
% closed-loop columns 1..k-1 (field names follow kpc_step_loop's traj)
r = S.res_kpc;
hist = 1:(k - 1);
traj = struct();
traj.t_offset = S.sc_start;
traj.t        = S.T_warm + hist * Ts;
traj.T_0s     = r.T_0s(hist);
traj.T_0r     = r.T_0r(hist);
traj.T_F0     = r.T_F0(hist);
traj.T_is     = r.T_is(:, hist);
traj.T_ir     = r.T_ir(:, hist);
traj.q_users  = r.q_users(:, hist);
traj.theta    = r.theta(:, hist);
traj.d        = r.d_i(:, hist);
[Z_hist, ~, ~] = candidate_library(traj, p);
z_0 = Z_hist(:, end);

% demand forecast over the horizon, same call kpc_step_loop makes
Np = pred.Np;
d_seq = zeros(n_user, Np);
for h = 1:Np
    for i = 1:n_user
        d_seq(i, h) = compute_prosumer_Q(S.sc_start + t_now + h * Ts, ...
                                         net.Nodes(con_idx(i)).params);
    end
end
% alignment check: the plant logged its realised demand at the first
% forecast instant, so a time-base error here would show up immediately
d_gap = max(abs(d_seq(:, 1) - r.d_i(:, k))) / max(1, max(abs(r.d_i(:, k))));
assert(d_gap <= 1e-6, ...
    'reconstruction fail: forecast vs logged demand differ by %.2e relative', d_gap);

% previous applied input. T_0s and r_q are logged verbatim; the T_ir
% slot of the previous optimum is not logged, so the plant-side return
% stands in (the same substitution the closed loop makes at its own
% first step, where u_prev comes from the warmup plant state).
u_prev = [r.T_0s(k - 1); r.r_q(:, k - 1); r.T_ir(:, k - 1)];

[u_opt, info] = kpc_v2_solve(z_0, d_seq, u_prev, pred, p, net, tune);
assert(info.exitflag > 0, 'reference solve failed: exitflag %d', info.exitflag);
x   = info.x_full;
n_u = size(info.U_seq, 1);
n_U = n_u * Np;
fprintf(['step k = %d (t = %.1f h), n = %d vars, %d ineq, %d eq, ' ...
         'fval = %.6g, %.0f ms\n'], ...
        k, (S.sc_start + t_now) / 3600, info.n_total, ...
        numel(info.b_ineq), numel(info.beq), info.fval, info.solve_ms);

%% Q1 feasibility residuals of the returned solution
scale_b = max(1, norm([info.b_ineq; info.beq], inf));
res_in  = max([info.A_ineq * x - info.b_ineq; 0]);
res_eq  = max([abs(info.Aeq * x - info.beq); 0]);
res_lb  = max([info.LB - x; 0]);
res_ub  = max([x - info.UB; 0]);
res_all = max([res_in, res_eq, res_lb, res_ub]);
assert(res_all <= 1e-6 * scale_b, ...
    'Q1 fail: max constraint residual %.3e exceeds 1e-6 * %.3g', res_all, scale_b);
fprintf(['Q1 ok: feasible. max residual ineq %.2e, eq %.2e, bounds %.2e ' ...
         '(tol %.1e)\n'], res_in, res_eq, max(res_lb, res_ub), 1e-6 * scale_b);

%% Q2 KKT conditions from the returned multipliers
lam = info.lambda;
g_stat = info.H * x + info.f + info.A_ineq' * lam.ineqlin ...
       + info.Aeq' * lam.eqlin + lam.upper - lam.lower;
rel_stat = norm(g_stat, inf) / max(1, norm(info.f, inf));
assert(rel_stat <= 1e-6, ...
    'Q2 fail: stationarity residual %.3e relative to ||f||', rel_stat);

% dual feasibility: inequality and bound multipliers must be >= 0
min_lam = min([lam.ineqlin; lam.lower; lam.upper]);
assert(min_lam >= -1e-9, 'Q2 fail: negative multiplier %.3e', min_lam);

% complementary slackness: lam .* slack = 0 on inequalities and finite
% bounds; multipliers on the +inf upper bounds must be identically zero
slack   = info.b_ineq - info.A_ineq * x;
comp_in = max(abs(lam.ineqlin .* slack));
comp_lb = max(abs(lam.lower .* (x - info.LB)));
iu      = isfinite(info.UB);
comp_ub = max([abs(lam.upper(iu) .* (info.UB(iu) - x(iu))); 0]);
assert(all(lam.upper(~iu) == 0), 'Q2 fail: nonzero multiplier on an inf bound');
comp_all = max([comp_in, comp_lb, comp_ub]);
assert(comp_all <= 1e-6 * max(1, abs(info.fval)), ...
    'Q2 fail: complementarity residual %.3e vs fval %.3g', comp_all, info.fval);
fprintf(['Q2 ok: KKT holds. stationarity %.2e rel, complementarity %.2e ' ...
         '(vs |fval| %.3g), min multiplier %.1e\n'], ...
        rel_stat, comp_all, abs(info.fval), min_lam);

%% Q3 cross-solver: active-set on the identical QP data
opts_as = optimoptions('quadprog', ...
    'Display', 'off', 'Algorithm', 'active-set', ...
    'OptimalityTolerance', 1e-9, 'ConstraintTolerance', 1e-9, ...
    'MaxIterations', 5000);
[x_as, fval_as, ef_as] = quadprog(full(info.H), info.f, ...
    full(info.A_ineq), info.b_ineq, full(info.Aeq), info.beq, ...
    info.LB, info.UB, x, opts_as);
assert(ef_as > 0, 'Q3 fail: active-set exitflag %d', ef_as);

rel_obj = abs(fval_as - info.fval) / max(1, abs(info.fval));
assert(rel_obj <= 1e-6, ...
    'Q3 fail: objectives differ by %.3e relative (ip %.8g, as %.8g)', ...
    rel_obj, info.fval, fval_as);

% the two solutions must agree on the input block U, where the move
% suppression makes H positive definite (D_U is invertible, W > 0);
% the slacks are linear in the cost so only U is compared pointwise
dU = norm(x_as(1:n_U) - x(1:n_U), inf);
dx = x_as - x;
h_gap = dx' * info.H * dx;
assert(dU <= 1e-3, ...
    'Q3 fail: input blocks differ by %.3e (K or kg/s)', dU);
assert(h_gap <= 1e-6 * max(1, abs(info.fval)), ...
    'Q3 fail: H-norm gap %.3e between the two solutions', h_gap);
fprintf(['Q3 ok: active-set reproduces the QP. objective gap %.2e rel, ' ...
         'max input gap %.2e, H-norm gap %.2e\n'], rel_obj, dU, h_gap);

fprintf('\n=== QP validation: all 3 checks passed ===\n');
