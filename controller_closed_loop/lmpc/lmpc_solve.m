function [u_opt, info] = lmpc_solve(z_0, d_seq, u_prev, pred, p, net, tune)
% LMPC_SOLVE  One MPC step of the iterated one-step Koopman LMPC. This is
% the same QP as the production KPC solver: same decision vector
% [U; S_hi; S_lo; S_Ts], same soft c <= d / c >= 0 / T_r <= T_s
% constraints, same rate caps, same demand-adequacy floor, same Kirchhoff
% equality, same energy term and move-blocking. The whole point is a
% single-variable swap, so only two things differ from the KPC solver:
%
%   1. The prediction is the iterated roll-out, linearised around the
%      previous step, so c_lin = cp (F z_0 + G U + b) carries the offset
%      pred.b and the T_r <= T_s constraint uses the predicted T_s
%      F_Ts z_0 + Ts_off (the demand-driven part of T_s lives in the
%      offset, exactly as KPC's demand-in-V extras put it there). KPC has
%      no b because its V_seq maps predict theta directly.
%
%   2. There is no mixing-residual block. The iterated A,B does not fit
%      the per-junction mixing residual, so the controller cannot enforce
%      it. Everything else is shared.
%
% The tune is the KPC tune with use_mixing off; every other knob behaves
% identically.

Np      = pred.Np;
n_user  = pred.n_user;
n_u     = pred.n_u;
n_U     = Np * n_u;
n_S     = Np * n_user;
n_edges = numel(net.Edges);

% predictor shape invariants
assert(size(pred.F, 1) == n_user * Np, 'lmpc_solve: pred.F row count mismatch');
assert(size(pred.G, 1) == n_user * Np, 'lmpc_solve: pred.G row count mismatch');
assert(size(pred.G, 2) == n_U,         'lmpc_solve: pred.G column count mismatch');
assert(numel(z_0) == size(pred.F, 2),  'lmpc_solve: z_0 length mismatch with pred.F');
assert(isequal(size(d_seq), [n_user, Np]), 'lmpc_solve: d_seq must be n_user x Np');

% x = [U; S_hi; S_lo; S_Ts]   (no mixing slack)
n_total = n_U + 3 * n_S;

%% per-step bounds on u = [T_0s; r_q (per edge); T_ir (per user)]
T_0s_lo = p.Tin_nom + p.u_min;
T_0s_hi = p.Tin_nom + p.u_max;
if isfield(tune, 'T_0s_min_floor') && ~isempty(tune.T_0s_min_floor)
    T_0s_lo = max(T_0s_lo, tune.T_0s_min_floor);
end
lo_factor = p.excite.r_q_lo_factor;
if isfield(tune, 'r_q_lo_factor') && ~isempty(tune.r_q_lo_factor)
    lo_factor = tune.r_q_lo_factor;
end
r_q_lo = lo_factor              * net.mdotEdges(:);
r_q_hi = p.excite.r_q_hi_factor * net.mdotEdges(:);
T_ir_lo = p.consumer.T_r_min    * ones(n_user, 1);
T_ir_hi = (p.Tin_nom + p.u_max) * ones(n_user, 1);

u_lb = [T_0s_lo; r_q_lo; T_ir_lo];
u_ub = [T_0s_hi; r_q_hi; T_ir_hi];
LB_U = repmat(u_lb, Np, 1);
UB_U = repmat(u_ub, Np, 1);

LB_S = zeros(n_S, 1);
UB_S = inf(n_S, 1);
LB = [LB_U; LB_S; LB_S; LB_S];
UB = [UB_U; UB_S; UB_S; UB_S];

%% c_lin <= d (soft, S_hi) and c_lin >= 0 (soft, S_lo)
% c_lin = cp (F z_0 + G U + b), so the constant part is cp (F z_0 + b).
% The b offset is what carries the op-point term and the known demand
% roll-out; KPC carries the demand the same way through G_d D.
D_vec = reshape(d_seq', [], 1);   % consumer-major, matches the predictor
F_z0  = pred.F * z_0 + pred.b;     % theta_lin without cp

A_cd_pos = sparse([ p.cp * pred.G,  -speye(n_S),       sparse(n_S, n_S), sparse(n_S, n_S)]);
A_cd_neg = sparse([-p.cp * pred.G,   sparse(n_S, n_S), -speye(n_S),      sparse(n_S, n_S)]);
A_ineq = [A_cd_pos; A_cd_neg];
b_ineq = [ D_vec - p.cp * F_z0;
                   p.cp * F_z0 ];

%% rate constraints |u_k - u_{k-1}| <= du_max
[D_U, D_pre] = build_diff_op(Np, n_u);
du_max_step = [tune.dT_0s_max;
               tune.dr_q_max  * ones(n_edges, 1);
               tune.dT_ir_max * ones(n_user, 1)];
DU_MAX = repmat(du_max_step, Np, 1);

A_rate_pos = sparse([ D_U,  sparse(n_U, 3 * n_S)]);
b_rate_pos = DU_MAX - D_pre * u_prev;
A_rate_neg = sparse([-D_U,  sparse(n_U, 3 * n_S)]);
b_rate_neg = DU_MAX + D_pre * u_prev;
A_ineq = [A_ineq; A_rate_pos; A_rate_neg];
b_ineq = [b_ineq; b_rate_pos; b_rate_neg];

%% T_r^i <= T_s^i_pred (soft, slack S_Ts)
% T_ir lives at U index (h-1)*n_u + 1 + n_edges + i. The predicted T_s is
% F_Ts z_0 + Ts_off, where Ts_off is the demand-driven part of the roll-out.
use_Tr_Ts = isfield(tune, 'use_Tr_le_Ts') && tune.use_Tr_le_Ts ...
            && pred.has_extras;
if use_Tr_Ts
    n_Tr = n_user * Np;
    rows_E = zeros(n_Tr, 1);
    cols_E = zeros(n_Tr, 1);
    for h = 1:Np
        for i = 1:n_user
            row = (i-1)*Np + h;
            col = (h-1)*n_u + 1 + n_edges + i;
            rows_E(row) = row;
            cols_E(row) = col;
        end
    end
    E_Tir = sparse(rows_E, cols_E, ones(n_Tr,1), n_Tr, n_U);
    A_Tr  = sparse([E_Tir - pred.G_Ts, ...
                    sparse(n_Tr, n_S), sparse(n_Tr, n_S), -speye(n_Tr)]);
    b_Tr  = pred.F_Ts * z_0 + pred.Ts_off;
    A_ineq = [A_ineq; A_Tr];
    b_ineq = [b_ineq; b_Tr];
end

%% cost
% delivery reward: -cp 1^T (F z_0 + G U + b). The constant cp (F z_0 + b)
% drops out of argmin, so the demand offset in b does not change the
% gradient; only G U enters. Same form as KPC.
f_track_U  = -p.cp * (sum(pred.G, 1)');

f_slack_hi = tune.rho_slack * ones(n_S, 1);
rho_neg_ratio = 10;
if isfield(tune, 'rho_neg_ratio') && ~isempty(tune.rho_neg_ratio)
    rho_neg_ratio = tune.rho_neg_ratio;
end
assert(rho_neg_ratio > 0, 'lmpc_solve: rho_neg_ratio must be strictly positive');
rho_neg = tune.rho_slack * rho_neg_ratio;
if isfield(tune, 'rho_slack_neg') && ~isempty(tune.rho_slack_neg)
    rho_neg = tune.rho_slack_neg;
end
f_slack_lo = rho_neg * ones(n_S, 1);
rho_Ts = rho_neg;
if isfield(tune, 'rho_slack_Ts') && ~isempty(tune.rho_slack_Ts)
    rho_Ts = tune.rho_slack_Ts;
end
f_slack_Ts = rho_Ts * ones(n_S, 1);

% plant-side energy term alpha (T_0s - T_0r) Q_net, linearised the same
% way as KPC around the previous-step (T_0s, T_0r, Q_net). q and T_0r feed
% it linearly in U; the constant demand offsets in their roll-out cancel
% in the gradient, so we use F_q/G_q and F_T0r/G_T0r without offsets.
use_energy = isfield(tune, 'alpha_energy') && tune.alpha_energy > 0 ...
             && pred.has_extras;
if use_energy
    alpha_E  = tune.alpha_energy;
    T0s_nom  = pred.opp.T_0s_0;
    T0r_nom  = pred.opp.T_0r_0;
    Qnet_nom = pred.opp.Q_net_0;

    e_T0s = zeros(n_U, 1);
    for h = 1:Np
        e_T0s((h-1)*n_u + 1) = 1;
    end
    sum_G_q   = sum(pred.G_q,   1)';
    sum_G_T0r = sum(pred.G_T0r, 1)';

    f_energy = alpha_E * (T0s_nom - T0r_nom) * sum_G_q ...
             + alpha_E *  Qnet_nom          * (e_T0s - sum_G_T0r);
else
    f_energy = zeros(n_U, 1);
end

% quadratic move-suppression on U
W      = diag([tune.R_du_T0s;
               tune.R_du_rq  * ones(n_edges, 1);
               tune.R_du_Tir * ones(n_user, 1)]);
W_full = kron(eye(Np), W);
D_U_s    = sparse(D_U);
D_pre_s  = sparse(D_pre);
W_full_s = sparse(W_full);
H_du_UU = 2 * (D_U_s' * W_full_s * D_U_s);
f_du_U  = 2 * (D_U_s' * W_full_s * D_pre_s * u_prev);

H = blkdiag(H_du_UU, sparse(3 * n_S, 3 * n_S)) + 1e-9 * speye(n_total);
f = [f_track_U + f_du_U + f_energy; f_slack_hi; f_slack_lo; f_slack_Ts];

%% demand-adequacy floor on the user-stub r_q (identical to KPC)
ei = edge_user_index(net);
dT_floor = p.Tin_nom - p.consumer.T_r_min;
if isfield(tune, 'dT_floor') && ~isempty(tune.dT_floor)
    dT_floor = tune.dT_floor;
end
assert(dT_floor > 0, 'lmpc_solve: dT_floor = Tin_nom - T_r_min must be > 0');
adequacy_safety = 1.2;
if isfield(tune, 'adequacy_safety') && ~isempty(tune.adequacy_safety)
    adequacy_safety = tune.adequacy_safety;
end
r_q_user_upper = p.excite.r_q_hi_factor * net.mdotEdges(ei.user);

n_dmin = n_user * Np;
rows_D = zeros(n_dmin, 1);
cols_D = zeros(n_dmin, 1);
b_dmin = zeros(n_dmin, 1);
for h = 1:Np
    for i = 1:n_user
        row = (i-1)*Np + h;
        col = (h-1)*n_u + 1 + ei.user(i);
        rows_D(row) = row;
        cols_D(row) = col;
        floor_req = adequacy_safety * d_seq(i, h) / (p.cp * dT_floor);
        floor_cap = 0.95 * r_q_user_upper(i);
        b_dmin(row) = -min(floor_req, floor_cap);
    end
end
E_rq_user = sparse(rows_D, cols_D, ones(n_dmin, 1), n_dmin, n_U);
A_dmin = sparse([-E_rq_user, sparse(n_dmin, 3 * n_S)]);
A_ineq = [A_ineq; A_dmin];
b_ineq = [b_ineq; b_dmin];

%% Kirchhoff equality on r_q (+ optional move-blocking on tune.Nc)
use_kirch = isfield(tune, 'use_kirchhoff') && tune.use_kirchhoff;
if use_kirch
    K = build_incidence_v2(net);
    n_eq_per_h = K.n_sn + K.n_rn;
    M_block = [K.M_supply; K.M_return];
    [Mi, Mj, Mv] = find(M_block);

    nnz_per_h = numel(Mv);
    II = zeros(Np * nnz_per_h, 1);
    JJ = zeros(Np * nnz_per_h, 1);
    VV = zeros(Np * nnz_per_h, 1);
    for h = 1:Np
        rq_col_offset = (h-1)*n_u + 1;
        row_offset    = (h-1)*n_eq_per_h;
        slice = (h-1)*nnz_per_h + (1:nnz_per_h);
        II(slice) = Mi + row_offset;
        JJ(slice) = Mj + rq_col_offset;
        VV(slice) = Mv;
    end
    Aeq = sparse(II, JJ, VV, n_eq_per_h * Np, n_total);
    assert(nnz(Aeq) == Np * nnz_per_h, 'lmpc_solve: Kirchhoff nnz mismatch');
    beq = zeros(n_eq_per_h * Np, 1);
else
    Aeq = [];
    beq = [];
end

Nc = Np;
if isfield(tune, 'Nc') && ~isempty(tune.Nc)
    Nc = max(1, min(Np, round(tune.Nc)));
end
if Nc < Np
    n_block_steps = Np - Nc;
    n_block_eq    = n_block_steps * n_u;
    II_b = (1:n_block_eq).';
    II_b = [II_b; II_b];
    JJ_b = zeros(2 * n_block_eq, 1);
    VV_b = zeros(2 * n_block_eq, 1);
    row = 0;
    for h_idx = Nc:Np-1
        for j = 1:n_u
            row = row + 1;
            JJ_b(row)             = h_idx*n_u + j;
            JJ_b(n_block_eq+row)  = (Nc-1)*n_u + j;
            VV_b(row)             = +1;
            VV_b(n_block_eq+row)  = -1;
        end
    end
    A_block = sparse(II_b, JJ_b, VV_b, n_block_eq, n_total);
    b_block = zeros(n_block_eq, 1);
    Aeq = [Aeq; A_block];
    beq = [beq; b_block];
end

%% solve (same solver settings as KPC)
qp_algo     = 'interior-point-convex';
qp_opt_tol  = 1e-9;
qp_con_tol  = 1e-9;
qp_step_tol = 1e-12;
qp_max_iter = 200;
if isfield(tune, 'qp_algorithm') && ~isempty(tune.qp_algorithm)
    qp_algo = tune.qp_algorithm;
end
if isfield(tune, 'qp_opt_tol') && ~isempty(tune.qp_opt_tol)
    qp_opt_tol = tune.qp_opt_tol;
end
if isfield(tune, 'qp_con_tol') && ~isempty(tune.qp_con_tol)
    qp_con_tol = tune.qp_con_tol;
end
if isfield(tune, 'qp_step_tol') && ~isempty(tune.qp_step_tol)
    qp_step_tol = tune.qp_step_tol;
end
if isfield(tune, 'qp_max_iter') && ~isempty(tune.qp_max_iter)
    qp_max_iter = tune.qp_max_iter;
end
assert(qp_opt_tol > 0 && qp_con_tol > 0 && qp_step_tol > 0, ...
    'lmpc_solve: quadprog tolerances must be strictly positive');
qp_opts = optimoptions('quadprog', ...
    'Display', 'off', ...
    'Algorithm', qp_algo, ...
    'OptimalityTolerance', qp_opt_tol, ...
    'ConstraintTolerance', qp_con_tol, ...
    'StepTolerance',       qp_step_tol, ...
    'MaxIterations',       qp_max_iter);

if isfield(tune, 'warm_start_U') && ~isempty(tune.warm_start_U) ...
        && numel(tune.warm_start_U) == n_U
    U_init = tune.warm_start_U(:);
else
    U_init = repmat(u_prev, Np, 1);
end
x0 = [U_init;
      1e6 * ones(n_S, 1);
      1e6 * ones(n_S, 1);
      1e6 * ones(n_S, 1)];
x0 = max(LB, min(UB, x0));

tic;
[x_opt, fval, exitflag, output] = quadprog(H, f, A_ineq, b_ineq, Aeq, beq, LB, UB, x0, qp_opts);
solve_ms = toc * 1000;

if exitflag <= 0
    warning('lmpc_solve: quadprog exitflag %d (%s)', exitflag, output.message);
    x_opt = x0;
end
if any(~isfinite(x_opt(1:n_U)))
    error('lmpc_solve: non-finite values in U after solve');
end

U_opt    = x_opt(1:n_U);
S_hi_opt = x_opt(n_U +     (1:n_S));
S_lo_opt = x_opt(n_U + 1*n_S + (1:n_S));
S_Ts_opt = x_opt(n_U + 2*n_S + (1:n_S));

u_opt = U_opt(1:n_u);
info.U_seq      = reshape(U_opt,    n_u,    Np);
info.S_seq      = reshape(S_hi_opt, Np, n_user).';
info.S_lo_seq   = reshape(S_lo_opt, Np, n_user).';
info.S_Ts_seq   = reshape(S_Ts_opt, Np, n_user).';
info.fval       = fval;
info.exitflag   = exitflag;
info.iterations = output.iterations;
info.solve_ms   = solve_ms;
info.theta_pred = pred.F * z_0 + pred.G * U_opt + pred.b;
info.c_pred     = p.cp * info.theta_pred;
end


function [D_U, D_pre] = build_diff_op(Np, n_u)
% difference operator: rate_k = u_k - u_{k-1} with u_{-1} = u_prev
n_rows = Np * n_u;
D_U   = zeros(n_rows, n_rows);
D_pre = zeros(n_rows, n_u);

D_U(1:n_u, 1:n_u) = eye(n_u);
D_pre(1:n_u, :)   = -eye(n_u);

for k = 1:Np-1
    rk = k*n_u + (1:n_u);
    D_U(rk, (k-1)*n_u + (1:n_u)) = -eye(n_u);
    D_U(rk,  k*n_u    + (1:n_u)) =  eye(n_u);
end
end
