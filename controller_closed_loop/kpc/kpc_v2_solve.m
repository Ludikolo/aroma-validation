function [u_opt, info] = kpc_v2_solve(z_0, d_seq, u_prev, pred, p, net, tune)
% KPC_V2_SOLVE  One MPC step of the v2 Koopman predictive controller.
%
% Solves the spec's eq:kpc_problem_dynamic_flow as a single QP. The
% inequalities c <= d, c >= 0, T_r <= T_s and the mixing law are
% softened with slack variables (the V_seq predictor cannot guarantee
% them at every step). Kirchhoff conservation is imposed on the
% reference r_q rather than on q itself; with q_0 conservative this
% keeps q conservative for all k.

Np      = pred.Np;
n_user  = pred.n_user;
n_u     = pred.n_u;
n_U     = Np * n_u;
n_S     = Np * n_user;
n_edges = numel(net.Edges);

% pred.F is (n_user*Np x n_z), pred.G is (n_user*Np x Np*n_u), and d_seq is
% n_user x Np stacked consumer-major to match the theta read-out below.

% mixing slack only if both tune and predictor support it
use_mix_constr = isfield(tune, 'use_mixing') && tune.use_mixing ...
                 && isfield(pred, 'has_mixing') && pred.has_mixing;
has_mix = use_mix_constr;
if has_mix
    n_mix_total = pred.n_mix * Np;
else
    n_mix_total = 0;
end

% x = [U; S_hi; S_lo; S_Ts; S_mix]
n_total = n_U + 3 * n_S + n_mix_total;

% per-step bounds
% u_k = [T_0s ; r_q (one per edge) ; T_ir (one per user)]
T_0s_lo = p.Tin_nom + p.u_min;
T_0s_hi = p.Tin_nom + p.u_max;

% Optional thermal-capacity floor: keeps T_0s well above T_r_min so a
% low-demand window cannot drain the network of useable T_s before the
% next ramp. Off by default; the smoke tests exercise the raw box.
if isfield(tune, 'T_0s_min_floor') && ~isempty(tune.T_0s_min_floor)
    T_0s_lo = max(T_0s_lo, tune.T_0s_min_floor);
end

% Box lower on r_q: default = the PRBS lower factor used during data
% generation. Tune can raise it to enforce minimum circulation.
kpc_lo_factor = p.excite.r_q_lo_factor;
if isfield(tune, 'r_q_lo_factor') && ~isempty(tune.r_q_lo_factor)
    kpc_lo_factor = tune.r_q_lo_factor;
end
r_q_lo = kpc_lo_factor          * net.mdotEdges(:);
r_q_hi = p.excite.r_q_hi_factor * net.mdotEdges(:);

% T_ir bounds: T_r >= T_r_min hard; T_r <= T_s is enforced softly via
% S_Ts (T_s is read off the lift, not part of u). Upper here is a safe
% box covering any T_s in the source range.
T_ir_lo = p.consumer.T_r_min    * ones(n_user, 1);
T_ir_hi = (p.Tin_nom + p.u_max) * ones(n_user, 1);

u_lb = [T_0s_lo; r_q_lo; T_ir_lo];
u_ub = [T_0s_hi; r_q_hi; T_ir_hi];

LB_U = repmat(u_lb, Np, 1);
UB_U = repmat(u_ub, Np, 1);

% all slacks non-negative, unbounded above
LB_S = zeros(n_S, 1);
UB_S = inf(n_S, 1);

LB = [LB_U; LB_S; LB_S; LB_S; zeros(n_mix_total, 1)];
UB = [UB_U; UB_S; UB_S; UB_S; inf(n_mix_total,  1)];

% inequalities: c <= d (soft, S_hi) and c >= 0 (soft, S_lo)
% c = cp * theta with theta = F z_0 + G U (+ G_d D in d-in-V mode)
% match THETA stacking: consumer-major (all horizons of C1, then C2, ...)
D_vec = reshape(d_seq', [], 1);

% Effective free term: F z_0 plus the demand contribution from G_d D
% when the predictor was fitted with demand inside V_seq. d_seq is
% known at solve time (passed in by the caller), so G_d * D_vec folds
% into the constant offset and the QP structure stays unchanged.
has_d_in_V = isfield(pred, 'has_demand_in_V') && pred.has_demand_in_V;
F_z0 = pred.F * z_0;
if has_d_in_V
    F_z0 = F_z0 + pred.G_d * D_vec;
end

% Optional bias correction: per-consumer additive offset on predicted
% c^(i). If tune.bias_correction is set (n_user x 1), the predicted heat
% becomes c_pred + bias_correction. In c <= d this tightens the upper
% bound by bias_correction per consumer. Stack consumer-major to match
% D_vec and F_z0.
if isfield(tune, 'bias_correction') && ~isempty(tune.bias_correction)
    bias_corr = tune.bias_correction(:);
    assert(numel(bias_corr) == n_user, 'tune.bias_correction must be n_user x 1');
    bias_vec = reshape(repmat(bias_corr, 1, Np)', [], 1);
else
    bias_vec = zeros(n_S, 1);
end

A_cd_pos = sparse([ p.cp * pred.G,  -speye(n_S),       sparse(n_S, n_S), sparse(n_S, n_S), sparse(n_S, n_mix_total)]);
A_cd_neg = sparse([-p.cp * pred.G,   sparse(n_S, n_S), -speye(n_S),      sparse(n_S, n_S), sparse(n_S, n_mix_total)]);
A_ineq = [A_cd_pos; A_cd_neg];
b_ineq = [ D_vec - bias_vec - p.cp * F_z0;
           bias_vec + p.cp * F_z0 ];

% rate constraints, |u_k - u_{k-1}| <= du_max
[D_U, D_pre] = build_diff_op(Np, n_u);

du_max_step = [tune.dT_0s_max;
               tune.dr_q_max  * ones(n_edges, 1);
               tune.dT_ir_max * ones(n_user, 1)];
DU_MAX = repmat(du_max_step, Np, 1);

% rates only act on U; slack columns are zero
A_rate_pos = sparse([ D_U,  sparse(n_U, 3 * n_S + n_mix_total)]);
b_rate_pos = DU_MAX - D_pre * u_prev;
A_rate_neg = sparse([-D_U,  sparse(n_U, 3 * n_S + n_mix_total)]);
b_rate_neg = DU_MAX + D_pre * u_prev;

A_ineq = [A_ineq; A_rate_pos; A_rate_neg];
b_ineq = [b_ineq; b_rate_pos; b_rate_neg];

% T_r^i <= T_s^i (soft, slack S_Ts)
% T_ir lives at U index (h-1)*n_u + 1 + n_edges + i.
% T_s prediction: row (i-1)*Np + h of F_Ts z_0 + G_Ts U.
% Hard form is infeasible at low-demand states where T_s_pred dips
% below T_ir; soften and let the plant substation enforce the physical
% bound (cannot extract heat with T_r > T_s).
use_Tr_Ts = isfield(tune, 'use_Tr_le_Ts') && tune.use_Tr_le_Ts ...
            && isfield(pred, 'has_extras') && pred.has_extras;
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
    % T_ir - T_s_pred - s_Ts <= 0,  with x = [U; S_hi; S_lo; S_Ts; S_mix]
    A_Tr  = sparse([E_Tir - pred.G_Ts, ...
                    sparse(n_Tr, n_S), sparse(n_Tr, n_S), -speye(n_Tr), ...
                    sparse(n_Tr, n_mix_total)]);
    b_Tr  = pred.F_Ts * z_0;
    if has_d_in_V && isfield(pred, 'G_d_Ts')
        b_Tr = b_Tr + pred.G_d_Ts * D_vec;
    end
    A_ineq = [A_ineq; A_Tr];
    b_ineq = [b_ineq; b_Tr];
end

% cost
% tracking on U: -cp * 1^T (F z_0 + G U); constant drops out of argmin.
% In ADMM mode (tune.admm_area_id set) only the rows of pred.G for the
% local consumer enter the tracking gradient. Other consumers' tracking
% costs live in their own area-QPs and are summed at consensus.
admm_active = isfield(tune, 'admm_area_id') && ~isempty(tune.admm_area_id);
if admm_active
    a_id = tune.admm_area_id;
    G_area_rows = (a_id - 1) * Np + (1:Np);
    f_track_U = -p.cp * (sum(pred.G(G_area_rows, :), 1)');
else
    f_track_U = -p.cp * (sum(pred.G, 1)');
end
% Slacks live consumer-major; in ADMM mode keep only the local consumer's
% penalties non-zero so the sum across areas reproduces the central cost.
f_slack_hi = tune.rho_slack * ones(n_S, 1);
if admm_active
    area_mask = zeros(n_S, 1);
    area_mask((a_id - 1) * Np + (1:Np)) = 1;
    f_slack_hi = f_slack_hi .* area_mask;
end

% Default penalty ratios for the three constraint families that are
% softer than c <= d. The 10x factor on rho_neg reflects that c < 0 is
% unphysical (heat must flow into the substation), so any breach should
% be much more expensive than the c <= d (over-delivery) slack. The
% T_r <= T_s (return below supply) and mixing residuals inherit this
% scale by default but each is independently overridable via tune.
rho_neg_ratio_default = 10;
if isfield(tune, 'rho_neg_ratio') && ~isempty(tune.rho_neg_ratio)
    rho_neg_ratio_default = tune.rho_neg_ratio;
end
assert(rho_neg_ratio_default > 0, ...
    'kpc_v2_solve: rho_neg_ratio must be strictly positive');
rho_neg = tune.rho_slack * rho_neg_ratio_default;
if isfield(tune, 'rho_slack_neg') && ~isempty(tune.rho_slack_neg)
    rho_neg = tune.rho_slack_neg;
end
f_slack_lo = rho_neg * ones(n_S, 1);
if admm_active
    f_slack_lo = f_slack_lo .* area_mask;
end

rho_Ts = rho_neg;
if isfield(tune, 'rho_slack_Ts') && ~isempty(tune.rho_slack_Ts)
    rho_Ts = tune.rho_slack_Ts;
end
f_slack_Ts = rho_Ts * ones(n_S, 1);
if admm_active
    f_slack_Ts = f_slack_Ts .* area_mask;
end

rho_mix = rho_Ts;
if isfield(tune, 'rho_slack_mix') && ~isempty(tune.rho_slack_mix)
    rho_mix = tune.rho_slack_mix;
end
f_slack_mix = rho_mix * ones(n_mix_total, 1);
if admm_active && n_mix_total > 0
    % Mixing residual is per-junction not per-consumer; split evenly
    % across K = n_user areas so the sum matches the central penalty.
    f_slack_mix = f_slack_mix / n_user;
end

% Plant-side energy term alpha (T_0s - T_0r) Q_net is bilinear; linearise
% around (T_0s_nom, T_0r_nom, Q_net_nom) read from u_prev and z_0.
use_energy = isfield(tune, 'alpha_energy') && tune.alpha_energy > 0 ...
             && isfield(pred, 'has_extras') && pred.has_extras;
if use_energy
    alpha_E  = tune.alpha_energy;
    T0s_nom  = u_prev(1);
    T0r_nom  = z_0(pred.meta.idx.T_0r);
    Qnet_nom = z_0(pred.meta.idx.Q_net);

    e_T0s = zeros(n_U, 1);
    for h = 1:Np
        e_T0s((h-1)*n_u + 1) = 1;
    end

    sum_G_q   = sum(pred.G_q,   1)';        % n_U x 1
    sum_G_T0r = sum(pred.G_T0r, 1)';        % n_U x 1

    f_energy = alpha_E * (T0s_nom - T0r_nom) * sum_G_q ...
             + alpha_E *  Qnet_nom          * (e_T0s - sum_G_T0r);
else
    f_energy = zeros(n_U, 1);
end
% Plant-side energy is shared across consumers; split it 1/K across
% areas in ADMM mode so the K-agent sum matches the central gradient.
if admm_active
    f_energy = f_energy / n_user;
end

% Quadratic on U: rate penalty R_du * (D_U U + D_pre u_prev)^T W (D_U U + D_pre u_prev)
W      = diag([tune.R_du_T0s;
               tune.R_du_rq  * ones(n_edges, 1);
               tune.R_du_Tir * ones(n_user, 1)]);
W_full = kron(eye(Np), W);

D_U_s   = sparse(D_U);
D_pre_s = sparse(D_pre);
W_full_s = sparse(W_full);
H_du_UU = 2 * (D_U_s' * W_full_s * D_U_s);             % n_U x n_U sparse
f_du_U  = 2 * (D_U_s' * W_full_s * D_pre_s * u_prev);
% Rate-cost is a shared U-only term (acts on every agent identically).
% Split it 1/K across areas in ADMM mode so the sum across agents
% reproduces the central rate penalty.
if admm_active
    H_du_UU = H_du_UU / n_user;
    f_du_U  = f_du_U  / n_user;
end

% slacks enter the cost linearly only; tiny diagonal regulariser on H
H = blkdiag(H_du_UU, sparse(3 * n_S + n_mix_total, 3 * n_S + n_mix_total)) + 1e-9 * speye(n_total);
f = [f_track_U + f_du_U + f_energy; f_slack_hi; f_slack_lo; f_slack_Ts; f_slack_mix];

% ADMM consensus penalty on U only (opt-in, no-op when absent)
% Adds (rho/2) * ||U - z_U + y_U/rho||^2 to the cost, with z_U and y_U
% length n_U (so only the input block enters the consensus; per-area
% slacks stay private). With admm_rho missing or zero this block is a
% no-op and kpc_v2_solve is bit-exact unchanged.
if isfield(tune, 'admm_rho') && ~isempty(tune.admm_rho) && tune.admm_rho > 0
    rho_a = tune.admm_rho;
    z_a   = tune.admm_z(:);
    y_a   = tune.admm_y(:);
    assert(numel(z_a) == n_U && numel(y_a) == n_U, ...
        'kpc_v2_solve: admm_z/y must have length n_U = %d (U-only consensus)', n_U);
    % Add rho * I to the U-block of H (top-left n_U x n_U).
    H_admm = sparse(1:n_U, 1:n_U, rho_a * ones(n_U, 1), n_total, n_total);
    H = H + H_admm;
    % Linear shift on the U block of f.
    f(1:n_U) = f(1:n_U) + (y_a - rho_a * z_a);
end

% demand-adequacy floor on the user-stub r_q
% Capacity safeguard: r_q^(i)_h >= safety * d^(i)_h / (cp * dT_floor),
% capped at 0.95 * r_q_max so peak demand stays feasible against the
% rate constraints. Defaults match the published scenario-suite values;
% tune can tighten them for longer-horizon runs.
ei = edge_user_index(net);
dT_floor = p.Tin_nom - p.consumer.T_r_min;
if isfield(tune, 'dT_floor') && ~isempty(tune.dT_floor)
    dT_floor = tune.dT_floor;
end
assert(dT_floor > 0, 'kpc_v2_solve: dT_floor = Tin_nom - T_r_min must be > 0');
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
A_dmin = sparse([-E_rq_user, sparse(n_dmin, 3 * n_S + n_mix_total)]);
A_ineq = [A_ineq; A_dmin];
b_ineq = [b_ineq; b_dmin];

% mixing law as soft equality (spec tilde g_N(z) = 0)
% V_seq predicts the per-junction residual; impose |residual| <= S_mix
% with high penalty on S_mix. Hard equality conflicts with the plant's
% first-order flow non-conservation during transients.
if has_mix
    F_mix_z0 = pred.F_mix * z_0;            % n_mix*Np x 1
    if has_d_in_V && isfield(pred, 'G_d_mix')
        F_mix_z0 = F_mix_z0 + pred.G_d_mix * D_vec;
    end
    A_mix_pos = sparse([ pred.G_mix, sparse(n_mix_total, 3 * n_S), -speye(n_mix_total)]);
    A_mix_neg = sparse([-pred.G_mix, sparse(n_mix_total, 3 * n_S), -speye(n_mix_total)]);
    b_mix_pos = -F_mix_z0;
    b_mix_neg =  F_mix_z0;
    A_ineq = [A_ineq; A_mix_pos; A_mix_neg];
    b_ineq = [b_ineq; b_mix_pos; b_mix_neg];
end

% Kirchhoff equality on r_q
% Imposed on the input rather than on q. With q_0 conservative and
% r_q conservative every step, q stays conservative for all k.
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
    % Sanity: every preallocated triplet should land at a non-zero entry.
    % A mismatch points at a stale nnz_per_h or duplicated (II,JJ) pairs.
    assert(nnz(Aeq) == Np * nnz_per_h, 'kpc_v2_solve: Kirchhoff nnz mismatch');
    beq = zeros(n_eq_per_h * Np, 1);
else
    Aeq = [];
    beq = [];
end

% move blocking: tune.Nc < Np freezes u_h = u_{Nc-1} for h >= Nc
% Standard MPC technique (Cagienard & Morari 2007; Maciejowski 2002
% §3.4). Reduces decision freedom: Nc free decisions u_0..u_{Nc-1},
% then u_{Nc}..u_{Np-1} are clamped to u_{Nc-1}. Default Nc = Np
% (no move blocking).
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
    % U is stacked [u_1; u_2; ...; u_Np]; u_k occupies cols (k-1)*n_u+1:k*n_u.
    % h_idx ranges Nc..Np-1 (so columns h_idx*n_u+j point to u_{Nc+1}..u_Np)
    % and the second column points to u_{Nc} (free slot). This gives Nc free
    % moves and clamps the remaining Np-Nc to the last free input.
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

% solve
% Defaults are the production-locked values. Each is overridable via
% the tune struct so sensitivity-tests (Bucket A4 / A6) can sweep them
% without touching the production path.
qp_algo  = 'interior-point-convex';
qp_opt_tol = 1e-9;
qp_con_tol = 1e-9;
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
% ADMM sub-QPs may need more interior-point iterations because the
% consensus penalty changes the conditioning. Production / non-ADMM
% paths keep the original 200 limit.
if isfield(tune, 'qp_max_iter') && ~isempty(tune.qp_max_iter)
    qp_max_iter = tune.qp_max_iter;
end
assert(qp_opt_tol > 0 && qp_con_tol > 0 && qp_step_tol > 0, ...
    'kpc_v2_solve: quadprog tolerances must be strictly positive');
qp_opts = optimoptions('quadprog', ...
    'Display', 'off', ...
    'Algorithm', qp_algo, ...
    'OptimalityTolerance', qp_opt_tol, ...
    'ConstraintTolerance', qp_con_tol, ...
    'StepTolerance',       qp_step_tol, ...
    'MaxIterations',       qp_max_iter);

% warm-start: previous-step U_seq if supplied, else repeat u_prev.
% kpc_step_loop passes tune.warm_start_U = shifted prior optimum
% to amortise QP iterations across MPC steps.
if isfield(tune, 'warm_start_U') && ~isempty(tune.warm_start_U) ...
        && numel(tune.warm_start_U) == n_U
    U_init = tune.warm_start_U(:);
else
    U_init = repmat(u_prev, Np, 1);
end
x0 = [U_init;
      1e6 * ones(n_S, 1);
      1e6 * ones(n_S, 1);
      1e6 * ones(n_S, 1);
      1e6 * ones(n_mix_total, 1)];
x0 = max(LB, min(UB, x0));

tic;
[x_opt, fval, exitflag, output] = quadprog(H, f, A_ineq, b_ineq, Aeq, beq, LB, UB, x0, qp_opts);
solve_ms = toc * 1000;

if exitflag <= 0
    % Solver failure: hold the warm-start (repeated u_prev). The 50-test
    % regression has never seen this branch fire; if it does in deployment
    % the warning surfaces it and the next step starts from a fresh x0.
    warning('kpc_v2_solve: quadprog exitflag %d (%s)', exitflag, output.message);
    x_opt = x0;
end
% Guard against any pathological NaN/Inf leaking into the applied input.
if any(~isfinite(x_opt(1:n_U)))
    error('kpc_v2_solve: non-finite values in U after solve');
end

U_opt    = x_opt(1:n_U);
S_hi_opt = x_opt(n_U +     (1:n_S));
S_lo_opt = x_opt(n_U + 1*n_S + (1:n_S));
S_Ts_opt = x_opt(n_U + 2*n_S + (1:n_S));
if has_mix
    S_mix_opt = x_opt(n_U + 3*n_S + (1:n_mix_total));
else
    S_mix_opt = [];
end

u_opt = U_opt(1:n_u);
info.U_seq      = reshape(U_opt,    n_u,    Np);
% Consumer-major stacking: S_*_opt((i-1)*Np + h) is consumer i, horizon h.
% reshape(., Np, n_user).' makes info.*_seq(i, h) = consumer i at horizon h.
info.S_seq      = reshape(S_hi_opt, Np, n_user).';   % c <= d slack
info.S_lo_seq   = reshape(S_lo_opt, Np, n_user).';   % c >= 0 slack
info.S_Ts_seq   = reshape(S_Ts_opt, Np, n_user).';   % T_r <= T_s slack
if has_mix
    info.S_mix_seq = reshape(S_mix_opt, Np, pred.n_mix).';
end
info.fval       = fval;
info.exitflag   = exitflag;
info.iterations = output.iterations;
info.solve_ms   = solve_ms;
info.theta_pred = pred.F * z_0 + pred.G * U_opt;
if has_d_in_V
    info.theta_pred = info.theta_pred + pred.G_d * D_vec;
end
info.c_pred     = p.cp * info.theta_pred;
% Full solution + size, used by the ADMM coordinator to keep its
% consensus z and dual y vectors in sync with the QP layout.
info.x_full     = x_opt;
info.n_total    = n_total;

end


function [D_U, D_pre] = build_diff_op(Np, n_u)
% Difference operator: rate_k = u_k - u_{k-1} with u_{-1} = u_prev.
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
