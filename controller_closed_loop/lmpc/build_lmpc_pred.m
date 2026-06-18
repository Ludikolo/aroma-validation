function pred = build_lmpc_pred(fit_iter, Np, z_0, u_prev, p, n_edges, n_user, d_seq)
% BUILD_LMPC_PRED  Per-step prediction matrices for the iterated one-step
% Koopman LMPC. This is the single-variable swap against KPC: same lift,
% same decision vector, same constraints and cost, only the predictor is
% different. KPC predicts theta directly over the whole horizon with the
% direct V_seq maps; this baseline rolls the one-step Koopman model
%
%     z_{k+1} = A z_k + B_ctrl u_k + B_d d_{k+1}
%
% forward Np steps. Stacked over the horizon that gives
%
%     z_h = A^h z_0 + sum_{k=0..h-1} A^{h-1-k} ( B_ctrl u_k + B_d d_{k+1} )
%
% so Z = F_z z_0 + G_z U + Zdem, where U = [u_0; ...; u_{Np-1}] is the
% control sequence (35 inputs per step, exactly KPC's u = [T_0s; r_q(29);
% T_ir(5)]) and Zdem is the contribution of the known demand forecast.
%
% Why split B. The iterated model was fitted with demand as an extra
% input (input vector [u_k; d_{k+1}], width 40), so it carries the same
% demand information KPC has. Refitting without demand would strip that
% information and make the comparison unfair, so we keep it: B_ctrl =
% B(:,1:35) drives the 35 controls, B_d = B(:,36:40) drives the 5-user
% demand. Demand is known at solve time (the same forecast KPC folds in
% via G_d D), so the demand roll-out Zdem is a constant the QP can fold
% into its offset terms. The decision vector stays 35*Np wide, identical
% to KPC.
%
% Bilinear delivered heat. The plant rule is c = cp (T_s - T_ir) q, which
% is bilinear in (T_s, T_ir, q). We linearise it around the previous-step
% operating point (T_s_0, T_ir_0, q_0):
%
%     c_lin/cp = (T_s_0 - T_ir_0) q + q_0 (T_s - T_ir) - q_0 (T_s_0 - T_ir_0)
%
% with T_s and q read off z_h and T_ir read off u_h. This is the standard
% tangent-linearisation MPC: the controller never sees the bilinear
% product, only its first-order approximation around the last step. That
% is the whole point of this baseline; everything downstream is shared
% with KPC.
%
% Inputs
%   fit_iter : struct with A (n_z x n_z), B (n_z x 40), meta.idx (T_is,
%              T_ir, q, T_0r, Q_net, ...)
%   Np       : prediction horizon (samples)
%   z_0      : current lifted state (n_z x 1)
%   u_prev   : previous applied input (35 x 1) = [T_0s; r_q; T_ir]
%   p        : params (uses p.cp)
%   n_edges  : number of network edges (29)
%   n_user   : number of consumers (5)
%   d_seq    : demand forecast over the horizon, n_user x Np. d_seq(:,h)
%              is the demand at horizon h, i.e. the d_{k+1} that drives
%              the step producing z_h. Known at solve time.
%
% Output pred: same fields lmpc_solve and kpc_v2_solve consume.
%   pred.F, pred.G   : c_lin/cp coefficients per (consumer, horizon)
%   pred.b           : per-(consumer, horizon) constant offset, includes
%                      the demand roll-out contribution to c_lin
%   pred.F_Ts,pred.G_Ts, pred.Ts_off : T_s prediction (linear in z) plus
%                      its demand-driven offset, used by T_r <= T_s
%   pred.F_q, pred.G_q, pred.F_T0r, pred.G_T0r : q and T_0r predictions
%                      (linear in z; demand offsets cancel in the cost
%                      gradient, see lmpc_solve, so no offset is carried)
%   pred.has_extras = true, pred.has_mixing = false (iterated A,B does not
%                      fit the per-junction mixing residual)
%   pred.is_lmpc = true : tells the solver to add the b / Ts_off offsets

A = fit_iter.A;
B = fit_iter.B;
idx = fit_iter.meta.idx;
n_z = size(A, 1);

% Split the fitted B: first 35 columns are the controls (same order as
% u = [T_0s(1); r_q(n_edges); T_ir(n_user)]), last n_user columns drive
% the demand forecast d_{k+1}.
n_u   = 1 + n_edges + n_user;     % 35 control inputs
n_dem = n_user;                   % 5 demand inputs
assert(size(B, 2) == n_u + n_dem, ...
    'build_lmpc_pred: B width %d, expected %d (35 controls + %d demand)', ...
    size(B, 2), n_u + n_dem, n_dem);
assert(isequal(size(d_seq), [n_user, Np]), ...
    'build_lmpc_pred: d_seq must be n_user x Np');
B_ctrl = B(:, 1:n_u);
B_d    = B(:, n_u + (1:n_dem));

%% Control prediction matrices F_z, G_z over Np horizons (B_ctrl only)
% z_h = A^h z_0 + sum_{k=0..h-1} A^{h-1-k} B_ctrl u_k  (+ demand part below)
F_z = zeros(Np * n_z, n_z);
A_pow_h = eye(n_z);                 % A^0
for h = 1:Np
    A_pow_h = A_pow_h * A;          % becomes A^h
    rows = (h-1)*n_z + (1:n_z);
    F_z(rows, :) = A_pow_h;
end
% G_z, column by column. Input u_k contributes to z_h for h > k with
% coefficient A^{h-1-k} B_ctrl. Pre-build A^{j-1} B_ctrl for j = 1..Np.
AjB = cell(Np, 1);
AjB{1} = B_ctrl;                    % A^0 B_ctrl
for j = 2:Np
    AjB{j} = A * AjB{j-1};
end
G_z = zeros(Np * n_z, Np * n_u);
for k = 0:Np-1
    cols = k*n_u + (1:n_u);
    for h = (k+1):Np
        rows = (h-1)*n_z + (1:n_z);
        G_z(rows, cols) = AjB{h - k};   % A^{h-1-k} B_ctrl, j = h-k
    end
end

%% Demand roll-out: known constant contribution to each z_h
% zdem_h = sum_{k=0..h-1} A^{h-1-k} B_d d_{k+1}, with d_{k+1} = d_seq(:,k+1).
% A^{h-1-k} B_d uses the same AjB-style powers; reuse the A^{j-1} ladder.
AjBd = cell(Np, 1);
AjBd{1} = B_d;                      % A^0 B_d
for j = 2:Np
    AjBd{j} = A * AjBd{j-1};
end
Zdem = zeros(Np * n_z, 1);          % stacked [zdem_1; ...; zdem_Np]
for h = 1:Np
    acc = zeros(n_z, 1);
    for k = 0:h-1
        acc = acc + AjBd{h - k} * d_seq(:, k+1);   % A^{h-1-k} B_d d_{k+1}
    end
    Zdem((h-1)*n_z + (1:n_z)) = acc;
end

%% Operating-point values (read off the previous step)
T_s_0   = z_0(idx.T_is);                       % n_user x 1
q_0     = z_0(idx.q);                          % n_user x 1
T_0r_0  = z_0(idx.T_0r);                        % scalar
Q_net_0 = z_0(idx.Q_net);                       % scalar
T_0s_0  = u_prev(1);                            % scalar
% T_ir_0 lives in u_prev at positions 1 + n_edges + (1:n_user)
T_ir_0  = u_prev(1 + n_edges + (1:n_user));

%% Per-(consumer, horizon) c_lin/cp coefficients
% c_lin^(i)_h / cp = q_0(i) T_s^(i)_h + (T_s_0(i)-T_ir_0(i)) q^(i)_h
%                    - q_0(i) T_ir^(i)_h - q_0(i)(T_s_0(i)-T_ir_0(i))
% with T_s^(i)_h, q^(i)_h read off z_h and T_ir^(i)_h read off u_h.
% z_h = F_z[h] z_0 + G_z[h] U + zdem_h, so the demand piece zdem_h enters
% the per-consumer offset b(row).
n_S = n_user * Np;
F = zeros(n_S, n_z);
G = zeros(n_S, Np * n_u);
b = zeros(n_S, 1);

for i = 1:n_user
    a_Ts = q_0(i);
    a_Tr = -q_0(i);
    a_q  = T_s_0(i) - T_ir_0(i);
    b_i  = -q_0(i) * (T_s_0(i) - T_ir_0(i));

    % row selectors: pull T_s^(i), q^(i) out of z
    e_Ts = zeros(1, n_z); e_Ts(idx.T_is(i)) = 1;
    e_q  = zeros(1, n_z); e_q (idx.q(i))    = 1;

    for h = 1:Np
        row  = (i-1)*Np + h;
        zrow = (h-1)*n_z + (1:n_z);
        Fz_h = F_z(zrow, :);
        Gz_h = G_z(zrow, :);
        zdem_h = Zdem(zrow);

        % c_lin/cp linear in z_h via a_Ts T_s + a_q q
        F(row, :) = a_Ts * (e_Ts * Fz_h) + a_q * (e_q * Fz_h);
        G(row, :) = a_Ts * (e_Ts * Gz_h) + a_q * (e_q * Gz_h);
        % linear in u_h via T_ir^(i)_h coefficient a_Tr
        u_h_offset = (h-1)*n_u + 1 + n_edges + i;
        G(row, u_h_offset) = G(row, u_h_offset) + a_Tr;
        % constant: op-point term + demand roll-out folded into T_s and q
        b(row) = b_i + a_Ts * (e_Ts * zdem_h) + a_q * (e_q * zdem_h);
    end
end

%% T_s, q, T_0r predictions (linear in z). T_s also carries a demand offset
% Ts_off so the T_r <= T_s constraint in lmpc_solve sees the same
% demand-driven T_s that KPC's demand-in-V extras see. q and T_0r feed the
% energy term, which is linear in U; their constant demand offsets drop
% out of the gradient (verified in lmpc_solve), so we do not carry them.
F_Ts = zeros(n_S, n_z);
G_Ts = zeros(n_S, Np * n_u);
F_q  = zeros(n_S, n_z);
G_q  = zeros(n_S, Np * n_u);
Ts_off = zeros(n_S, 1);
for i = 1:n_user
    e_Ts = zeros(1, n_z); e_Ts(idx.T_is(i)) = 1;
    e_q  = zeros(1, n_z); e_q (idx.q(i))    = 1;
    for h = 1:Np
        row  = (i-1)*Np + h;
        zrow = (h-1)*n_z + (1:n_z);
        Fz_h = F_z(zrow, :);
        Gz_h = G_z(zrow, :);
        zdem_h = Zdem(zrow);
        F_Ts(row, :) = e_Ts * Fz_h;
        G_Ts(row, :) = e_Ts * Gz_h;
        Ts_off(row)  = e_Ts * zdem_h;
        F_q (row, :) = e_q  * Fz_h;
        G_q (row, :) = e_q  * Gz_h;
    end
end

% T_0r prediction: one row per horizon (not per consumer).
F_T0r = zeros(Np, n_z);
G_T0r = zeros(Np, Np * n_u);
e_T0r = zeros(1, n_z); e_T0r(idx.T_0r) = 1;
for h = 1:Np
    zrow = (h-1)*n_z + (1:n_z);
    F_T0r(h, :) = e_T0r * F_z(zrow, :);
    G_T0r(h, :) = e_T0r * G_z(zrow, :);
end

%% Assemble pred struct
pred.F      = F;
pred.G      = G;
pred.b      = b;
pred.F_Ts   = F_Ts;
pred.G_Ts   = G_Ts;
pred.Ts_off = Ts_off;
pred.F_q    = F_q;
pred.G_q    = G_q;
pred.F_T0r  = F_T0r;
pred.G_T0r  = G_T0r;
pred.Np     = Np;
pred.n_user = n_user;
pred.n_z    = n_z;
pred.n_u    = n_u;          % 35: control-only decision vector, same as KPC
pred.horizons = (1:Np)';
pred.meta   = fit_iter.meta;
pred.has_extras = true;
pred.has_mixing = false;    % iterated A,B does not fit the mixing residual
pred.is_lmpc = true;
% Op-point cached for the energy-term linearisation in lmpc_solve.
pred.opp.T_s_0   = T_s_0;
pred.opp.T_ir_0  = T_ir_0;
pred.opp.q_0     = q_0;
pred.opp.T_0s_0  = T_0s_0;
pred.opp.T_0r_0  = T_0r_0;
pred.opp.Q_net_0 = Q_net_0;
end
