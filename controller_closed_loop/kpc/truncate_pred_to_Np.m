function pred_new = truncate_pred_to_Np(pred, Np_new)
% TRUNCATE_PRED_TO_NP  Slice an MPC predictor to a shorter horizon.
% Used by the (Np, Nc) sweep and any orchestrator that fits once at
% the maximum horizon and re-uses the same predictor at smaller Np.
%
% Stacking conventions (must match build_kpc_v2_matrices):
%   F, G, F_Ts, G_Ts, F_q, G_q : rows are consumer-major (Np rows per consumer)
%   F_mix, G_mix               : rows are junction-major (Np rows per junction)
%   F_T0r, G_T0r               : single Np-row block
%   G* columns                 : Np u-blocks of size n_u (or Np demand-blocks of size n_user for G_d*)
%   G_d* demand columns        : consumer-major (Np columns per consumer)

assert(Np_new >= 1 && Np_new <= pred.Np, ...
    'truncate_pred_to_Np: Np_new must be in 1..pred.Np');
if Np_new == pred.Np
    pred_new = pred;
    return;
end

Np     = pred.Np;
n_user = pred.n_user;
n_u    = pred.n_u;

% row selectors for the consumer- and junction-major blocks
rows_consumer = block_rows(n_user, Np, Np_new);
col_U_new     = 1 : (Np_new * n_u);

pred_new          = pred;
pred_new.Np       = Np_new;
pred_new.horizons = pred.horizons(1:Np_new);

pred_new.F = pred.F(rows_consumer, :);
pred_new.G = pred.G(rows_consumer, col_U_new);

if pred.has_demand_in_V
    cols_D_new   = demand_cols(n_user, Np, Np_new);
    pred_new.G_d = pred.G_d(rows_consumer, cols_D_new);
end

if pred.has_extras
    pred_new.F_Ts  = pred.F_Ts (rows_consumer, :);
    pred_new.G_Ts  = pred.G_Ts (rows_consumer, col_U_new);
    pred_new.F_q   = pred.F_q  (rows_consumer, :);
    pred_new.G_q   = pred.G_q  (rows_consumer, col_U_new);
    pred_new.F_T0r = pred.F_T0r(1:Np_new, :);
    pred_new.G_T0r = pred.G_T0r(1:Np_new, col_U_new);
    if pred.has_demand_in_V
        cols_D_new       = demand_cols(n_user, Np, Np_new);
        pred_new.G_d_Ts  = pred.G_d_Ts (rows_consumer, cols_D_new);
        pred_new.G_d_q   = pred.G_d_q  (rows_consumer, cols_D_new);
        pred_new.G_d_T0r = pred.G_d_T0r(1:Np_new,      cols_D_new);
    end
end

if pred.has_mixing
    rows_mix       = block_rows(pred.n_mix, Np, Np_new);
    pred_new.F_mix = pred.F_mix(rows_mix, :);
    pred_new.G_mix = pred.G_mix(rows_mix, col_U_new);
    if pred.has_demand_in_V
        cols_D_new      = demand_cols(n_user, Np, Np_new);
        pred_new.G_d_mix = pred.G_d_mix(rows_mix, cols_D_new);
    end
end

end


function rows = block_rows(n_blocks, Np_old, Np_new)
% For a matrix stacked as n_blocks blocks of Np_old rows, pick the first
% Np_new rows from every block. Returns a column vector of row indices.
rows = zeros(n_blocks * Np_new, 1);
for b = 1:n_blocks
    dst = (b-1)*Np_new + (1:Np_new);
    src = (b-1)*Np_old + (1:Np_new);
    rows(dst) = src;
end
end


function cols = demand_cols(n_user, Np_old, Np_new)
% G_d demand columns are stacked consumer-major (Np_old columns per
% consumer). Pick the first Np_new columns from every consumer block.
cols = block_rows(n_user, Np_old, Np_new);
end
