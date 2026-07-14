function pred = build_kpc_v2_matrices(fits)
% BUILD_KPC_V2_MATRICES  Stack per-(consumer, horizon) V_seq fits into
% MPC prediction matrices F, G (and G_d in d-in-V mode).
%
% Two fit modes are supported transparently:
%
%   Legacy:    theta^(i)_h = c_{h,i}^T z_0 + d_{h,i}^T V_seq_h
%              V_seq_h = [u_0; ...; u_{h-1}],   dp length n_u*h.
%
%   d-in-V:    V_seq_h = [u_0; ...; u_{h-1}; d_1; ...; d_h],
%              dp length (n_u + n_user)*h, fits.includes_d_in_V = true.
%              Adds pred.G_d (n_user*Np x Np*n_user) so the QP can
%              compute c_pred(U, d_seq) = cp*(F z + G U + G_d D).
%
% Stacking is consumer-major: all Np horizons of C1, then C2, ...
%
% THEORY index:
%   line 38: F
%   line 41: G

n_user = size(fits.cp, 2);
Np     = numel(fits.horizons);
n_z    = fits.n_z;
n_u    = fits.n_u;
has_d  = isfield(fits, 'includes_d_in_V') && fits.includes_d_in_V;

F = zeros(n_user * Np, n_z);
G = zeros(n_user * Np, Np * n_u);
if has_d
    G_d = zeros(n_user * Np, Np * n_user);
end

for i = 1:n_user
    for hi = 1:Np
        h   = fits.horizons(hi);
        row = (i - 1) * Np + hi;

        % THEORY (F): row (i,h) = the state weights cp of that cell, acts on z0
        F(row, :) = fits.cp{hi, i}';

        % THEORY (G): row (i,h) = the input weights dp, only the first h inputs filled, block lower triangular
        dh = fits.dp{hi, i};
        if has_d
            expected = (n_u + n_user) * h;
            if numel(dh) ~= expected
                error('build_kpc_v2_matrices: d-in-V d_{%d,C%d} has length %d, expected %d', ...
                    h, i, numel(dh), expected);
            end
            dh_u = dh(1 : n_u * h);
            dh_d = dh(n_u * h + 1 : end);
            G  (row, 1:numel(dh_u)) = dh_u';
            % dh_d is horizon-major in fit order: [d_1(C1..Cn); d_2(C1..Cn); ...]
            % D_vec at solve time is consumer-major: [C1_h1..C1_hNp; C2_h1..C2_hNp; ...]
            % Re-index to put each fit coefficient at the right D_vec slot.
            for tt = 0 : h - 1
                for cc = 1 : n_user
                    src_idx = tt * n_user + cc;            % position in dh_d
                    dst_idx = (cc - 1) * Np + tt + 1;       % position in D_vec
                    G_d(row, dst_idx) = dh_d(src_idx);
                end
            end
        else
            if numel(dh) ~= n_u * h
                error('build_kpc_v2_matrices: d_{%d,C%d} has length %d, expected %d', ...
                    h, i, numel(dh), n_u * h);
            end
            G(row, 1:numel(dh)) = dh';
        end
    end
end

pred.F        = F;
pred.G        = G;
pred.Np       = Np;
pred.n_user   = n_user;
pred.n_z      = n_z;
pred.n_u      = n_u;
pred.horizons = fits.horizons(:);
pred.meta     = fits.meta;          % keep idx.T_0r, idx.Q_net etc. handy
pred.has_demand_in_V = has_d;
if has_d
    pred.G_d = G_d;
end

% extras: T_s, q^(i), T_0r predictions (from extend_vseq_extras)
% Optional: if absent the QP just skips T_r <= T_s and the energy term.
if isfield(fits, 'Ts_cp') && isfield(fits, 'q_cp') && isfield(fits, 'T0r_cp')
    F_Ts = zeros(n_user * Np, n_z);
    G_Ts = zeros(n_user * Np, Np * n_u);
    F_q  = zeros(n_user * Np, n_z);
    G_q  = zeros(n_user * Np, Np * n_u);
    if has_d
        G_d_Ts = zeros(n_user * Np, Np * n_user);
        G_d_q  = zeros(n_user * Np, Np * n_user);
    end
    for i = 1:n_user
        for hi = 1:Np
            h   = fits.horizons(hi);
            row = (i - 1) * Np + hi;
            F_Ts(row, :) = fits.Ts_cp{hi, i}';
            F_q (row, :) = fits.q_cp {hi, i}';
            d_Ts = fits.Ts_dp{hi, i};
            d_q  = fits.q_dp {hi, i};
            if has_d
                d_Ts_u = d_Ts(1 : n_u * h);
                d_Ts_d = d_Ts(n_u * h + 1 : end);
                d_q_u  = d_q (1 : n_u * h);
                d_q_d  = d_q (n_u * h + 1 : end);
                G_Ts(row, 1:numel(d_Ts_u)) = d_Ts_u';
                G_q (row, 1:numel(d_q_u )) = d_q_u';
                for tt = 0 : h - 1
                    for cc = 1 : n_user
                        src_idx = tt * n_user + cc;
                        dst_idx = (cc - 1) * Np + tt + 1;
                        G_d_Ts(row, dst_idx) = d_Ts_d(src_idx);
                        G_d_q (row, dst_idx) = d_q_d (src_idx);
                    end
                end
            else
                G_Ts(row, 1:numel(d_Ts)) = d_Ts';
                G_q (row, 1:numel(d_q )) = d_q';
            end
        end
    end

    F_T0r = zeros(Np, n_z);
    G_T0r = zeros(Np, Np * n_u);
    if has_d, G_d_T0r = zeros(Np, Np * n_user); end
    for hi = 1:Np
        h = fits.horizons(hi);
        F_T0r(hi, :) = fits.T0r_cp{hi}';
        d_T0r = fits.T0r_dp{hi};
        if has_d
            d_T0r_u = d_T0r(1 : n_u * h);
            d_T0r_d = d_T0r(n_u * h + 1 : end);
            G_T0r(hi, 1:numel(d_T0r_u)) = d_T0r_u';
            for tt = 0 : h - 1
                for cc = 1 : n_user
                    src_idx = tt * n_user + cc;
                    dst_idx = (cc - 1) * Np + tt + 1;
                    G_d_T0r(hi, dst_idx) = d_T0r_d(src_idx);
                end
            end
        else
            G_T0r(hi, 1:numel(d_T0r)) = d_T0r';
        end
    end

    pred.F_Ts  = F_Ts;
    pred.G_Ts  = G_Ts;
    pred.F_q   = F_q;
    pred.G_q   = G_q;
    pred.F_T0r = F_T0r;
    pred.G_T0r = G_T0r;
    pred.has_extras = true;
    if has_d
        pred.G_d_Ts  = G_d_Ts;
        pred.G_d_q   = G_d_q;
        pred.G_d_T0r = G_d_T0r;
    end
else
    pred.has_extras = false;
end

% mixing residual prediction matrices
% Per-junction residual mix^(j)_h = c_{h,j}^T z_0 + d_{h,j}^T V_seq_h.
% Used as a soft equality |mix^(j)_h| <= s_mix in kpc_v2_solve, since
% the plant only respects mass conservation asymptotically.
if isfield(fits, 'mix_cp') && isfield(fits, 'mix_M')
    n_mix = fits.mix_M.n_junctions;
    F_mix = zeros(n_mix * Np, n_z);
    G_mix = zeros(n_mix * Np, Np * n_u);
    if has_d, G_d_mix = zeros(n_mix * Np, Np * n_user); end
    for j = 1:n_mix
        for hi = 1:Np
            h   = fits.horizons(hi);
            row = (j - 1) * Np + hi;
            F_mix(row, :) = fits.mix_cp{hi, j}';
            d_m = fits.mix_dp{hi, j};
            if has_d
                d_m_u = d_m(1 : n_u * h);
                d_m_d = d_m(n_u * h + 1 : end);
                G_mix(row, 1:numel(d_m_u)) = d_m_u';
                for tt = 0 : h - 1
                    for cc = 1 : n_user
                        src_idx = tt * n_user + cc;
                        dst_idx = (cc - 1) * Np + tt + 1;
                        G_d_mix(row, dst_idx) = d_m_d(src_idx);
                    end
                end
            else
                G_mix(row, 1:numel(d_m)) = d_m';
            end
        end
    end
    pred.F_mix     = F_mix;
    pred.G_mix     = G_mix;
    pred.n_mix     = n_mix;
    pred.has_mixing = true;
    if has_d, pred.G_d_mix = G_d_mix; end
else
    pred.has_mixing = false;
end

end
