% EXTEND_VSEQ_EXTRAS  Add V_seq fits for T_s^i, q^i, T_0r and the
% mixing residuals on top of the existing theta fits. Needed for the
% T_r <= T_s slack and the energy term in kpc_v2_solve.

clear; clc;
startup;
p = params();

datadir = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'data_v2_ts900');
train_files = cell(p.data.n_train, 1);
for j = 1:p.data.n_train
    train_files{j} = fullfile(datadir, sprintf('train_%02d.mat', j));
end
val_files = cell(p.data.n_val, 1);
for j = 1:p.data.n_val
    val_files{j} = fullfile(datadir, sprintf('val_%02d.mat', j));
end

% load existing theta fits
fits_file = fullfile(fileparts(mfilename('fullpath')), '..', 'results', 'koopman_v2_ts900', 'vseq_fits_full.mat');
S = load(fits_file);
fits = S.fits;

% reload trajectories with extra signals attached (T_is, q_users, T_0r,
% plus mixing residuals at each junction for Stage 9d)
[net_for_M, ~] = build_plant(p);
M = build_mixing_v2(net_for_M);
tr = load_trajs_full(train_files, p, M);
va = load_trajs_full(val_files,   p, M);

n_user = tr{1}.meta.n_user;
n_z = size(tr{1}.Z, 1);
n_u = size(tr{1}.U, 1);
horizons = fits.horizons;

% smaller lambda grid (these targets are easier to fit than theta).
% lam_v grid extended during T5 to match run_O1_v2_vseq.m's extension.
lam_z_grid = [3 30];
lam_v_grid = [3 30 300 3000 30000 300000];

n_mix = M.n_junctions;

fits.Ts_cp   = cell(numel(horizons), n_user);
fits.Ts_dp   = cell(numel(horizons), n_user);
fits.q_cp    = cell(numel(horizons), n_user);
fits.q_dp    = cell(numel(horizons), n_user);
fits.T0r_cp  = cell(numel(horizons), 1);
fits.T0r_dp  = cell(numel(horizons), 1);
fits.mix_cp  = cell(numel(horizons), n_mix);   % per-junction mixing residual (Stage 9d)
fits.mix_dp  = cell(numel(horizons), n_mix);
fits.mix_M   = M;                              % keep junction structure for build_kpc_v2_matrices

fprintf('\n=== Extending V_seq fits with T_s_i, q_i, T_0r, mixing residuals ===\n');
fprintf('horizons %d..%d, %d users + 1 plant temp + %d mixing junctions\n', ...
    horizons(1), horizons(end), n_user, n_mix);
echo_h = [1 4 8 12 16];

for hi = 1:numel(horizons)
    h = horizons(hi);

    for i = 1:n_user
        % T_s_i target
        [Phi_tr, Y_tr] = stack_target(tr, @(t) t.T_is(i,:),    h);
        [Phi_va, Y_va] = stack_target(va, @(t) t.T_is(i,:),    h);
        [cp, dp, ~, ~, rmse_Ts] = sweep_split_reg(Phi_tr, Y_tr, Phi_va, Y_va, n_z, n_u, h, lam_z_grid, lam_v_grid);
        fits.Ts_cp{hi, i} = cp;
        fits.Ts_dp{hi, i} = dp;

        % q_i target
        [Phi_tr, Y_tr] = stack_target(tr, @(t) t.q_users(i,:), h);
        [Phi_va, Y_va] = stack_target(va, @(t) t.q_users(i,:), h);
        [cp, dp, ~, ~, rmse_q] = sweep_split_reg(Phi_tr, Y_tr, Phi_va, Y_va, n_z, n_u, h, lam_z_grid, lam_v_grid);
        fits.q_cp{hi, i} = cp;
        fits.q_dp{hi, i} = dp;
    end

    % T_0r (single output, no consumer)
    [Phi_tr, Y_tr] = stack_target(tr, @(t) t.T_0r,           h);
    [Phi_va, Y_va] = stack_target(va, @(t) t.T_0r,           h);
    [cp, dp, ~, ~, rmse_T0r] = sweep_split_reg(Phi_tr, Y_tr, Phi_va, Y_va, n_z, n_u, h, lam_z_grid, lam_v_grid);
    fits.T0r_cp{hi} = cp;
    fits.T0r_dp{hi} = dp;

    % per-junction mixing residual (Stage 9d, the tilde g_N(z) = 0 form)
    rmse_mix_max = 0;
    for j = 1:n_mix
        [Phi_tr, Y_tr] = stack_target(tr, @(t) t.mix(j,:), h);
        [Phi_va, Y_va] = stack_target(va, @(t) t.mix(j,:), h);
        [cp, dp, ~, ~, rmse_mix] = sweep_split_reg(Phi_tr, Y_tr, Phi_va, Y_va, n_z, n_u, h, lam_z_grid, lam_v_grid);
        fits.mix_cp{hi, j} = cp;
        fits.mix_dp{hi, j} = dp;
        rmse_mix_max = max(rmse_mix_max, rmse_mix);
    end

    if any(h == echo_h)
        fprintf('  h=%d:  Ts_C2 RMSE=%.4f K   q_C2 RMSE=%.4f kg/s   T_0r RMSE=%.4f K   mix max RMSE=%.4f K*kg/s\n', ...
            h, rmse_Ts, rmse_q, rmse_T0r, rmse_mix_max);
    end
end

save(fits_file, 'fits');
fprintf('\nSaved extended V_seq fits to %s\n', fits_file);


%% helpers (same shape as in run_O1_v2_vseq.m)
function trs = load_trajs_full(files, p, M)
trs = cell(numel(files), 1);
n_mix = M.n_junctions;
for j = 1:numel(files)
    S = load(files{j});
    traj = S.traj;
    [Z, ~, meta] = candidate_library(traj, p);
    U = [traj.T_0s(:)'; traj.r_q; traj.T_ir];

    % per-junction mixing residual: sum_in (T_u - T_v) * q_e per sample
    N = size(Z, 2);
    mix = zeros(n_mix, N);
    for kj = 1:n_mix
        J = M.junction(kj);
        T_v = traj.Tout(J.node_idx, :);
        for ii = 1:numel(J.in_edges)
            T_u = traj.Tout(J.upstream(ii), :);
            q_e = traj.q_edges(J.in_edges(ii), :);
            mix(kj, :) = mix(kj, :) + (T_u - T_v) .* q_e;
        end
    end

    trs{j} = struct('Z', Z, 'U', U, ...
                    'theta',   traj.theta, ...
                    'T_is',    traj.T_is, ...
                    'q_users', traj.q_users, ...
                    'T_0r',    traj.T_0r, ...
                    'mix',     mix, ...
                    'meta',    meta);
end
end

function [Phi, Yh] = stack_target(trs, target_fn, h)
% Generic V_seq stacker: Phi = [z_k; V_seq_k], Yh = target_fn(traj)_{k+h}.
Phi_chunks = cell(numel(trs), 1);
Yh_chunks  = cell(numel(trs), 1);
for j = 1:numel(trs)
    Z = trs{j}.Z;
    U = trs{j}.U;
    target = target_fn(trs{j});
    N = size(Z, 2);
    n_z = size(Z, 1);
    n_u = size(U, 1);

    ks = max(1, trs{j}.meta.valid_start);
    ke = N - h;
    if ke < ks
        Phi_chunks{j} = zeros(n_z + n_u*h, 0);
        Yh_chunks{j}  = zeros(1, 0);
        continue;
    end

    n_s = ke - ks + 1;
    cphi = zeros(n_z + n_u*h, n_s);
    cyh  = zeros(1, n_s);
    col = 0;
    for k = ks:ke
        col = col + 1;
        V = U(:, k:k+h-1);
        cphi(:, col) = [Z(:, k); V(:)];
        cyh(col)     = target(k + h);
    end
    Phi_chunks{j} = cphi;
    Yh_chunks{j}  = cyh;
end
Phi = horzcat(Phi_chunks{:});
Yh  = horzcat(Yh_chunks{:});
end

function [cp, dp, lz_best, lv_best, rmse_best] = ...
    sweep_split_reg(Phi_tr, Y_tr, Phi_va, Y_va, n_z, n_u, h, lam_z_grid, lam_v_grid)
G_tr     = Phi_tr * Phi_tr';
Phi_tr_y = Phi_tr * Y_tr';
rmse_best = inf;
cp = []; dp = []; lz_best = 0; lv_best = 0;
for lz = lam_z_grid
    for lv = lam_v_grid
        Lambda = blkdiag(lz*eye(n_z), lv*eye(n_u*h));
        coef   = (G_tr + Lambda) \ Phi_tr_y;
        err    = coef' * Phi_va - Y_va;
        rmse   = sqrt(mean(err.^2));
        if rmse < rmse_best
            rmse_best = rmse;
            cp = coef(1:n_z);
            dp = coef(n_z+1:end);
            lz_best = lz;
            lv_best = lv;
        end
    end
end
end
