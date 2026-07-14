% NOTE: reproduces the spatial dictionary comparison in dictionary_study.pdf; not on the run_all_tests path.
% EXTEND_VSEQ_EXTRAS_SPATIAL  Same extras fit as extend_vseq_extras (T_s^i,
% q^i, T_0r and the mixing residuals on top of the theta fits), on the
% spatial full lift. Reads and resaves vseq_fits_spatial_full.mat only; the
% production vseq_fits_full.mat is not touched. Needed so the spatial
% closed-loop leg runs the controller with the same production features
% (T_r <= T_s slack, energy term, mixing slack).

clear; clc;
startup;
p   = params();
net = build_plant(p);
variant = 'full';

datadir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
train_files = cell(p.data.n_train, 1);
for j = 1:p.data.n_train
    train_files{j} = fullfile(datadir, sprintf('train_%02d.mat', j));
end
op = dir(fullfile(datadir, 'op_*.mat'));
for j = 1:numel(op)
    train_files{end+1} = fullfile(datadir, op(j).name); %#ok<SAGROW>
end
val_files = cell(p.data.n_val, 1);
for j = 1:p.data.n_val
    val_files{j} = fullfile(datadir, sprintf('val_%02d.mat', j));
end

fits_file = fullfile(datadir, 'vseq_fits_spatial_full.mat');
S = load(fits_file);
fits = S.fits;
assert(strcmp(fits.dict, 'spatial_full'), 'expected the spatial full fits');

M = build_mixing(net);
tr = load_trajs_full(train_files, p, net, variant, M);
va = load_trajs_full(val_files,   p, net, variant, M);

n_user = tr{1}.meta.n_user;
n_z = size(tr{1}.Z, 1);
n_u = size(tr{1}.U, 1);
horizons = fits.horizons;

% same grids as fit_vseq_spatial (one decade wider than the production fit)
lam_z_grid = [3 30 300];
lam_v_grid = [3 30 300 3000 30000 300000 3000000];

n_mix = M.n_junctions;

fits.Ts_cp   = cell(numel(horizons), n_user);
fits.Ts_dp   = cell(numel(horizons), n_user);
fits.q_cp    = cell(numel(horizons), n_user);
fits.q_dp    = cell(numel(horizons), n_user);
fits.T0r_cp  = cell(numel(horizons), 1);
fits.T0r_dp  = cell(numel(horizons), 1);
fits.mix_cp  = cell(numel(horizons), n_mix);
fits.mix_dp  = cell(numel(horizons), n_mix);
fits.mix_M   = M;

fprintf('\n=== Extending spatial V_seq fits (n_z=%d) with T_s_i, q_i, T_0r, mixing ===\n', n_z);
echo_h = [1 8 16];

for hi = 1:numel(horizons)
    h = horizons(hi);

    for i = 1:n_user
        [Phi_tr, Y_tr] = stack_target(tr, @(t) t.T_is(i,:),    h);
        [Phi_va, Y_va] = stack_target(va, @(t) t.T_is(i,:),    h);
        [cp, dp, ~, ~, rmse_Ts] = sweep_split_reg(Phi_tr, Y_tr, Phi_va, Y_va, n_z, n_u + n_user, h, lam_z_grid, lam_v_grid);
        fits.Ts_cp{hi, i} = cp;
        fits.Ts_dp{hi, i} = dp;

        [Phi_tr, Y_tr] = stack_target(tr, @(t) t.q_users(i,:), h);
        [Phi_va, Y_va] = stack_target(va, @(t) t.q_users(i,:), h);
        [cp, dp, ~, ~, rmse_q] = sweep_split_reg(Phi_tr, Y_tr, Phi_va, Y_va, n_z, n_u + n_user, h, lam_z_grid, lam_v_grid);
        fits.q_cp{hi, i} = cp;
        fits.q_dp{hi, i} = dp;
    end

    [Phi_tr, Y_tr] = stack_target(tr, @(t) t.T_0r,           h);
    [Phi_va, Y_va] = stack_target(va, @(t) t.T_0r,           h);
    [cp, dp, ~, ~, rmse_T0r] = sweep_split_reg(Phi_tr, Y_tr, Phi_va, Y_va, n_z, n_u + n_user, h, lam_z_grid, lam_v_grid);
    fits.T0r_cp{hi} = cp;
    fits.T0r_dp{hi} = dp;

    rmse_mix_max = 0;
    for j = 1:n_mix
        [Phi_tr, Y_tr] = stack_target(tr, @(t) t.mix(j,:), h);
        [Phi_va, Y_va] = stack_target(va, @(t) t.mix(j,:), h);
        [cp, dp, ~, ~, rmse_mix] = sweep_split_reg(Phi_tr, Y_tr, Phi_va, Y_va, n_z, n_u + n_user, h, lam_z_grid, lam_v_grid);
        fits.mix_cp{hi, j} = cp;
        fits.mix_dp{hi, j} = dp;
        rmse_mix_max = max(rmse_mix_max, rmse_mix);
    end

    if any(h == echo_h)
        fprintf('  h=%d:  Ts_C%d RMSE=%.4f K   q_C%d RMSE=%.4f kg/s   T_0r RMSE=%.4f K   mix max RMSE=%.4f K*kg/s\n', ...
            h, n_user, rmse_Ts, n_user, rmse_q, rmse_T0r, rmse_mix_max);
    end
end

save(fits_file, 'fits');
fprintf('\nSaved extended spatial V_seq fits to %s\n', fits_file);


%% helpers
function trs = load_trajs_full(files, p, net, variant, M)
trs = cell(numel(files), 1);
n_mix = M.n_junctions;
for j = 1:numel(files)
    S = load(files{j});
    traj = S.traj;
    [Z, ~, meta] = spatial_library(traj, p, net, variant);
    U = [traj.T_0s(:)'; traj.r_q; traj.T_ir];

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

    [~, nm] = fileparts(files{j});
    w = 1; if startsWith(nm, 'op_'), w = p.o1.op_weight; end
    trs{j} = struct('Z', Z, 'U', U, 'D', traj.d, 'w', w, ...
                    'theta',   traj.theta, ...
                    'T_is',    traj.T_is, ...
                    'q_users', traj.q_users, ...
                    'T_0r',    traj.T_0r, ...
                    'mix',     mix, ...
                    'meta',    meta);
end
end

function [Phi, Yh] = stack_target(trs, target_fn, h)
Phi_chunks = cell(numel(trs), 1);
Yh_chunks  = cell(numel(trs), 1);
for j = 1:numel(trs)
    Z = trs{j}.Z;
    U = trs{j}.U;
    target = target_fn(trs{j});
    D = trs{j}.D;
    N = size(Z, 2);
    n_z = size(Z, 1);
    n_u = size(U, 1);
    n_d = size(D, 1);

    ks = max(1, trs{j}.meta.valid_start);
    ke = N - h;
    if ke < ks
        Phi_chunks{j} = zeros(n_z + (n_u+n_d)*h, 0);
        Yh_chunks{j}  = zeros(1, 0);
        continue;
    end

    n_s = ke - ks + 1;
    cphi = zeros(n_z + (n_u+n_d)*h, n_s);
    cyh  = zeros(1, n_s);
    col = 0;
    for k = ks:ke
        col = col + 1;
        V  = U(:, k:k+h-1);
        Vd = D(:, k+1:k+h);
        cphi(:, col) = [Z(:, k); V(:); Vd(:)];
        cyh(col)     = target(k + h);
    end
    sw = sqrt(trs{j}.w);
    Phi_chunks{j} = cphi * sw;
    Yh_chunks{j}  = cyh * sw;
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
