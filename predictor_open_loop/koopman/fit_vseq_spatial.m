% NOTE: reproduces the spatial dictionary comparison in dictionary_study.pdf; not on the run_all_tests path.
% FIT_VSEQ_SPATIAL  V_seq fit on the spatial dictionary variants.
%
% Same direct multi-step fit as fit_vseq (one ridge fit per consumer per
% horizon, split regularisation, val-selected lambdas), with the lift
% swapped for spatial_library. Four variants are fitted in one run:
%   none    production 47 lift, refit under this protocol (study reference)
%   full    production 47 + supply-ring block          n_z = 78
%   nolag   as full but without the T_F0 delay block   n_z = 66
%   return  as full plus T_R1..T_R8                    n_z = 86
% Each variant is saved to its own file, vseq_fits_spatial_<variant>.mat;
% the production fit vseq_fits_full.mat is not touched.
%
% Both lambda grids get one extra decade compared to fit_vseq: lam_z because
% the state block is larger, lam_v because the first run picked the 3e5
% boundary for several (consumer, horizon) pairs, in every variant including
% the 47 reference; the val optimum must be interior to the grid. Selection
% stays on val RMSE as before.

clear; clc;
startup;
p   = params();
net = build_plant(p);

variants = {'none', 'full', 'nolag', 'return'};

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

lam_z_grid = [3 30 300];
lam_v_grid = [3 30 300 3000 30000 300000 3000000];
horizons   = 1:p.o1.H_max;   % full sweep, as in fit_vseq
echo_h     = [1 8 16];       % fewer echo horizons than fit_vseq: four variants in one run

for vi = 1:numel(variants)
    variant = variants{vi};

    tr = load_trajs(train_files, p, net, variant);
    va = load_trajs(val_files, p, net, variant);

    n_user = tr{1}.meta.n_user;
    n_z = size(tr{1}.Z, 1);
    n_u = size(tr{1}.U, 1);

    fprintf('\n=== V_seq spatial fit, variant %s: n_z=%d, n_u=%d, %d train + %d val ===\n', ...
        variant, n_z, n_u, numel(tr), numel(va));

    fits = struct();
    fits.horizons = horizons;
    fits.cp     = cell(numel(horizons), n_user);
    fits.dp     = cell(numel(horizons), n_user);
    fits.lam_z  = zeros(numel(horizons), n_user);
    fits.lam_v  = zeros(numel(horizons), n_user);
    fits.meta   = tr{1}.meta;
    fits.n_z    = n_z;
    fits.n_u    = n_u;
    fits.dict   = tr{1}.meta.dict;

    fprintf('%-3s %-3s  %-7s %-7s  %-12s %-12s  %s\n', ...
        'h', 'C', 'lam_z', 'lam_v', 'V_seq c [W]', 'ZOH c [W]', 'KPC/ZOH');
    for hi = 1:numel(horizons)
        h = horizons(hi);
        for i = 1:n_user
            [Phi_tr, Yth_tr]            = stack_horizon(tr, i, h);
            [Phi_va, Yth_va, theta0_va] = stack_horizon(va, i, h);

            [cp, dp, lz_best, lv_best, ~] = ...
                sweep_split_reg(Phi_tr, Yth_tr, Phi_va, Yth_va, n_z, n_u + n_user, h, lam_z_grid, lam_v_grid);

            fits.cp{hi, i}    = cp;
            fits.dp{hi, i}    = dp;
            fits.lam_z(hi, i) = lz_best;
            fits.lam_v(hi, i) = lv_best;

            th_pred = cp' * Phi_va(1:n_z, :) + dp' * Phi_va(n_z+1:end, :);
            c_pred  = p.cp * th_pred;
            c_true  = p.cp * Yth_va;
            c_zoh   = p.cp * theta0_va;

            rmse_kpc = sqrt(mean((c_pred - c_true).^2));
            rmse_zoh = sqrt(mean((c_zoh  - c_true).^2));
            ratio = rmse_kpc / max(rmse_zoh, 1e-9);
            if any(h == echo_h)
                flag = '  ok';
                if rmse_kpc >= rmse_zoh, flag = '  WORSE'; end
                fprintf('%-3d C%d  %-7g %-7g  %-12.1f %-12.1f  %.3f%s\n', ...
                    h, i, lz_best, lv_best, rmse_kpc, rmse_zoh, ratio, flag);
            end
        end
    end

    fits.includes_d_in_V = true;
    fits.d_alignment = 'horizon_major';

    outfile = fullfile(datadir, sprintf('vseq_fits_spatial_%s.mat', variant));
    save(outfile, 'fits');
    fprintf('Saved %s fits to %s\n', variant, outfile);
end

%% --- helpers ---------------------------------------------------------
function trs = load_trajs(files, p, net, variant)
trs = cell(numel(files), 1);
for j = 1:numel(files)
    S = load(files{j});
    traj = S.traj;
    [Z, ~, meta] = spatial_library(traj, p, net, variant);
    U = [traj.T_0s(:)'; traj.r_q; traj.T_ir];
    [~, nm] = fileparts(files{j});
    w = 1; if startsWith(nm, 'op_'), w = p.o1.op_weight; end
    trs{j} = struct('Z',Z,'U',U,'D',traj.d,'theta',traj.theta,'meta',meta,'w',w);
end
end

function [Phi, Yh, Y0] = stack_horizon(trs, consumer_idx, h)
% identical stacking to fit_vseq: Phi = [z_k; u_k..u_{k+h-1}; d_{k+1}..d_{k+h}]
Phi_chunks = cell(numel(trs), 1);
Yh_chunks  = cell(numel(trs), 1);
Y0_chunks  = cell(numel(trs), 1);

for j = 1:numel(trs)
    Z = trs{j}.Z;
    U = trs{j}.U;
    D = trs{j}.D;
    theta = trs{j}.theta(consumer_idx, :);
    N = size(Z, 2);
    n_z = size(Z, 1);
    n_row = n_z + (size(U,1) + size(D,1)) * h;

    ks = max(1, trs{j}.meta.valid_start);
    ke = N - h;
    if ke < ks
        Phi_chunks{j} = zeros(n_row, 0);
        Yh_chunks{j}  = zeros(1, 0);
        Y0_chunks{j}  = zeros(1, 0);
        continue;
    end

    n_s = ke - ks + 1;
    chunk_phi = zeros(n_row, n_s);
    chunk_yh  = zeros(1, n_s);
    chunk_y0  = zeros(1, n_s);
    col = 0;
    for k = ks:ke
        col = col + 1;
        V  = U(:, k:k+h-1);
        Vd = D(:, k+1:k+h);
        chunk_phi(:, col) = [Z(:, k); V(:); Vd(:)];
        chunk_yh(col) = theta(k + h);
        chunk_y0(col) = theta(k);
    end
    sw = sqrt(trs{j}.w);
    Phi_chunks{j} = chunk_phi * sw;
    Yh_chunks{j}  = chunk_yh * sw;
    Y0_chunks{j}  = chunk_y0;
end

Phi = horzcat(Phi_chunks{:});
Yh  = horzcat(Yh_chunks{:});
Y0  = horzcat(Y0_chunks{:});
end

function [cp, dp, lz_best, lv_best, rmse_best] = ...
    sweep_split_reg(Phi_tr, Y_tr, Phi_va, Y_va, n_z, n_u, h, lam_z_grid, lam_v_grid)
% split-lambda ridge, identical to fit_vseq
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
