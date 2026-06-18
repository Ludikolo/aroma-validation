% FIT_VSEQ  Direct multi-step (V_seq) predictor fit with split
% regularisation, per consumer per horizon.
%
% theta^(i)_{k+h} = c_{h,i}^T z_k + d_{h,i}^T V_seq_k
% V_seq_k = [u_k; ...; u_{k+h-1};  d_{k+1}; ...; d_{k+h}],  u = [T_0s; r_q; T_ir]
% The demand forecast d is a known input, so it sits
% in the regressor alongside the applied controls.
%
% Separate lambda_z on c and lambda_v on d, picked by val MSE.
% Evaluation: convert predicted theta to c, compare against ZOH.

clear; clc;
startup;
p = params();

datadir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
train_files = cell(p.data.n_train, 1);
for j = 1:p.data.n_train
    train_files{j} = fullfile(datadir, sprintf('train_%02d.mat', j));
end
% augment with smooth operational days when present (calibrate the operating regime)
op = dir(fullfile(datadir, 'op_*.mat'));
for j = 1:numel(op)
    train_files{end+1} = fullfile(datadir, op(j).name); %#ok<SAGROW>
end
val_files = cell(p.data.n_val, 1);
for j = 1:p.data.n_val
    val_files{j} = fullfile(datadir, sprintf('val_%02d.mat', j));
end

% Load all trajectories once
tr = load_trajs(train_files, p);
va = load_trajs(val_files, p);

n_user = tr{1}.meta.n_user;
n_z = size(tr{1}.Z, 1);
n_u = size(tr{1}.U, 1);

fprintf('\n=== V_seq split-reg fit ===\n');
fprintf('n_z=%d (incl. delays), n_u=%d, %d train + %d val trajs\n', ...
    n_z, n_u, numel(tr), numel(va));

% lambda grids: state and input penalised separately. Trimmed to keep
% the full 1..H_max horizon sweep fast. The lam_v grid was extended
% to 3e5 because the original 3e3 upper bound was binding
% for C3/C4/C5 at long horizon - the ridge optimum for those
% consumers is actually 1-2 orders of magnitude higher (Hastie,
% Tibshirani & Friedman 2009 sec 3.4: extend the cross-validation
% grid until val performance stops improving).
lam_z_grid = [3 30];
lam_v_grid = [3 30 300 3000 30000 300000];

horizons = 1:p.o1.H_max;  % full horizon sweep (= 1:16 at Ts=900); the MPC horizon is picked from this

% per-(consumer, horizon) fit + report. Save the coefficients for the
% MPC: it will stack F = [c_1; c_2; ...; c_Np] and the input
% mapping G derived from d_h.
fits.horizons = horizons;
fits.cp     = cell(numel(horizons), n_user);
fits.dp     = cell(numel(horizons), n_user);
fits.lam_z  = zeros(numel(horizons), n_user);
fits.lam_v  = zeros(numel(horizons), n_user);
fits.meta   = tr{1}.meta;
fits.n_z    = n_z;
fits.n_u    = n_u;

% only echo the headline horizons (less terminal noise during the full sweep)
echo_h = [1 4 8 12 16];
fprintf('\n%-3s %-3s  %-7s %-7s  %-12s %-12s  %s\n', ...
    'h', 'C', 'lam_z', 'lam_v', 'V_seq c [W]', 'ZOH c [W]', 'KPC/ZOH');
% THEORY (direct multi-step): one independent fit per (consumer, horizon), no rollout, so errors do not compound
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

% V_seq includes the demand forecast over the horizon as a regressor block,
% so flag it for the controller to fold the demand contribution G_d * D
% into the prediction.
fits.includes_d_in_V = true;
% documents the demand-block ordering in the fit (horizon-major); the actual re-index to consumer-major order is done in build_kpc_v2_matrices.m
fits.d_alignment = 'horizon_major';

% save the full 1..Np horizon sweep for the controller
outdir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end
save(fullfile(outdir, 'vseq_fits_full.mat'), 'fits');
fprintf('\nSaved V_seq fits to %s\n', fullfile(outdir, 'vseq_fits_full.mat'));

%% --- helpers ---------------------------------------------------------
function trs = load_trajs(files, p)
trs = cell(numel(files), 1);
for j = 1:numel(files)
    S = load(files{j});
    traj = S.traj;
    [Z, ~, meta] = candidate_library(traj, p);
    U = [traj.T_0s(:)'; traj.r_q; traj.T_ir];
    [~, nm] = fileparts(files{j});
    w = 1; if startsWith(nm, 'op_'), w = p.o1.op_weight; end
    trs{j} = struct('Z',Z,'U',U,'D',traj.d,'theta',traj.theta,'meta',meta,'w',w);
end
end

function [Phi, Yh, Y0] = stack_horizon(trs, consumer_idx, h)
% Phi: [z_k; u_k..u_{k+h-1}; d_{k+1}..d_{k+h}]  x N_samples
% The demand forecast over the horizon is part of the regressor:
% demand is a known input.
% Yh : theta(consumer_idx)_{k+h}                     1 x N_samples
% Y0 : theta(consumer_idx)_{k}     (for ZOH baseline) 1 x N_samples
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
        V  = U(:, k:k+h-1);          % the h applied inputs
        Vd = D(:, k+1:k+h);          % the demand forecast over the horizon
        chunk_phi(:, col) = [Z(:, k); V(:); Vd(:)];
        chunk_yh(col) = theta(k + h);
        chunk_y0(col) = theta(k);
    end
    sw = sqrt(trs{j}.w);          % operational days weigh more (sqrt in LS)
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
% Ridge with split lambda on (state, input) blocks. Returns the (cp, dp)
% that minimises val RMSE on the target signal.

G_tr     = Phi_tr * Phi_tr';   % Gram matrix on the training regressors
Phi_tr_y = Phi_tr * Y_tr';

rmse_best = inf;
cp = []; dp = []; lz_best = 0; lv_best = 0;
for lz = lam_z_grid
    for lv = lam_v_grid
        % THEORY (split ridge): separate lambda on the state block (n_z columns) and the
        % input block ((n_u+n_user)*h columns = the applied controls plus the demand
        % forecast over the horizon), picked on validation RMSE
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
