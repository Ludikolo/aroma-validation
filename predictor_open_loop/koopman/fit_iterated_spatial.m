% NOTE: reproduces the spatial dictionary comparison in dictionary_study.pdf; not on the run_all_tests path.
% FIT_ITERATED_SPATIAL  Iterated one-step fit on the spatial lift (full variant).
%
% Same fit as fit_iterated (ridge on one-step pairs, lambda picked by the
% 1-step val RMSE on theta, stability filter rho(A) < 1), with the lift
% swapped for spatial_library('full'). Saved to iterated_AB_spatial.mat;
% the production iterated_AB.mat is not touched. Needed for the iterated
% leg of the dictionary study, open and closed loop.

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

[Z_tr, U_tr, Zn_tr, meta] = stack_pairs(train_files, p, net, variant);
[~, ~, ~, ~, val_data] = stack_pairs(val_files, p, net, variant);
n_z = size(Z_tr, 1);
n_u = size(U_tr, 1);
fprintf('Iterated (A,B) spatial fit: n_z = %d, n_u = %d, %d training pairs\n', ...
        n_z, n_u, size(Z_tr, 2));
theta_idx = meta.idx.theta;

lambda_grid = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 300000];
rmse_per_lam = zeros(size(lambda_grid));

X = [Z_tr; U_tr];
G_tr   = X * X';
ZnX_tr = Zn_tr * X';

best.rmse = inf;
for li = 1:numel(lambda_grid)
    lam = lambda_grid(li);
    AB = ZnX_tr / (G_tr + lam * eye(n_z + n_u));
    A  = AB(:, 1:n_z);
    B  = AB(:, n_z+1:end);
    rho_li = max(abs(eig(A)));

    rmse_per_lam(li) = mean_theta_rollout_rmse(val_data, A, B, theta_idx, 1);
    fprintf('  lambda = %5g  ->  1-step val theta RMSE = %.4f K*kg/s   rho(A) = %.4f\n', ...
        lam, rmse_per_lam(li), rho_li);
    if rho_li < 1 && rmse_per_lam(li) < best.rmse
        best.rmse   = rmse_per_lam(li);
        best.lambda = lam;
        best.A      = A;
        best.B      = B;
    end
end
fprintf('\nBest lambda = %g (1-step val theta RMSE = %.4f K*kg/s)\n', ...
        best.lambda, best.rmse);

rho_A = max(abs(eig(best.A)));
fprintf('rho(A) at chosen lambda = %.4f\n', rho_A);

fit_iter = struct();
fit_iter.A          = best.A;
fit_iter.B          = best.B;
fit_iter.lambda     = best.lambda;
fit_iter.rho_A      = rho_A;
fit_iter.lambda_grid = lambda_grid;
fit_iter.rmse_per_lambda = rmse_per_lam;
fit_iter.n_z = n_z; fit_iter.n_u = n_u;
fit_iter.meta = meta;
fit_iter.dict = meta.dict;
save(fullfile(datadir, 'iterated_AB_spatial.mat'), 'fit_iter');
fprintf('\nSaved iterated_AB_spatial.mat (lambda*=%g, rho(A)=%.4f)\n', best.lambda, rho_A);


%% --- helpers ---------------------------------------------------------
function [Z, U, Zn, meta_first, traj_data] = stack_pairs(files, p, net, variant)
% identical pair stacking to fit_iterated, with the spatial lift
n_traj = numel(files);
traj_data = cell(n_traj, 1);
chunks_z  = cell(n_traj, 1);
chunks_u  = cell(n_traj, 1);
chunks_zn = cell(n_traj, 1);
for j = 1:n_traj
    S = load(files{j});
    traj = S.traj;
    [Zj, ~, m] = spatial_library(traj, p, net, variant);
    Uj = [traj.T_0s(:)'; traj.r_q; traj.T_ir];

    N  = size(Zj, 2);
    Wj = [Uj(:, 1:N-1); traj.d(:, 2:N)];
    ks = max(1, m.valid_start);
    ke = N - 1;
    if ke < ks
        chunks_z{j} = []; chunks_u{j} = []; chunks_zn{j} = [];
        traj_data{j} = struct('Z', Zj, 'U', Wj, 'meta', m);
        continue;
    end
    [~, nm] = fileparts(files{j});
    sw = 1; if startsWith(nm, 'op_'), sw = sqrt(p.o1.op_weight); end
    chunks_z{j}  = Zj(:, ks:ke)     * sw;
    chunks_u{j}  = Wj(:, ks:ke)     * sw;
    chunks_zn{j} = Zj(:, ks+1:ke+1) * sw;
    traj_data{j} = struct('Z', Zj, 'U', Wj, 'meta', m);
    if j == 1, meta_first = m; end
end
Z  = horzcat(chunks_z{:});
U  = horzcat(chunks_u{:});
Zn = horzcat(chunks_zn{:});
end

function rmse_mean = mean_theta_rollout_rmse(traj_data, A, B, theta_idx, H)
% identical rollout metric to fit_iterated
n_user = numel(theta_idx);
err_sq_sum = zeros(n_user, 1);
n_count    = 0;
for j = 1:numel(traj_data)
    if isempty(traj_data{j}), continue; end
    Z = traj_data{j}.Z;
    U = traj_data{j}.U;
    m = traj_data{j}.meta;
    N  = size(Z, 2);
    ks = max(1, m.valid_start);
    ke = N - H;
    if ke < ks, continue; end

    for k = ks:ke
        z = Z(:, k);
        for s = 0:H-1
            z = A * z + B * U(:, k + s);
        end
        theta_pred = z(theta_idx);
        theta_true = Z(theta_idx, k + H);
        err_sq_sum = err_sq_sum + (theta_pred - theta_true).^2;
        n_count    = n_count + 1;
    end
end
if n_count == 0
    rmse_mean = inf;
    return;
end
rmse_per_user = sqrt(err_sq_sum / n_count);
rmse_mean = mean(rmse_per_user);
end
