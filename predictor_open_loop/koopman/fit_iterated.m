% FIT_ITERATED  Fit the iterated one-step predictor
%
%   z_{k+1} = A z_k + B [u_k; d_{k+1}],   u = [T_0s; r_q; T_ir]
%
% The input carries the next-step demand forecast d_{k+1} as well, so the
% iterated model sees the same demand information as the direct V_seq fit
% and the comparison isolates the predictor structure, not the inputs.
%
% via ridge regression on the same training data as the V_seq fit. The
% lambda is picked by the 1-step val RMSE on theta (the standard one-step
% EDMDc fit). The iterated model is then rolled out h steps, where its
% error accumulates; that rollout gap against the direct V_seq fit is the
% multi-step result, not a crippled one-step.
%
% This is the "iterated (A, B)" leg of the three-way comparison
% (V_seq vs iterated vs ZOH); written as a fresh fit rather than
% reusing an old helper because the dictionary, input vector, and
% sample rate are all different.

clear; clc;
startup;
p = params();

datadir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
train_files = cell(p.data.n_train, 1);
for j = 1:p.data.n_train
    train_files{j} = fullfile(datadir, sprintf('train_%02d.mat', j));
end
% same augmented training set as the V_seq fit, for a fair comparison
op = dir(fullfile(datadir, 'op_*.mat'));
for j = 1:numel(op)
    train_files{end+1} = fullfile(datadir, op(j).name); %#ok<SAGROW>
end
val_files = cell(p.data.n_val, 1);
for j = 1:p.data.n_val
    val_files{j} = fullfile(datadir, sprintf('val_%02d.mat', j));
end

% Load every trajectory and build (z_k, u_k, z_{k+1}) tuples.
% Same dictionary builder as the V_seq fit so the iterated and direct
% predictors see the same lift; the comparison isolates the predictor
% structure (iterated vs direct multi-step) rather than the lift.
[Z_tr, U_tr, Zn_tr, meta] = stack_pairs(train_files, p);
[~, ~, ~, ~, val_data] = stack_pairs(val_files, p);
n_z = size(Z_tr, 1);
n_u = size(U_tr, 1);
fprintf('Iterated (A,B) fit: n_z = %d, n_u = %d, %d training pairs\n', ...
        n_z, n_u, size(Z_tr, 2));
theta_idx = meta.idx.theta;     % 1 x n_user, where theta_C_i sits in z

% Ridge sweep. Extended grid (log-spaced) to make sure the validation
% optimum is interior to the grid; an earlier [3 30 300 3000] grid
% capped the search at lambda = 3000, a boundary choice where RMSE was
% still decreasing. Pick by 1-step val RMSE on theta over
% the 5 consumers.
lambda_grid = [1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 300000];
rmse_per_lam = zeros(size(lambda_grid));

X = [Z_tr; U_tr];
G_tr     = X * X';
ZnX_tr   = Zn_tr * X';

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
    % only keep stable fits (rho(A) < 1): an unstable rollout blows up at long
    % horizon, which would make iterated an unfairly weak baseline
    if rho_li < 1 && rmse_per_lam(li) < best.rmse
        best.rmse   = rmse_per_lam(li);
        best.lambda = lam;
        best.A      = A;
        best.B      = B;
    end
end
fprintf('\nBest lambda = %g (1-step val theta RMSE = %.4f K*kg/s)\n', ...
        best.lambda, best.rmse);

% Report rho(A) for the chosen fit. The selection above already rejects
% rho(A) >= 1, so the rollout is stable and the comparison is fair;
% this just records the value that ended up being used.
rho_A = max(abs(eig(best.A)));
fprintf('rho(A) at chosen lambda = %.4f (< 1 means iterated rollout is stable)\n', rho_A);

% Save in the same directory as V_seq fits so the comparison plot can load both.
outdir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end
fit_iter = struct();
fit_iter.A          = best.A;
fit_iter.B          = best.B;
fit_iter.lambda     = best.lambda;
fit_iter.rho_A      = rho_A;
fit_iter.lambda_grid = lambda_grid;
fit_iter.rmse_per_lambda = rmse_per_lam;
fit_iter.n_z = n_z; fit_iter.n_u = n_u;
fit_iter.meta = meta;
save(fullfile(outdir, 'iterated_AB.mat'), 'fit_iter');
fprintf('\nSaved iterated_AB.mat (lambda*=%g, rho(A)=%.4f) to %s\n', ...
        best.lambda, rho_A, outdir);


%% --- helpers ---------------------------------------------------------
function [Z, U, Zn, meta_first, traj_data] = stack_pairs(files, p)
% Aggregate (z_k, u_k, z_{k+1}) over all training trajectories. Also
% return the per-traj structures for val rollout.
n_traj = numel(files);
traj_data = cell(n_traj, 1);
chunks_z  = cell(n_traj, 1);
chunks_u  = cell(n_traj, 1);
chunks_zn = cell(n_traj, 1);
for j = 1:n_traj
    S = load(files{j});
    traj = S.traj;
    [Zj, ~, m] = candidate_library(traj, p);
    Uj = [traj.T_0s(:)'; traj.r_q; traj.T_ir];

    N  = size(Zj, 2);
    % the input also carries the next-step demand forecast, so the iterated
    % model sees the same demand the V_seq fit does
    Wj = [Uj(:, 1:N-1); traj.d(:, 2:N)];
    ks = max(1, m.valid_start);
    ke = N - 1;
    if ke < ks
        chunks_z{j} = []; chunks_u{j} = []; chunks_zn{j} = [];
        traj_data{j} = struct('Z', Zj, 'U', Wj, 'meta', m);
        continue;
    end
    % operational days weigh more, same as the V_seq fit, so the comparison
    % is on equal footing (sqrt weight in least squares)
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
% For each val traj, roll the iterated predictor forward H steps from
% every valid sample and compute theta RMSE per consumer; return the
% mean across consumers and trajectories.
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
    n_z = size(Z, 1);

    for k = ks:ke
        z = Z(:, k);
        for s = 0:H-1
            z = A * z + B * U(:, k + s);
        end
        % predicted theta at k+H
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
