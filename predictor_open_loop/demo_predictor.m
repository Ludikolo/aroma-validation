% DEMO_PREDICTOR  Open-loop demonstration and validation of the
%   47-feature exergy Koopman lift and the V_seq direct multi-step predictor.
%   Loads the saved training + test trajectories, evaluates the
%   production V_seq fit against the iterated A,B baseline and the
%   ZOH null model on the held-out test set, saves the comparison
%   plus the per-consumer parity data, plots the horizon figure, and then
%   verifies five predictor invariants (T1-T5).
%
%   See predictor_validation.pdf, Section 3 (the figures) and Section 4 (the checks).
%
%   The plant lives in shared/ and is verified separately by demo_plant.
%   This demo loads the predictor artefacts that fit_vseq and
%   fit_iterated produce; both fit scripts can be re-run from scratch
%   on a clean clone (saved .mat files ship so the demo runs in ~30 s).

clear; clc;
startup;
p = params();

root        = fileparts(mfilename('fullpath'));
results_dir = fullfile(root, 'results');
figs_dir    = fullfile(root, 'figures');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
if ~exist(figs_dir,    'dir'), mkdir(figs_dir);    end

fprintf('\n=== Predictor demo (V_seq vs iterated A,B vs ZOH) ===\n');
fprintf('Sample rate Ts = %g s, prediction horizon up to %d steps\n', ...
        p.Ts, p.o1.H_max);

%% Load fits and test trajectories
F = load(fullfile(results_dir, 'vseq_fits_full.mat'));
G = load(fullfile(results_dir, 'iterated_AB.mat'));
fits     = F.fits;
fit_iter = G.fit_iter;

horizons = fits.horizons;
n_h      = numel(horizons);
n_user   = size(fits.cp, 2);

test_files = cell(p.data.n_test, 1);
for j = 1:p.data.n_test
    test_files{j} = fullfile(results_dir, sprintf('test_%02d.mat', j));
end

% Pre-build the lifted state Z and the input + demand vectors for each
% test trajectory so the three predictors all start from identical
% inputs (the comparison isolates the predictor's structural choice).
test = cell(numel(test_files), 1);
for j = 1:numel(test_files)
    S = load(test_files{j});
    [Z, ~, meta] = candidate_library(S.traj, p);
    U = [S.traj.T_0s(:)'; S.traj.r_q; S.traj.T_ir];
    test{j} = struct('Z', Z, 'U', U, 'D', S.traj.d, 'theta', S.traj.theta, ...
                     'meta', meta);
end
n_z = size(test{1}.Z, 1);
n_u = size(test{1}.U, 1);
theta_idx = test{1}.meta.idx.theta;
has_d_in_V = isfield(fits, 'includes_d_in_V') && fits.includes_d_in_V;

% Peak demand per consumer is used to normalise RMSE.
d_max = zeros(1, n_user);
for j = 1:numel(test)
    d_max = max(d_max, max(test{j}.D, [], 2)');
end

%% Evaluate all three predictors per (horizon, consumer)
rmse_vseq = zeros(n_h, n_user);
rmse_iter = zeros(n_h, n_user);
rmse_zoh  = zeros(n_h, n_user);
parity_h1.y_true = cell(n_user, 1);
parity_h1.y_pred = cell(n_user, 1);
R2_vseq_h1 = zeros(1, n_user);

for hi = 1:n_h
    h = horizons(hi);
    err_v = cell(n_user, 1);
    err_i = cell(n_user, 1);
    err_z = cell(n_user, 1);
    y_t   = cell(n_user, 1);
    y_v   = cell(n_user, 1);
    for i = 1:n_user
        err_v{i} = []; err_i{i} = []; err_z{i} = [];
        y_t{i}   = []; y_v{i}   = [];
    end

    for j = 1:numel(test)
        T = test{j};
        N = size(T.Z, 2);
        ks = max(1, T.meta.valid_start);
        ke = N - h;
        if ke < ks, continue; end
        for k = ks:ke
            phi_z = T.Z(:, k);

            % iterated one-step rollout from (A, B): apply h times
            z_it = phi_z;
            for s = 0:h-1
                z_it = fit_iter.A * z_it + fit_iter.B * [T.U(:, k+s); T.D(:, k+s+1)];
            end
            theta_iter = z_it(theta_idx);

            % V_seq direct fit at horizon h
            V_u = T.U(:, k : k+h-1);
            V_d = T.D(:, k+1 : k+h);
            for i = 1:n_user
                cp_h = fits.cp{hi, i};
                dp_h = fits.dp{hi, i};
                if has_d_in_V
                    n_v_u = n_u * h;
                    dp_u = dp_h(1 : n_v_u);
                    dp_d = dp_h(n_v_u + 1 : end);
                    theta_v = cp_h(:)' * phi_z + dp_u(:)' * V_u(:) + dp_d(:)' * V_d(:);
                else
                    theta_v = cp_h(:)' * phi_z + dp_h(:)' * V_u(:);
                end
                theta_true = T.theta(i, k + h);
                err_v{i}(end+1) = p.cp * (theta_v       - theta_true);
                err_i{i}(end+1) = p.cp * (theta_iter(i) - theta_true);
                err_z{i}(end+1) = p.cp * (T.theta(i, k) - theta_true);   % ZOH
                if h == 1
                    y_t{i}(end+1) = p.cp * theta_true;
                    y_v{i}(end+1) = p.cp * theta_v;
                end
            end
        end
    end

    for i = 1:n_user
        rmse_vseq(hi, i) = sqrt(mean(err_v{i}.^2));
        rmse_iter(hi, i) = sqrt(mean(err_i{i}.^2));
        rmse_zoh (hi, i) = sqrt(mean(err_z{i}.^2));
        if h == 1
            parity_h1.y_true{i} = y_t{i};
            parity_h1.y_pred{i} = y_v{i};
            ss_res = sum((y_t{i} - y_v{i}).^2);
            ss_tot = sum((y_t{i} - mean(y_t{i})).^2);
            R2_vseq_h1(i) = 1 - ss_res / max(ss_tot, eps);
        end
    end
    fprintf('  h = %2d   suite RMSE  ZOH = %5.0f W   iter = %5.0f W   V_seq = %5.0f W\n', ...
            h, mean(rmse_zoh(hi, :)), mean(rmse_iter(hi, :)), mean(rmse_vseq(hi, :)));
end

% Normalised RMSE (suite mean across consumers, percent of peak demand)
demo.horizons     = horizons;
demo.nrmse_zoh    = 100 * mean(rmse_zoh,  2)' ./ mean(d_max);
demo.nrmse_iter   = 100 * mean(rmse_iter, 2)' ./ mean(d_max);
demo.nrmse_vseq   = 100 * mean(rmse_vseq, 2)' ./ mean(d_max);
demo.rmse_vseq    = rmse_vseq;
demo.rmse_iter    = rmse_iter;
demo.rmse_zoh     = rmse_zoh;
demo.parity_h1    = parity_h1;
demo.R2_vseq_h1   = R2_vseq_h1;
demo.d_max        = d_max;
% provenance: validate_predictor checks n_z and horizons against the fits, so a
% demo.mat generated from different fits is rejected there instead of passing stale
demo.n_z          = fits.n_z;
% checksum over all fitted coefficients: a refit with the same shape but
% different numbers must invalidate this cache too
demo.fit_checksum = sum(cellfun(@(c) sum(c(:)), fits.cp(:))) + ...
                    sum(cellfun(@(c) sum(c(:)), fits.dp(:)));

%% Save
save(fullfile(results_dir, 'demo.mat'), '-struct', 'demo');
fprintf('\nSaved %s/demo.mat\n', results_dir);

%% Plot
plot_predictor_demo(demo, p);

%% Inline validation: five predictor invariants
fprintf('\n=== Predictor invariants ===\n');

% T1: V_seq beats the ZOH null model at every horizon. The point of
% fitting a predictor is to beat predict-current-value.
assert(all(mean(rmse_vseq, 2) <= mean(rmse_zoh, 2)), ...
       'T1 fail: V_seq does not beat ZOH at every horizon');
fprintf('T1 pass: V_seq RMSE <= ZOH RMSE at every horizon (mean across consumers)\n');

% T2: V_seq beats the iterated A,B baseline at the maximum horizon.
% This is the headline claim of the part: direct multi-step does not
% accumulate one-step errors.
hi_max = n_h;
gap_pct = 100 * (mean(rmse_iter(hi_max, :)) - mean(rmse_vseq(hi_max, :))) / mean(rmse_iter(hi_max, :));
assert(mean(rmse_vseq(hi_max, :)) < mean(rmse_iter(hi_max, :)), ...
       'T2 fail: V_seq does not beat iterated at h = %d', horizons(hi_max));
fprintf('T2 pass: V_seq beats iterated at h = %d by %.1f %% RMSE (suite mean)\n', ...
        horizons(hi_max), gap_pct);

% T3: 1-step suite R^2 >= 0.7. A linear lift over a 47-feature exergy
% dictionary should explain at least 70 % of the next-step variance
% on the test set, so 0.7 is a strict lower bound for a usable predictor.
suite_R2_h1 = mean(R2_vseq_h1);
assert(suite_R2_h1 >= 0.7, ...
       'T3 fail: suite R^2 at h = 1 is %.3f, below 0.7', suite_R2_h1);
fprintf('T3 pass: V_seq 1-step suite R^2 = %.3f (>= 0.7)\n', suite_R2_h1);

% T4: 1-step NRMSE <= 30 %. The acceptance
% bar for the V_seq fit: per-consumer RMSE at h = 1 below 30 % of
% peak demand. Stricter than T3 and per-consumer rather than suite.
max_nrmse_h1 = max(100 * rmse_vseq(1, :) ./ d_max);
assert(max_nrmse_h1 <= 30, ...
       'T4 fail: worst-consumer NRMSE at h = 1 is %.1f %%, above 30 %%', max_nrmse_h1);
fprintf('T4 pass: worst-consumer NRMSE at h = 1 = %.1f %% (<= 30 %%)\n', max_nrmse_h1);

% T5: data integrity. The first-order flow dynamics
% q_{k+1} = a q_k + (1 - a) r_q,k must hold on every test
% trajectory at machine precision. If this fails the data
% pipeline is broken and any predictor fit on top is suspect.
fprintf('  data-integrity check ');
max_resid = 0;
[net, ~] = build_plant(p);
Fq_a = net.flow_dyn.a;
Fq_b = net.flow_dyn.b;
for j = 1:p.data.n_test
    S = load(test_files{j});
    Q  = S.traj.q_edges;
    Rq = S.traj.r_q;
    pred_q  = Fq_a .* Q(:, 1:end-1) + Fq_b .* Rq(:, 1:end-1);
    resid_j = max(abs(Q(:, 2:end) - pred_q), [], 'all');
    max_resid = max(max_resid, resid_j);
end
fprintf('(max residual = %.2e kg/s)\n', max_resid);
assert(max_resid < 1e-6, ...
       'T5 fail: flow dynamics residual %.2e on test data', max_resid);
fprintf('T5 pass: q_{k+1} = a q_k + (1 - a) r_q,k on every test trajectory\n');

fprintf('\nAll five predictor invariants verified.\n');
