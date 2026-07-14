% VALIDATE_PREDICTOR  Check the V_seq direct multi-step predictor
%   against five invariants on the saved demo trajectories.
%
%   Run demo_predictor first: this script reads its demo.mat and would
%   otherwise silently check a stale cache, so it checks below that
%   demo.mat records the same n_z and horizons as the fits and errors
%   if not. It does not refit; it operates on the saved data + the
%   committed fit artefacts and prints a per-block summary.
%
%   T1  V_seq RMSE <= ZOH RMSE at every horizon (the predictor beats
%       predict-current-value).
%   T2  V_seq RMSE < iterated A,B RMSE at the longest horizon (the
%       direct multi-step claim).
%   T3  Suite R^2 at h = 1 >= 0.7 (the lift carries enough information
%       to be a usable controller predictor), and the weakest single
%       consumer reaches R^2 >= 0.5 so the suite mean hides no one.
%   T4  Worst-consumer NRMSE at h = 1 <= 30 % of peak demand (the
%       per-consumer acceptance bar).
%   T5  Flow dynamics on every test trajectory satisfy
%       q_{k+1} = a q_k + (1 - a) r_q,k at machine precision (data
%       pipeline is internally consistent).

clear; clc;
startup;
p = params();

results_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                       'predictor_open_loop', 'results');

% guard against a stale demo cache: demo.mat records the n_z and horizons of
% the fits it was generated from, so a demo produced against different fits
% fails loudly here (file dates are useless for this check: on a fresh clone
% the checkout order sets the mtimes, not the generation order)
S  = load(fullfile(results_dir, 'demo.mat'));
Fv = load(fullfile(results_dir, 'vseq_fits_full.mat'));
chk = sum(cellfun(@(c) sum(c(:)), Fv.fits.cp(:))) + ...
      sum(cellfun(@(c) sum(c(:)), Fv.fits.dp(:)));
assert(isfield(S, 'fit_checksum') && abs(S.fit_checksum - chk) < 1e-9 * max(1, abs(chk)), ...
    'stale demo cache: fit coefficients changed since demo.mat was generated');
assert(isfield(S, 'n_z') && S.n_z == Fv.fits.n_z && isequal(S.horizons, Fv.fits.horizons), ...
       'demo.mat provenance does not match vseq_fits_full.mat; run demo_predictor to refresh it');

fprintf('\n=== Predictor validation: five invariants ===\n');
fprintf('provenance ok: demo.mat matches the fits (n_z = %d, %d horizons)\n', ...
        S.n_z, numel(S.horizons));

%% T1 V_seq beats ZOH at every horizon
assert(all(mean(S.rmse_vseq, 2) <= mean(S.rmse_zoh, 2)), ...
       'T1 fail: V_seq does not beat ZOH at every horizon');
fprintf('T1 ok: V_seq RMSE <= ZOH RMSE at every horizon\n');

%% T2 V_seq beats iterated A,B at the longest horizon
hi_max = numel(S.horizons);
gap_pct = 100 * (mean(S.rmse_iter(hi_max, :)) - mean(S.rmse_vseq(hi_max, :))) ...
                / mean(S.rmse_iter(hi_max, :));
assert(mean(S.rmse_vseq(hi_max, :)) < mean(S.rmse_iter(hi_max, :)), ...
       'T2 fail: V_seq does not beat iterated at h = %d', S.horizons(hi_max));
fprintf('T2 ok: V_seq beats iterated at h = %d by %.1f %% RMSE (suite mean)\n', ...
        S.horizons(hi_max), gap_pct);

%% T3 Suite R^2 at h = 1 >= 0.7, weakest consumer >= 0.5
suite_R2 = mean(S.R2_vseq_h1);
assert(suite_R2 >= 0.7, 'T3 fail: suite R^2 at h = 1 is %.3f, below 0.7', suite_R2);
fprintf('T3 ok: V_seq 1-step suite R^2 = %.3f\n', suite_R2);

% per-consumer floor: the suite mean can hide one badly fitted consumer, and
% the parity figure prints the per-consumer values, so pin the weakest one too
min_R2 = min(S.R2_vseq_h1);
assert(min_R2 >= 0.5, ...
       'T3 floor fail: weakest consumer R^2 at h = 1 is %.3f, below 0.5', min_R2);
fprintf('T3 floor ok: weakest consumer 1-step R^2 = %.3f\n', min_R2);

%% T4 Worst-consumer NRMSE at h = 1 <= 30 %
max_nrmse = max(100 * S.rmse_vseq(1, :) ./ S.d_max);
assert(max_nrmse <= 30, ...
       'T4 fail: worst-consumer NRMSE at h = 1 is %.1f %%, above 30 %%', max_nrmse);
fprintf('T4 ok: worst-consumer NRMSE at h = 1 = %.1f %%\n', max_nrmse);

%% T5 Flow-dynamics consistency on every test trajectory
[net, ~] = build_plant(p);
Fq_a = net.flow_dyn.a;
Fq_b = net.flow_dyn.b;
max_resid = 0;
for j = 1:p.data.n_test
    T = load(fullfile(results_dir, sprintf('test_%02d.mat', j)));
    Q  = T.traj.q_edges;
    Rq = T.traj.r_q;
    pred_q = Fq_a .* Q(:, 1:end-1) + Fq_b .* Rq(:, 1:end-1);
    resid_j = max(abs(Q(:, 2:end) - pred_q), [], 'all');
    max_resid = max(max_resid, resid_j);
end
assert(max_resid < 1e-6, ...
       'T5 fail: flow-dynamics residual %.2e on a test trajectory', max_resid);
fprintf('T5 ok: q_{k+1} = a q_k + (1 - a) r_q,k holds (max residual %.2e kg/s)\n', max_resid);

fprintf('\nAll predictor-validation invariants verified.\n');
