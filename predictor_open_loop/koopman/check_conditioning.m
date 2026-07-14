% CHECK_CONDITIONING  Rank and conditioning of the production lift on the
%   training set.
%
%   Builds the 47-feature lift on the same train_* + op_* trajectories the
%   V_seq fit uses (from valid_start onward, as the fit sees them) and
%   checks:
%
%   C1  n_z = 47 (guards the dictionary flags in params.m).
%   C2  rank(Z) = 41: exactly 6 rank deficiencies, all exact by
%       construction and named:
%         Q_net   = sum_i q^i                       (1)
%         theta^i = exergy_s^i - exergy_r^i         (5, one per consumer)
%       Each identity is verified at machine precision. The ridge in
%       fit_vseq absorbs these 6 exact
%       null directions, so the conditioning that matters for the fit is
%       sigma_1/sigma_41 over the numerical rank, printed below.
%
%   Informational (no assert): the same line for the 78-feature spatial
%   lift of the dictionary study; it wraps rows 1..47 unchanged, so it
%   inherits the same 6 dependencies.

clear; clc;
startup;
p = params();

results_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results');
[net, ~] = build_plant(p);

files = {};
for j = 1:p.data.n_train
    files{end+1} = fullfile(results_dir, sprintf('train_%02d.mat', j)); %#ok<SAGROW>
end
op = dir(fullfile(results_dir, 'op_*.mat'));
for j = 1:numel(op)
    files{end+1} = fullfile(results_dir, op(j).name); %#ok<SAGROW>
end

Z = []; Zs = []; idx = [];
for j = 1:numel(files)
    S = load(files{j});
    [Zj, ~, meta] = candidate_library(S.traj, p);
    Z = [Z Zj(:, meta.valid_start:end)]; %#ok<AGROW>
    if j == 1, idx = meta.idx; end
    [Zsj, ~, ms] = spatial_library(S.traj, p, net, 'full');
    Zs = [Zs Zsj(:, ms.valid_start:end)]; %#ok<AGROW>
end

fprintf('\n=== Lift conditioning: %d trajectories (%d train + %d op), %d samples ===\n', ...
        numel(files), p.data.n_train, numel(op), size(Z, 2));

%% C1 production lift size
n_z = size(Z, 1);
assert(n_z == 47, 'C1 fail: n_z = %d, expected the 47-feature production lift', n_z);
fprintf('C1 ok: production lift n_z = %d\n', n_z);

%% C2 rank 41 with the 6 deficiencies named
s   = svd(Z);
tol = max(size(Z)) * eps(max(s));
rnk = sum(s > tol);
assert(rnk == 41, 'C2 fail: rank(Z) = %d (deficiency %d), expected 41 (deficiency 6)', ...
       rnk, n_z - rnk);

n_dep  = 1 + numel(idx.theta);
res    = zeros(n_dep, 1);
labels = cell(n_dep, 1);
res(1)    = max(abs(Z(idx.Q_net, :) - sum(Z(idx.q, :), 1)));
labels{1} = 'Q_net = sum_i q_i';
for i = 1:numel(idx.theta)
    res(1+i)    = max(abs(Z(idx.theta(i), :) ...
                  - (Z(idx.exergy_s(i), :) - Z(idx.exergy_r(i), :))));
    labels{1+i} = sprintf('theta_C%d = exergy_s_C%d - exergy_r_C%d', i, i, i);
end
for m = 1:n_dep
    assert(res(m) < 1e-9, 'C2 fail: dependency "%s" has residual %.2e', labels{m}, res(m));
    fprintf('   dependency %d: %-38s (max residual %.1e)\n', m, labels{m}, res(m));
end
fprintf('C2 ok: rank(Z) = %d = 47 - 6, all %d deficiencies exact by construction\n', ...
        rnk, n_dep);
fprintf('   conditioning over the numerical rank: sigma_1/sigma_%d = %.2e (Gram %.2e)\n', ...
        rnk, s(1) / s(rnk), (s(1) / s(rnk))^2);

%% spatial 78-lift (informational, no assert)
ss   = svd(Zs);
tols = max(size(Zs)) * eps(max(ss));
rs   = sum(ss > tols);
fprintf('info: spatial lift n_z = %d, rank %d (deficiency %d), sigma_1/sigma_%d = %.2e\n', ...
        size(Zs, 1), rs, size(Zs, 1) - rs, rs, ss(1) / ss(rs));

fprintf('\nAll conditioning checks passed.\n');
