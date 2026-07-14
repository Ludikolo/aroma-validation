% CHECK_EXCITATION  Data-excitation check on the saved PRBS training set.
%
%   Reads the saved PRBS training trajectories (the op_* days are smooth
%   operational data, not excitation; check_conditioning covers them) and
%   verifies, without re-simulating anything, that the excitation channels
%   are rich enough for the V_seq fit:
%
%   E1  Per-channel persistency of excitation: the mosaic Hankel matrix
%       (one Hankel block per trajectory, concatenated) of T_0s and of
%       every per-edge r_q command has full rank at depth 60, well past
%       the fit horizon H_max. Trajectory means are removed first so the
%       rank measures the excitation, not the operating point.
%   E2  Joint PE at the fit horizon: the block Hankel of [T_0s; r_q] at
%       depth H_max has full row rank, which is the excitation the
%       multi-step regressor actually needs.
%   E3  Range coverage: every edge's r_q command reaches both PRBS bounds
%       [0.3, 1.5] x nominal edge flow, and T_0s spans its full 2 x amp
%       excitation band.
%   E4  Demand-phase coverage: the trajectory day-offsets are all distinct
%       and leave no gap larger than 24 h / n_train in the diurnal demand
%       cycle, so every demand regime (morning ramp, midday, evening peak,
%       night) is in the training set.

clear; clc;
startup;
p = params();

results_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results');
[net, ~] = build_plant(p);

tr = cell(p.data.n_train, 1);
for j = 1:p.data.n_train
    S = load(fullfile(results_dir, sprintf('train_%02d.mat', j)));
    tr{j} = S.traj;
end
n_edges = size(tr{1}.r_q, 1);

fprintf('\n=== Excitation check: %d PRBS trajectories, T_0s + %d r_q edges ===\n', ...
        numel(tr), n_edges);

%% E1 per-channel PE to depth 60 (mosaic Hankel over the trajectories)
L = 60;
rk_T0s = mosaic_hankel_rank(tr, @(t) t.T_0s(:)', L);
assert(rk_T0s == L, 'E1 fail: T_0s Hankel rank %d < %d at depth %d', rk_T0s, L, L);
rk_rq = zeros(n_edges, 1);
for e = 1:n_edges
    rk_rq(e) = mosaic_hankel_rank(tr, @(t) t.r_q(e, :), L);
end
assert(all(rk_rq == L), 'E1 fail: r_q edge %d Hankel rank %d < %d at depth %d', ...
       find(rk_rq < L, 1), min(rk_rq), L, L);
fprintf('E1 ok: full Hankel rank %d at depth %d for T_0s and all %d r_q edges\n', ...
        L, L, n_edges);

%% E2 joint PE at the fit horizon
Lh   = p.o1.H_max;
n_ch = 1 + n_edges;
blocks = cell(numel(tr), 1);
for j = 1:numel(tr)
    U  = [tr{j}.T_0s(:)'; tr{j}.r_q];
    U  = U - mean(U, 2);
    Nc = size(U, 2) - Lh + 1;
    B  = zeros(n_ch * Lh, Nc);
    for i = 1:Lh
        B((i-1)*n_ch + (1:n_ch), :) = U(:, i:i+Nc-1);
    end
    blocks{j} = B;
end
rk_joint = rank(horzcat(blocks{:}));
assert(rk_joint == n_ch * Lh, ...
       'E2 fail: joint block Hankel rank %d < %d at depth %d', rk_joint, n_ch * Lh, Lh);
fprintf('E2 ok: joint [T_0s; r_q] block Hankel full rank %d at depth H_max = %d\n', ...
        rk_joint, Lh);

%% E3 range coverage
lo   = p.excite.r_q_lo_factor;
hi   = p.excite.r_q_hi_factor;
mdot = net.mdotEdges(:);
Rq_all  = [];
T0s_all = [];
for j = 1:numel(tr)
    Rq_all  = [Rq_all  tr{j}.r_q];      %#ok<AGROW>
    T0s_all = [T0s_all tr{j}.T_0s(:)']; %#ok<AGROW>
end
rq_lo_dev = max(abs(min(Rq_all, [], 2) - lo * mdot) ./ mdot);
rq_hi_dev = max(abs(max(Rq_all, [], 2) - hi * mdot) ./ mdot);
assert(rq_lo_dev < 1e-9 && rq_hi_dev < 1e-9, ...
       'E3 fail: r_q misses the PRBS bounds (worst rel dev %.2e lo / %.2e hi)', ...
       rq_lo_dev, rq_hi_dev);
span = max(T0s_all) - min(T0s_all);
assert(span >= 2 * p.excite.T0s_amp - 1e-9, ...
       'E3 fail: T_0s span %.3f K < %.1f K excitation band', span, 2 * p.excite.T0s_amp);
fprintf('E3 ok: every edge hits [%.1f, %.1f] x nominal flow; T_0s spans %.1f K (band %.1f K)\n', ...
        lo, hi, span, 2 * p.excite.T0s_amp);

%% E4 demand-phase coverage (staggered day offsets)
T_day = 24 * 3600;
offs  = sort(cellfun(@(t) t.t_offset, tr));
assert(numel(uniquetol(offs, 1, 'DataScale', 1)) == numel(tr), ...
       'E4 fail: duplicate demand-phase offsets in the training set');
gaps = diff([offs; offs(1) + T_day]);   % circular gaps over the 24 h cycle
assert(max(gaps) <= T_day / numel(tr) + 1, ...
       'E4 fail: largest demand-phase gap %.2f h > %.2f h', ...
       max(gaps) / 3600, T_day / numel(tr) / 3600);
fprintf('E4 ok: %d distinct day-offsets cover the demand cycle, max gap %.1f h\n', ...
        numel(tr), max(gaps) / 3600);

fprintf('\nAll excitation checks passed.\n');


%% ---- local functions ----
function rk = mosaic_hankel_rank(tr, chan, L)
% Rank of the depth-L Hankel blocks of one channel, concatenated over the
% trajectories (mean removed per trajectory).
blocks = cell(numel(tr), 1);
for j = 1:numel(tr)
    u = chan(tr{j});
    u = u - mean(u);
    N = numel(u);
    H = zeros(L, N - L + 1);
    for i = 1:L
        H(i, :) = u(i:N-L+i);
    end
    blocks{j} = H;
end
rk = rank(horzcat(blocks{:}));
end
