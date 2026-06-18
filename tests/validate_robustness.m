% VALIDATE_ROBUSTNESS  Part 4 acceptance checks on the head-to-head comparison and
% the 90-day longevity. Run run_comparison and run_longevity first; this loads
% their saved results and asserts six invariants. It does not re-simulate.
%
%   C1  KPC's QP is feasible on every comparison step (convex, reliable); each
%       baseline's convergence rate is reported.
%   C2  KPC stays above the 95 % comfort bar and is the best worst-consumer met %.
%   C3  KPC solves at least 10x faster than the nonlinear NMPC.
%   C4  KPC median solve time stays under 5 % of the Ts budget (real-time).
%   L1  90-day worst-consumer met % drift is small (<= 0.5 pp).
%   L2  90-day infeasible-step rate is negligible (<= 0.1 %).

clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
rob  = fullfile(root, 'robustness_scalability_sustainability', 'results');
cmp  = fullfile(root, 'visualizations', '4_comparison', 'comparison.mat');

C  = load(cmp);
LO = load(fullfile(rob, 'longevity.mat'));
fld = @(nm) C.M.(matlab.lang.makeValidName(nm));   % access M by controller name

fprintf('\n=== Part 4 validation: six invariants ===\n');

%% C1 KPC's convex QP is feasible on every step (reliability); report baseline rates
kpc = fld('KPC');
assert(kpc.feas == kpc.n, 'C1 fail: KPC feasible only %d / %d steps', kpc.feas, kpc.n);
fprintf('C1 ok: KPC feasible on every one of the %d steps; baseline convergence:\n', kpc.n);
for j = 1:numel(C.names)
    m = fld(C.names{j});
    fprintf('       %-15s %d / %d steps converged\n', C.names{j}, m.feas, m.n);
end

%% C2 KPC stays above the 95 % comfort bar and is the best of all controllers
others = cellfun(@(nm) fld(nm).worst, setdiff(C.names, {'KPC'}));
assert(kpc.worst >= 95, 'C2 fail: KPC worst-consumer met = %.3f %% < 95 %% bar', kpc.worst);
assert(kpc.worst >= max(others), 'C2 fail: KPC (%.3f %%) is not the best worst-consumer met %%', kpc.worst);
fprintf('C2 ok: KPC worst-consumer met = %.3f %% (>= 95 %% bar and best; next best %.3f %%)\n', ...
        kpc.worst, max(others));

%% C3 KPC at least 10x faster than NMPC
nmpc = fld('NMPC');
speedup = nmpc.med / max(kpc.med, eps);
assert(speedup >= 10, 'C3 fail: KPC speedup vs NMPC = %.1fx < 10x', speedup);
fprintf('C3 ok: KPC median solve %.0f ms vs NMPC %.0f ms (%.0fx faster)\n', kpc.med, nmpc.med, speedup);

%% C4 KPC real-time: median solve < 5 % of the Ts budget
p = params(); budget_ms = p.Ts * 1000;
assert(kpc.med < 0.05 * budget_ms, 'C4 fail: KPC median solve %.0f ms >= 5 %% of %g s budget', kpc.med, p.Ts);
fprintf('C4 ok: KPC median solve %.0f ms = %.3f %% of the %g s budget (real-time)\n', ...
        kpc.med, 100*kpc.med/budget_ms, p.Ts);

%% L1 longevity drift <= 0.5 pp
assert(abs(LO.drift_total) <= 0.5, 'L1 fail: 90-day drift = %+.4f pp (|.| > 0.5)', LO.drift_total);
fprintf('L1 ok: %d-day worst-consumer drift = %+.4f pp (small)\n', LO.n_days_actual, LO.drift_total);

%% L2 longevity infeasible-step rate <= 0.1 %
assert(LO.infeas_rate <= 0.1, 'L2 fail: infeasible-step rate = %.4f %% > 0.1 %%', LO.infeas_rate);
fprintf('L2 ok: %d-day infeasible-step rate = %.4f %% (negligible)\n', LO.n_days_actual, LO.infeas_rate);

fprintf('\n=== Part 4 validation: all 6 invariants passed ===\n');
