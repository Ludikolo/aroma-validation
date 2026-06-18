% VALIDATE_CONTROLLER  Check the KPC closed-loop controller
%   against seven invariants on the saved demo trajectories.
%
%   Run demo_controller first (or rely on the precomputed demo.mat in
%   controller_closed_loop/results) so the comparison is on disk. This
%   script does not refit anything; it operates on the saved data and
%   the committed tune artefact, prints a per-block summary, and exits
%   non-zero if any check fails.
%
%   T1  Every QP solve in the demo returned exitflag > 0 (feasible).
%   T2  KPC delivers 100 %% met %% on every (consumer, scenario) cell
%       of the 6-scenario suite at design capacity.
%   T3  KPC suite-mean met %% at stressed capacity >= 95 %% (above bar).
%   T4  KPC beats hold-nominal on stressed worst-consumer met %%.
%   T5  Substation contract c <= d on every (consumer, step) within
%       the QP's slack tolerance (no hard violation by the plant).
%   T6  Non-negative delivered heat c >= 0 on every (consumer, step).
%   T7  Per-step solve time stays well under the sample budget
%       (mean <= 5 %% of Ts).

clear; clc;
startup;
p = params();

results_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                       'controller_closed_loop', 'results');

S = load(fullfile(results_dir, 'demo.mat'));

fprintf('\n=== Controller validation: seven invariants ===\n');

%% T1 every QP solve feasible
ef = S.res_kpc.exitflag(:);
assert(all(ef > 0), ...
    'T1 fail: %d / %d steps had exitflag <= 0', sum(ef <= 0), numel(ef));
fprintf('T1 ok: every QP returned exitflag > 0 (%d steps)\n', numel(ef));

%% T2 100 %% met %% on every (consumer, scenario) cell at design
assert(all(S.met_design_kpc(:) >= 99.99), ...
    'T2 fail: KPC drops to %.3f %% on at least one design cell', ...
    min(S.met_design_kpc(:)));
fprintf('T2 ok: KPC met%% = 100 on all %d design cells\n', ...
        numel(S.met_design_kpc));

%% T3 stressed suite mean above the 95 %% bar
suite_mean_stress = mean(S.met_kpc);
assert(suite_mean_stress >= 95, ...
    'T3 fail: KPC suite-mean at stressed = %.3f %% < 95 %%', suite_mean_stress);
fprintf('T3 ok: KPC suite-mean at stressed = %.3f %% >= 95 %% bar\n', suite_mean_stress);

%% T4 KPC beats hold on stressed worst consumer
worst_kpc  = min(S.met_kpc);
worst_hold = min(S.met_hold);
assert(worst_kpc > worst_hold, ...
    'T4 fail: KPC worst (%.3f) not above hold worst (%.3f)', worst_kpc, worst_hold);
fprintf('T4 ok: KPC worst-consumer = %.3f %% > hold worst = %.3f %% (+%.2f pp)\n', ...
        worst_kpc, worst_hold, worst_kpc - worst_hold);

%% T5 substation contract c <= d (plant clipping must hold)
%   The plant clips c at the substation; the QP's S_hi slack only
%   describes the predicted gap. The validation asserts plant-side
%   contract on the realised trajectory.
gap = max(S.res_kpc.c_i - S.res_kpc.d_i, [], 'all');
assert(gap <= 1e-6, 'T5 fail: max(c - d) = %.3g W on the realised trace', gap);
fprintf('T5 ok: max(c - d) = %.3g W on the realised trace (plant clip holds)\n', gap);

%% T6 non-negative delivered heat
c_min = min(S.res_kpc.c_i, [], 'all');
assert(c_min >= -1e-6, 'T6 fail: min(c) = %.3g W < 0', c_min);
fprintf('T6 ok: min(c) = %.3g W >= 0 on the realised trace\n', c_min);

%% T7 per-step solve time well under Ts budget
mean_solve_ms = mean(S.res_kpc.solve_ms);
max_solve_ms  = max(S.res_kpc.solve_ms);
budget_ms = 1000 * p.Ts;
assert(mean_solve_ms <= 0.05 * budget_ms, ...
    'T7 fail: mean solve %.0f ms > 5 %% of Ts budget %.0f ms', mean_solve_ms, budget_ms);
fprintf('T7 ok: mean solve %.0f ms, max %.0f ms (Ts budget %.0f ms)\n', ...
        mean_solve_ms, max_solve_ms, budget_ms);

fprintf('\n=== Controller validation: all 7 invariants passed ===\n');
