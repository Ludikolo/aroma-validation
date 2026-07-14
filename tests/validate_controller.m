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
%   T2  KPC delivers 100 % met % on every (consumer, scenario) cell
%       of the 6-scenario suite at design capacity.
%   T3  KPC suite-mean met % at stressed capacity >= 95 % (above bar).
%   T4  KPC beats hold-nominal on stressed worst-consumer met %.
%   T5  Delivery contract c <= d on every (consumer, step) of the
%       realised trace. The plant enforces this by clipping at the
%       substation, so this checks the closed-loop bookkeeping.
%   T6  Every applied flow reference stays within the QP's per-edge box
%       and rate limit dr_q_max.
%   T7  Per-step solve time stays well under the sample budget
%       (mean <= 5 % of Ts).

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

%% T5 delivery contract c <= d (plant-side bookkeeping check)
%   The plant enforces c <= d by clipping at the substation, so this
%   cannot catch a controller fault; it confirms the closed-loop logs
%   respect the delivery contract on the realised trajectory.
gap = max(S.res_kpc.c_i - S.res_kpc.d_i, [], 'all');
assert(gap <= 1e-6, 'T5 fail: max(c - d) = %.3g W on the realised trace', gap);
fprintf('T5 ok: max(c - d) = %.3g W on the realised trace (plant clip holds)\n', gap);

%% T6 applied flow reference within the QP rate limit
%   res_kpc.r_q logs the flow reference the controller applied each
%   step (simulate_plant stores r_q_fun, not the lagged edge flow), so
%   consecutive columns must respect the hard rate bound dr_q_max the
%   QP puts on every edge.
%   The rate cap alone rarely binds at stressed capacity, so the box
%   bounds (which DO bind there) are asserted as well: every applied
%   reference must lie inside the QP's hard per-edge box, the same
%   factors kpc_v2_solve puts on r_q.
dr_step = max(abs(diff(S.res_kpc.r_q, 1, 2)), [], 'all');
assert(dr_step <= S.tune.dr_q_max + 1e-6, ...
    'T6 fail: applied flow step %.4f kg/s exceeds rate limit %.3f kg/s', ...
    dr_step, S.tune.dr_q_max);
net_nom = build_plant(p);
md_s   = S.mdot_scale * net_nom.mdotEdges(:);
lo_f   = p.excite.r_q_lo_factor;            % kpc_v2_solve lets the tune override the floor
if isfield(S.tune, 'r_q_lo_factor') && ~isempty(S.tune.r_q_lo_factor)
    lo_f = S.tune.r_q_lo_factor;
end
rq_lo  = lo_f * md_s;
rq_hi  = p.excite.r_q_hi_factor * md_s;
box_lo = min(S.res_kpc.r_q - rq_lo, [], 'all');
box_hi = min(rq_hi - S.res_kpc.r_q, [], 'all');
assert(box_lo >= -1e-6 && box_hi >= -1e-6, ...
    'T6 fail: applied r_q leaves the QP box (lo margin %.2e, hi margin %.2e)', box_lo, box_hi);
fprintf(['T6 ok: applied r_q inside the QP box on every edge and step, and max\n' ...
         '       flow step %.4f kg/s <= rate limit %.2f kg/s per step\n'], ...
        dr_step, S.tune.dr_q_max);

%% T7 per-step solve time well under Ts budget
mean_solve_ms = mean(S.res_kpc.solve_ms);
max_solve_ms  = max(S.res_kpc.solve_ms);
budget_ms = 1000 * p.Ts;
assert(mean_solve_ms <= 0.05 * budget_ms, ...
    'T7 fail: mean solve %.0f ms > 5 %% of Ts budget %.0f ms', mean_solve_ms, budget_ms);
fprintf('T7 ok: mean solve %.0f ms, max %.0f ms (Ts budget %.0f ms)\n', ...
        mean_solve_ms, max_solve_ms, budget_ms);

fprintf('\n=== Controller validation: all 7 invariants passed ===\n');
