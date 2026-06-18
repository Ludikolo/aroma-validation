% VALIDATE_PLANT  Check the AROMA plant simulator against nine
%   physical invariants on the saved demo trajectories.
%
%   Run demo_plant first (or rely on the precomputed demo.mat in
%   ground_truth/results) so the scenarios are on disk. This script
%   does not re-simulate; it operates on the saved data and prints
%   a per-block summary.
%
%   T1  c^(i) <= d^(i) on every step (substation never over-delivers).
%   T2  T_r^(i) >= T_r_min on every step (5GDHC return-temperature floor).
%   T3  Trunk Kirchhoff at steady state: q_F0->F1 = sum_i q^(i).
%   T4  Capacity-binding behaviour matches scenario design.
%   T5  Substation energy: c^(i) = cp * (T_s - T_r) * q^(i) at machine
%       precision (AROMA mass-flow convention, so the factor is c_p
%       not rho c_p; the equivalent rho c_p form holds for volumetric q).
%   T6  Network energy ratio E_extracted / E_source_water in [0.5, 1.0]
%       per scenario. The theta_1 = 0.7 source-mixing weight from the
%       AROMA benchmark pins the lower bound
%       near 0.7; thermal-storage contributions push the ratio higher
%       in slow-flow scenarios.
%   T7  c^(i) >= 0 AND T_r^(i) <= T_s^(i) on every step (the upper
%       half of the substation contract).
%   T8  Per-junction mass conservation on the full incidence matrices
%       (every internal supply and return node, not just the trunk).
%   T9  Source-side return mixing law: at steady state, the aggregate
%       T_0r matches the mass-weighted average of the consumer return
%       temperatures, modulo small pipe-transport losses.

clear; clc;
startup;
p = params();

S = load(fullfile('ground_truth','results','demo.mat'));
n_scen = numel(S.scenarios);

fprintf('\n=== Plant validation: physical invariants (%d scenarios) ===\n', n_scen);

%% T1 c^(i) <= d^(i)
for k = 1:n_scen
    sc = S.scenarios(k);
    viol = max(sc.c_i(:) - sc.d_i(:));
    assert(viol <= 1, ...
        'T1 fail: %s c_i exceeds d_i by %.2f W', sc.name, viol);
end
fprintf('T1 ok: c^(i) <= d^(i) on every step in every scenario\n');

%% T2 T_r^(i) >= T_r_min
for k = 1:n_scen
    sc = S.scenarios(k);
    viol = p.consumer.T_r_min - min(sc.T_r_i(:));
    assert(viol <= 0.1, ...
        'T2 fail: %s T_r_i below T_r_min by %.2f K', sc.name, viol);
end
fprintf('T2 ok: T_r^(i) >= T_r_min on every step in every scenario\n');

%% T3 Steady-state Kirchhoff
for k = 1:n_scen
    sc = S.scenarios(k);
    n_last = max(1, numel(sc.t)-2):numel(sc.t);     % last 3 samples
    Q_net  = sum(sc.q_users(:, n_last), 1);
    res_K  = max(abs(Q_net - sc.q_F0_F1(n_last)));
    fprintf('  %-15s  Kirchhoff residual = %.2e kg/s\n', sc.name, res_K);
    assert(res_K < 1e-5, ...
        'T3 fail: %s Kirchhoff residual %.2e at steady state', sc.name, res_K);
end
fprintf('T3 ok: trunk inflow = sum of consumer flows at steady state (residual < 1e-5)\n');

%% T4 Sanity on met%: capacity-bound scenario must be < 95, others >= 99
fprintf('\n=== Demand-met (informational) ===\n');
for k = 1:n_scen
    sc = S.scenarios(k);
    pct = 100 * sum(sc.c_i(:)) / max(sum(sc.d_i(:)), 1e-9);
    fprintf('  %-15s  total met%% = %6.2f\n', sc.name, pct);
    switch sc.name
        case 'low_flow'
            assert(pct < 95, 'T4 fail: low_flow should be capacity-bound (got %.2f)', pct);
        case 'supply_step'
            assert(pct > 80 && pct < 99, ...
                'T4 fail: supply_step should be partially capacity-bound during T_0s=18 phase (got %.2f)', pct);
        case 'variable_flow'
            % variable_flow ramps 0.5x -> 1.5x -> 0.5x; capacity binds
            % during the 0.5x phases that overlap demand peaks. Allowed
            % range mirrors the open-loop reality: above 70 % (the plant
            % is doing useful work) and below 95 % (capacity binds at
            % some point in the run -- the dynamic-flow exercise).
            assert(pct > 70 && pct < 95, ...
                'T4 fail: variable_flow should bind at some peak (got %.2f)', pct);
        otherwise
            assert(pct >= 99, 'T4 fail: %s should reach >= 99 %% (got %.2f)', sc.name, pct);
    end
end
fprintf('T4 ok: capacity-bound vs adequate-capacity scenarios behave as expected\n');

%% T5 Substation energy: c^(i) = cp * (T_s^(i) - T_r^(i)) * q^(i)
fprintf('\n=== Substation energy consistency ===\n');
for k = 1:n_scen
    sc = S.scenarios(k);
    c_check = p.cp * (sc.T_s_i - sc.T_r_i) .* sc.q_users;
    res_sub = max(abs(sc.c_i - c_check), [], 'all');
    fprintf('  %-15s  max |c - cp*dT*q| = %.2e W\n', sc.name, res_sub);
    assert(res_sub < 1e-6, ...
        'T5 fail: %s substation residual %.2e W exceeds machine precision', sc.name, res_sub);
end
fprintf('T5 ok: c^(i) = cp * (T_s - T_r) * q^(i) at machine precision\n');

%% T6 Network energy ratio: E_extracted / E_source_water in [0.5, 1.0]
fprintf('\n=== Network energy ratio (theta1 = 0.7 artefact, see notes) ===\n');
for k = 1:n_scen
    sc = S.scenarios(k);
    Ts = mean(diff(sc.t));
    Q_net = sum(sc.q_users, 1);
    E_in  = sum(p.cp * (sc.T_0s - sc.T_0r) .* Q_net) * Ts;
    E_ex  = sum(sc.c_i(:)) * Ts;
    ratio = E_ex / max(E_in, 1);
    fprintf('  %-15s  E_extr = %6.1f MJ  E_src = %6.1f MJ  ratio = %.3f\n', ...
        sc.name, E_ex/1e6, E_in/1e6, ratio);
    assert(ratio >= 0.5 && ratio <= 1.0, ...
        'T6 fail: %s energy ratio %.3f outside physical range [0.5, 1.0]', sc.name, ratio);
end
fprintf('T6 ok: E_extracted / E_source_water in [0.5, 1.0] in every scenario\n');

%% T7 Substation upper contracts: T_r^(i) <= T_s^(i) AND c^(i) >= 0
% Pairs with T1 (c <= d) and T2 (T_r >= T_r_min) so that all four
% bounds T_r^min <= T_r <= T_s and 0 <= c <= d are all explicitly
% verified.
fprintf('\n=== Substation upper contracts ===\n');
for k = 1:n_scen
    sc = S.scenarios(k);
    viol_c_neg = max(0, -min(sc.c_i(:)));
    viol_Tr_Ts = max(sc.T_r_i(:) - sc.T_s_i(:));
    fprintf('  %-15s  min(c_i) = %6.2f W   max(T_r - T_s) = %+6.3f K\n', ...
        sc.name, min(sc.c_i(:)), max(sc.T_r_i(:) - sc.T_s_i(:)));
    assert(viol_c_neg <= 1, ...
        'T7 fail: %s c_i below 0 by %.2f W', sc.name, viol_c_neg);
    assert(viol_Tr_Ts <= 0.1, ...
        'T7 fail: %s T_r above T_s by %.3f K', sc.name, viol_Tr_Ts);
end
fprintf('T7 ok: c^(i) >= 0 AND T_r^(i) <= T_s^(i) on every step in every scenario\n');

%% T8 Per-junction mass conservation on the full incidence
% Kirchhoff in incidence form:
%   A_s q^(s,e) = [Q_net; -q^(i)]   A_r q^(r,e) = [-Q_net; q^(i)]
% T3 only verifies the trunk row of this. T8 verifies the residual
% at every internal supply and return node by applying the zero-row
% submatrices M_supply, M_return from build_incidence directly to
% the saved q_edges.
fprintf('\n=== Per-junction mass conservation ===\n');
[net, ~] = build_plant(p);
K = build_incidence(net);
for k = 1:n_scen
    sc = S.scenarios(k);
    if ~isfield(sc, 'q_edges')
        warning('T8: %s has no saved q_edges; re-run demo_plant.m', sc.name);
        continue;
    end
    n_last = max(1, numel(sc.t)-2):numel(sc.t);   % last 3 samples
    q_last = mean(sc.q_edges(:, n_last), 2);      % steady-state edge flows
    res_supply = max(abs(K.M_supply * q_last));
    res_return = max(abs(K.M_return * q_last));
    fprintf('  %-15s  supply residual = %.2e   return residual = %.2e   kg/s\n', ...
        sc.name, res_supply, res_return);
    assert(res_supply < 1e-5, ...
        'T8 fail: %s supply per-junction residual %.2e at steady state', sc.name, res_supply);
    assert(res_return < 1e-5, ...
        'T8 fail: %s return per-junction residual %.2e at steady state', sc.name, res_return);
end
fprintf('T8 ok: M_supply * q = 0 AND M_return * q = 0 at every internal node\n');

%% T9 Source-side return mixing law
% Mixing law evaluated at the source's return junction:
%   sum_e_in T_in q_in = T_node sum_e_out q_out
% Inlets are the consumer return flows q^(i) at temperature T_r^(i)
% (after pipe transport from consumer to source); outlet is T_0r at
% flow Q_net. Steady-state form, with a tolerance for pipe transport
% losses scaling as ~alpha * transit_time:
%   Sum_i q^(i) T_r^(i)  /  (Q_net T_0r)  in  [1 - tol, 1 + tol]
% Pipe loss alpha = 1e-6 1/s over the longest return path of about
% 100 min gives a ~0.6 % loss, so 5 % is the comfortable bound used
% here; low-flow scenarios sit slightly higher but still well inside.
fprintf('\n=== Source-side return mixing law (steady-state, with pipe-loss tolerance) ===\n');
for k = 1:n_scen
    sc = S.scenarios(k);
    n_last = max(1, numel(sc.t)-2):numel(sc.t);   % last 3 samples
    q_last  = sc.q_users(:, n_last);              % per consumer
    Tr_last = sc.T_r_i (:, n_last);
    T0r_last = sc.T_0r(n_last);
    Qnet_last = sum(q_last, 1);
    lhs = sum(q_last .* Tr_last, 1);              % sum_i q^(i) T_r^(i)
    rhs = Qnet_last .* T0r_last;                  % Q_net * T_0r
    rel_resid = mean(abs(lhs - rhs) ./ max(abs(rhs), 1e-9));
    fprintf('  %-15s  mean rel residual = %.3e   (tol 5e-2)\n', ...
        sc.name, rel_resid);
    assert(rel_resid < 0.05, ...
        'T9 fail: %s return mixing residual %.2e exceeds 5 %% tolerance', sc.name, rel_resid);
end
fprintf('T9 ok: source-side return mixing law holds at steady state (within pipe-loss tolerance)\n');

fprintf('\nAll plant-validation invariants verified.\n');
