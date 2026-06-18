% VIZ_COMPARE  KPC vs three baselines on the severe stressed day.
%
% Renders the head-to-head saved by run_comparison.m (run that first; it executes the
% four controllers on the same plant, scenario, warmup and actuator limits, so only
% the controller differs). Two time-domain panels: comfort (delivered heat vs demand
% for the binding consumer) and speed (solve time per step against the Ts = 900 s
% budget). Under the severe bind KPC tracks the demand best of all four and stays
% above the comfort bar, while the iterated Koopman-LMPC and the Jacobian LMPC fall
% below it; the exact-model NMPC is close on comfort but re-integrates the nonlinear
% ODE every solve, so it is an order of magnitude slower and sometimes fails to
% converge in the budget here.

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here))); startup();

S = load(fullfile(here, 'comparison.mat'));
p = S.p; names = S.names;
R  = {S.res_kpc, S.res_klmpc, S.res_jlmpc, S.res_nmpc};
cc = {[0.10 0.30 0.80], [0.20 0.60 0.20], [0.85 0.50 0.10], [0.60 0.10 0.60]};  % KPC, K-LMPC, J-LMPC, NMPC
N  = size(S.res_kpc.c_i, 2); th = (0:N-1) * p.Ts / 3600;

% binding consumer = the one with the lowest demand-met across the controllers
pcmet  = @(r) 100 * sum(r.c_i, 2) ./ max(sum(r.d_i, 2), 1e-9);
allmet = cell2mat(cellfun(@(r) pcmet(r), R, 'UniformOutput', false));   % n_user x 4
[~, bc] = min(min(allmet, [], 2));

figure('Color', 'w', 'Position', [100 100 1000 700]);

% --- comfort: binding-consumer delivered heat vs demand ---
subplot(2, 1, 1);
plot(th, S.res_kpc.d_i(bc, :)/1e3, 'k:', 'LineWidth', 2.4); hold on;     % demand (dotted)
for j = 1:4
    plot(th, R{j}.c_i(bc, :)/1e3, 'Color', cc{j}, 'LineWidth', 1.3 + 0.9*(j==1));   % KPC thicker
end
grid on; xlabel('time [h]'); ylabel('delivered heat [kW]'); xlim([0 24]);
legend(['demand', names], 'Location', 'best');
title(sprintf('Delivered vs demand, binding consumer C%d (stressed day, mdot = %.2f)', bc, S.mdot_scale));

% --- speed: solve time per step, against the Ts budget ---
subplot(2, 1, 2);
for j = 1:4
    semilogy(th, max(R{j}.solve_ms, 1e-1), 'Color', cc{j}, 'LineWidth', 1.2 + 0.8*(j==1)); hold on;
end
yline(p.Ts * 1000, 'r--', 'Ts = 900 s budget', 'LineWidth', 1.0);
grid on; xlabel('time [h]'); ylabel('solve time [ms]'); xlim([0 24]);
legend(names, 'Location', 'best');
title('Solve time per step: Koopman QP vs plant-rolling baselines');

exportgraphics(gcf, fullfile(here, 'viz_compare.pdf'), 'ContentType', 'vector');

% --- headline ---
fprintf('\nComparison (stressed day, mdot = %.2f, Np = %d):\n', S.mdot_scale, S.Np);
fprintf('%-15s  worst-met%%  unmet[kWh]  med-solve\n', 'controller');
for j = 1:4
    m = S.M.(matlab.lang.makeValidName(names{j}));
    fprintf('%-15s  %8.3f  %9.2f  %8.0f ms\n', names{j}, m.worst, m.unmet, m.med);
end
fprintf('Saved viz_compare.pdf\n');
