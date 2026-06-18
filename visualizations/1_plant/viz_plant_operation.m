% VIZ_PLANT_OPERATION  A normal day, plus two conservation checks.
%
% I run the plant for a full 24 h under the diurnal demand at nominal flow and plot
% the main signals (source and consumer supply/return temperatures, the consumer
% flows, total delivered vs demanded heat). Then I check two conservation laws on the
% trajectory: substation energy balance c = cp (T_s - T_r) q at every step, and
% Kirchhoff mass conservation M q = 0 at every junction. Both residuals come out at
% machine precision, and the signals look sensible (return below supply, return above
% the 15 C floor, delivery tracking demand).
%
% Runs live against build_plant + simulate_plant; nothing is loaded from a .mat.

clear; clc;

% --- put the aroma-build code on the path (same setup the demos use) ---
here = fileparts(mfilename('fullpath'));   % .../visualizations/1_plant
root = fileparts(fileparts(here));         % .../aroma-build (repo root)
addpath(root); startup();

p = params();
p.Ts = 60;                 % read-out cadence [s]

[net, z0_cold] = build_plant(p);
ei     = edge_user_index(net);
n_user = ei.n_user;
R0_idx = find(strcmp({net.Nodes.name}, 'R0'));   % source return node, for T_0r
mdotE  = net.mdotEdges(:);                        % nominal per-edge flow [kg/s]

t_start = 12.5 * 3600;   % start 12:30 [s]
T_warm  = 90  * 60;      % warmup [s]
T_sim   = 24  * 3600;    % one full day [s]
r_q_fun = @(t) mdotE;    % nominal flow throughout

% warmup, then a 24 h run at nominal supply (T_0s = Tin_nom, so offset u = 0)
p_wu = p; p_wu.t_offset = t_start; p_wu.r_q_fun = r_q_fun;
res_wu = simulate_plant(net, z0_cold, p_wu, @(t) 0, @(t) 1.0, T_warm);
net_run = net; net_run.q0 = res_wu.q_edges(:, end);
p_run = p; p_run.t_offset = t_start + T_warm; p_run.r_q_fun = r_q_fun;
res = simulate_plant(net_run, res_wu.z_final, p_run, @(t) 0, @(t) 1.0, T_sim);

t_h  = res.t / 3600;
T_0s = p.Tin_nom + res.u;        % source supply temp [C]
T_0r = res.Tout(R0_idx, :);      % source return temp [C]

% --- figure: the main signals over the day ---
figure('Color', 'w');
subplot(2, 2, 1);
plot(t_h, T_0s, 'k', 'LineWidth', 1.3); hold on; plot(t_h, res.T_s_i', '-');
grid on; xlabel('time [h]'); ylabel('[\circC]'); title('supply temperatures (source + consumers)');
subplot(2, 2, 2);
plot(t_h, T_0r, 'k', 'LineWidth', 1.3); hold on; plot(t_h, res.T_r_i', '-');
grid on; xlabel('time [h]'); ylabel('[\circC]'); title('return temperatures (source + consumers)');
subplot(2, 2, 3);
plot(t_h, res.q_users', '-');
grid on; xlabel('time [h]'); ylabel('[kg/s]'); title('consumer flows q^{(i)}');
subplot(2, 2, 4);
plot(t_h, sum(res.d_i, 1)/1e3, 'k--', 'LineWidth', 1.6); hold on; plot(t_h, sum(res.c_i, 1)/1e3, 'b', 'LineWidth', 1.1);
grid on; xlabel('time [h]'); ylabel('total heat [kW]'); legend('demand', 'delivered', 'Location', 'best');
title('demand vs delivered');

% WHAT YOU SEE: temperatures and flows follow the diurnal demand, every return
% temperature stays below its supply and above the 15 C floor, and at nominal
% flow the delivered heat tracks the demand.

exportgraphics(gcf, fullfile(here, 'viz_plant_operation.pdf'), 'ContentType', 'vector');

% --- conservation checks (printed; both should be about machine precision) ---
% (1) substation energy balance: c = cp (T_s - T_r) q   (q in kg/s, AROMA convention)
c_check    = p.cp * (res.T_s_i - res.T_r_i) .* res.q_users;
res_energy = max(abs(res.c_i - c_check), [], 'all');
% (2) Kirchhoff mass conservation at every junction (supply and return networks)
K = build_incidence_v2(net);
res_supply = max(abs(K.M_supply * res.q_edges), [], 'all');
res_return = max(abs(K.M_return * res.q_edges), [], 'all');
fprintf('\nConservation checks (should be ~0):\n');
fprintf('  energy balance   max|c - cp*dT*q|  = %.2e W\n', res_energy);
fprintf('  mass, supply     max|M_supply * q| = %.2e kg/s\n', res_supply);
fprintf('  mass, return     max|M_return * q| = %.2e kg/s\n', res_return);
fprintf('Saved viz_plant_operation.pdf\n');
