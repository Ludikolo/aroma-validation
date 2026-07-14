% PLOT_FLOW_BALANCE  Kirchhoff balance of the applied flow commands.
%
% The QP carries the hard equality [M_supply; M_return] r_q = 0 on every
% horizon step (lmpc_solve.m / kpc_v2_solve.m). This script re-applies those
% same incidence matrices to the flow commands the controller actually sent
% to the plant on the shipped stressed day (run demo_controller first) and
% plots the worst junction residual per step. A flat trace at solver
% tolerance, orders below the negligibility line, shows the constraint holds
% on what was applied, not only inside the optimizer.

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here)); startup();

S = load(fullfile(here, 'results', 'demo.mat'));
p = params(); [net, ~] = build_plant(p);
K = build_incidence_v2(net);
M_block = [K.M_supply; K.M_return];

rq   = S.res_kpc.r_q;                          % applied commands, n_edges x N
res  = max(abs(M_block * rq), [], 1);          % worst junction per step
th   = (0:size(rq, 2)-1) * p.Ts / 3600;

figure('Color', 'w', 'Position', [100 100 900 380]);
semilogy(th, max(res, 1e-16), '-', 'LineWidth', 1.6, 'Color', [0.10 0.30 0.80]);
hold on;
yline(1e-6, 'r--', 'negligible vs \sim0.1 kg/s edge flows', 'LineWidth', 1.0);
grid on; xlim([0 24]);
xlabel('time [h]'); ylabel('max junction residual [kg/s]');
title('Kirchhoff balance of the applied flow commands (stressed day)');
exportgraphics(gcf, fullfile(here, 'figures', 'flow_balance.pdf'), 'ContentType', 'vector');

fprintf('worst junction residual over the day: %.3g kg/s (%d steps)\n', ...
        max(res), numel(res));
fprintf('Saved figures/flow_balance.pdf\n');
