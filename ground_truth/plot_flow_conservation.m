% PLOT_FLOW_CONSERVATION  Trunk flow vs the consumer sum through a flow ramp.
%
% One panel from the saved variable_flow scenario (run demo_plant first): the
% trunk flow q_{F0->F1} overlaid on the sum of the five consumer flows. At
% steady state the two coincide to solver precision (the T3/T8 invariants
% assert this at < 1e-9 kg/s); for a few samples after a flow move the
% independent per-edge first-order lags let them separate by a few percent,
% which is visible at the ramp corners and documented in demo_plant.m. The
% conservation claim is a steady-state claim, and this figure shows both the
% equality and its honest transient scope in one glance.

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here)); startup();

S  = load(fullfile(here, 'results', 'demo.mat'));
sc = S.scenarios(strcmp({S.scenarios.name}, 'variable_flow'));
th = (sc.t - sc.t(1)) / 3600;
qsum = sum(sc.q_users, 1);

figure('Color', 'w', 'Position', [100 100 900 420]);
plot(th, sc.q_F0_F1, '-',  'LineWidth', 2.0); hold on;
plot(th, qsum,       '--', 'LineWidth', 1.6);
grid on; xlim([0 24]);
xlabel('time [h]'); ylabel('flow [kg/s]');
legend({'trunk q_{F_0 \rightarrow F_1}', 'sum of consumer flows'}, 'Location', 'best');
title('Mass conservation: the trunk carries the consumer sum (variable-flow ramp)');
exportgraphics(gcf, fullfile(here, 'figures', 'flow_conservation.pdf'), 'ContentType', 'vector');

% console evidence: steady-state residual vs worst transient
res = abs(qsum - sc.q_F0_F1);
n_last = max(1, numel(sc.t)-4):numel(sc.t);
fprintf('steady-state residual (last samples): %.3g kg/s\n', max(res(n_last)));
fprintf('worst transient residual (ramp corners): %.3g kg/s (%.1f%% of trunk)\n', ...
        max(res), 100 * max(res ./ max(sc.q_F0_F1, 1e-9)));
fprintf('Saved figures/flow_conservation.pdf\n');
