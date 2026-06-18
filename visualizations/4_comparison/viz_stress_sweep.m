% VIZ_STRESS_SWEEP  Render the stress sweep saved by run_stress_sweep.
% Worst-consumer demand-met for KPC and the iterated Koopman-LMPC as the flow
% capacity is tightened: the two track together at easy capacity and fan apart
% as it gets harder, with KPC holding above the comfort bar longer.

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here))); startup();

S = load(fullfile(here, 'stress_sweep.mat'));

figure('Color', 'w', 'Position', [100 100 760 520]);
plot(S.mdots, S.met_kpc,  '-o', 'Color', [0.10 0.30 0.80], 'LineWidth', 2.2, 'MarkerFaceColor', [0.10 0.30 0.80]); hold on;
plot(S.mdots, S.met_lmpc, '-s', 'Color', [0.20 0.60 0.20], 'LineWidth', 2.0, 'MarkerFaceColor', [0.20 0.60 0.20]);
yline(95, 'k--', '95% comfort bar', 'LineWidth', 1.0, 'LabelHorizontalAlignment', 'left');
set(gca, 'XDir', 'reverse'); grid on;
xlabel('flow capacity (mdot scale)   [tighter to the right]');
ylabel('worst-consumer demand-met [%]');
legend({'KPC (direct multi-step)', 'Koopman-LMPC (iterated one-step)'}, 'Location', 'southwest');
title('Direct multi-step pulls ahead as the network is squeezed harder');
exportgraphics(gcf, fullfile(here, 'stress_sweep.pdf'), 'ContentType', 'vector');
fprintf('Saved stress_sweep.pdf\n');
