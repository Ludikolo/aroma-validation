% PLOT_DELAY_PARITY  Measured pulse arrival vs the geometry-computed delay.
%
% Parity view of the strongest quantitative claim of the pulse experiment
% (run demo_pulse first): for each consumer, the arrival time of the +5 K
% source pulse (measured) against the transport delay computed independently
% from the network itself (pipe volume over volumetric flow, summed along the
% supply path). demo_pulse asserts every ratio inside [0.7, 1.2]; this figure
% shows the same check at a glance. Points on the diagonal mean the heat
% front travels at the bulk velocity the geometry predicts.

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here)); startup();

S = load(fullfile(here, 'results', 'pulse.mat'));
pu = S.pulse;
tc = pu.tau_comp(:) / 3600;   % computed delay [h]
ta = pu.t_arr(:)    / 3600;   % measured arrival after pulse onset [h]

figure('Color', 'w', 'Position', [100 100 560 520]);
lim = [0, 1.15 * max([tc; ta])];
fill([lim, fliplr(lim)], [0.7*lim, fliplr(1.2*lim)], [0.92 0.95 1.00], ...
     'EdgeColor', 'none'); hold on;
plot(lim, lim, 'k-', 'LineWidth', 1.0);
plot(tc, ta, 'o', 'MarkerSize', 8, 'LineWidth', 1.6, ...
     'MarkerFaceColor', [0.10 0.30 0.80], 'Color', [0.10 0.30 0.80]);
for i = 1:numel(tc)
    text(tc(i), ta(i), sprintf('  %s', pu.consumers{i}), 'FontSize', 10);
end
grid on; axis equal; xlim(lim); ylim(lim);
xlabel('computed delay from geometry [h]'); ylabel('measured pulse arrival [h]');
legend({'accepted band [0.7, 1.2] \times computed', 'parity', 'consumers'}, ...
       'Location', 'northwest');
title('Transport delay: measurement vs geometry');
exportgraphics(gcf, fullfile(here, 'figures', 'delay_parity.pdf'), 'ContentType', 'vector');

fprintf('arrival / computed ratios: %s (assert band [0.7, 1.2])\n', ...
        mat2str(round(ta ./ tc, 3)'));
fprintf('Saved figures/delay_parity.pdf\n');
