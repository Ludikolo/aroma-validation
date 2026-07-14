% PLOT_DICTIONARY_GAP  Delivered-heat RMSE gap: the spatial lift minus the
%   47-feature production lift, per regime and horizon. Below zero the spatial
%   lift predicts better; above zero the 47 does.
%
% Values come from results/dict_gap.mat, which eval_dictionary_gap.m writes by
% scoring both lifts on the same trajectories. Run that first if the file is
% missing. Nothing here is typed in by hand.

clear; clc;
here = fileparts(mfilename('fullpath'));
figs = fullfile(here, 'figures');
if ~exist(figs, 'dir'), mkdir(figs); end

G = load(fullfile(here, 'results', 'dict_gap.mat'));
gap = G.gap;  H = G.H;  names = G.names;

% Group the five trajectories into the three regimes the study reports. Each
% regime holds more than one trajectory, so the bar is their mean and the
% whisker is their spread. That spread is the honest error bar here: it is the
% range across trajectories, not a bootstrap interval.
grp = {'PRBS (held out)',        find(contains(names, 'PRBS')); ...
       'operational (held out)', find(contains(names, 'op-day')); ...
       'on-policy (stressed)',   find(contains(names, 'on-policy'))};

show_h = [1 8 12];                     % one step, and the two deployed horizons
nG = size(grp, 1);  nH = numel(show_h);

val = zeros(nG, nH);  lo = zeros(nG, nH);  hi = zeros(nG, nH);
for g = 1:nG
    rows = grp{g, 2};
    for j = 1:nH
        c = find(H == show_h(j));
        v = gap(rows, c);
        val(g, j) = mean(v);  lo(g, j) = min(v);  hi(g, j) = max(v);
    end
end

fig = figure('Visible', 'off', 'Position', [100 100 780 400]);
b = bar(val, 'grouped');
cols = [0.35 0.55 0.75; 0.55 0.70 0.45; 0.78 0.32 0.22];
for j = 1:nH, b(j).FaceColor = cols(j, :); end
hold on;
gw = min(0.8, nH / (nH + 1.5));
for j = 1:nH
    x = (1:nG) - gw/2 + (2*j - 1) * gw / (2*nH);
    errorbar(x, val(:, j), val(:, j) - lo(:, j), hi(:, j) - val(:, j), ...
             'k', 'LineStyle', 'none', 'LineWidth', 1.0);
end
yline(0, 'k-', 'LineWidth', 1.0);
set(gca, 'XTickLabel', grp(:, 1), 'FontSize', 10);
ylabel('RMSE(spatial) - RMSE(47)   [W]');
legend(arrayfun(@(h) sprintf('h = %d', h), show_h, 'UniformOutput', false), ...
       'Location', 'northwest');
title('The spatial lift wins only one step ahead on PRBS; the 47 wins on-policy');

exportgraphics(fig, fullfile(figs, 'dict_gap.pdf'), 'ContentType', 'vector');
close(fig);
fprintf('Saved figures/dict_gap.pdf\n');
for g = 1:nG
    fprintf('  %-24s', grp{g, 1});
    for j = 1:nH
        fprintf('   h=%2d: %+4.0f [%+4.0f, %+4.0f]', show_h(j), val(g, j), lo(g, j), hi(g, j));
    end
    fprintf('\n');
end
