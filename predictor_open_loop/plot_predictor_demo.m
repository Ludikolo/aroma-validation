function plot_predictor_demo(demo, p)
% PLOT_PREDICTOR_DEMO  Two figures for the predictor demo:
%   parity.pdf       h = 1 parity scatter, per consumer, with R^2 in the legend.
%   rmse_horizon.pdf NRMSE vs horizon for ZOH, iterated A,B, V_seq.

outdir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(outdir, 'dir'), mkdir(outdir); end

cmap   = lines(5);
labels = {'C1', 'C2', 'C3', 'C4', 'C5'};

%% Figure 1: parity at h = 1 (V_seq prediction vs true c on test trajectories)
fig = figure('Visible', 'off', 'Position', [100 100 720 700]);
hold on; grid on; box on;
yt_all = []; yp_all = [];
for i = 1:5
    yt = demo.parity_h1.y_true{i} / 1e3;
    yp = demo.parity_h1.y_pred{i} / 1e3;
    yt_all = [yt_all, yt];
    yp_all = [yp_all, yp];
    scatter(yt, yp, 26, cmap(i, :), 'filled', 'MarkerFaceAlpha', 0.55, ...
        'DisplayName', sprintf('%s (R^2 = %.3f)', labels{i}, demo.R2_vseq_h1(i)));
end
lo = min([yt_all, yp_all]);
hi = max([yt_all, yp_all]);
pad = 0.05 * (hi - lo);
plot([lo - pad, hi + pad], [lo - pad, hi + pad], 'k-', ...
     'LineWidth', 1.0, 'HandleVisibility', 'off');
xlabel('true heat extraction c^{(i)} [kW]', 'FontSize', 13);
ylabel('predicted c^{(i)} [kW]', 'FontSize', 13);
title(sprintf('V_{seq} 1-step prediction on test trajectories, suite R^2 = %.3f', ...
              mean(demo.R2_vseq_h1)), 'FontSize', 12);
legend('Location', 'northwest', 'FontSize', 11);
axis equal;
xlim([lo - pad, hi + pad]);
ylim([lo - pad, hi + pad]);
exportgraphics(fig, fullfile(outdir, 'parity.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(outdir, 'parity.png'), 'Resolution', 200);
close(fig);

%% Figure 2: NRMSE vs horizon, three predictors
fig = figure('Visible', 'off', 'Position', [100 100 880 480]);
hold on; grid on; box on;
plot(demo.horizons, demo.nrmse_zoh,  'o-', 'Color', [0.50 0.50 0.50], ...
    'MarkerFaceColor', [0.50 0.50 0.50], 'MarkerSize', 8, ...
    'DisplayName', 'ZOH (no model)');
plot(demo.horizons, demo.nrmse_iter, 's-', 'Color', [0.85 0.40 0.10], ...
    'MarkerFaceColor', [0.85 0.40 0.10], 'MarkerSize', 8, ...
    'DisplayName', 'iterated A,B');
plot(demo.horizons, demo.nrmse_vseq, 'd-', 'Color', [0.10 0.40 0.85], ...
    'MarkerFaceColor', [0.10 0.40 0.85], 'MarkerSize', 10, 'LineWidth', 2.2, ...
    'DisplayName', 'V_{seq} direct multi-step (this work)');
xlabel(sprintf('horizon h [steps of %g s]', p.Ts), 'FontSize', 13);
ylabel('NRMSE / d_{max}  [%]', 'FontSize', 13);
title('Predictor accuracy with horizon', 'FontSize', 12);
legend('Location', 'northwest', 'FontSize', 11);
xticks(demo.horizons(1:2:end));
exportgraphics(fig, fullfile(outdir, 'rmse_horizon.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(outdir, 'rmse_horizon.png'), 'Resolution', 200);
close(fig);

fprintf('Saved parity + rmse_horizon (.pdf + .png) to %s\n', outdir);

end
