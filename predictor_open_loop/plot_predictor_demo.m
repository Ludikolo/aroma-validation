function plot_predictor_demo(demo, p)
% PLOT_PREDICTOR_DEMO  The predictor demo figure:
%   rmse_horizon.pdf NRMSE vs horizon for ZOH, iterated A,B, V_seq.
% The per-consumer h = 1 parity data stays available in demo.mat
% (demo.parity_h1); the deployment-regime parity figure lives in
% visualizations/2_predictor.

outdir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(outdir, 'dir'), mkdir(outdir); end

%% NRMSE vs horizon, three predictors
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
ylabel('NRMSE  [% of d_{max}]', 'FontSize', 13);
title('Predictor accuracy with horizon', 'FontSize', 12);
legend('Location', 'northwest', 'FontSize', 11);
xticks(demo.horizons(1:2:end));
exportgraphics(fig, fullfile(outdir, 'rmse_horizon.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(outdir, 'rmse_horizon.png'), 'Resolution', 200);
close(fig);

fprintf('Saved rmse_horizon (.pdf + .png) to %s\n', outdir);

end
