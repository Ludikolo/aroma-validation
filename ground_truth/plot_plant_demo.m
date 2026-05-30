function plot_plant_demo(scenarios, p)
% PLOT_PLANT_DEMO  One 2x2-panel figure per scenario plus a summary
%   figure across all scenarios. Saved to ground_truth/figures/.
%
%   2x2 layout (sized so two scenario figures sit side by side on a
%   16:9 slide):
%     [1,1] Source (T_0s commanded, T_0r aggregate return)
%     [1,2] Per-consumer supply temperature T_s^(i)
%     [2,1] Per-consumer demand vs delivered (d^(i) solid, c^(i) dashed)
%     [2,2] Per-consumer flow q^(i) and the trunk q_{F0->F1}

outdir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
                  'ground_truth', 'figures');
if ~exist(outdir, 'dir'), mkdir(outdir); end
n_user = size(scenarios(1).T_s_i, 1);
cmap   = lines(n_user);

for k = 1:numel(scenarios)
    sc = scenarios(k);
    t_min = (sc.t - sc.t(1)) / 60;

    fig = figure('Visible', 'off', 'Position', [100 100 960 640]);
    tl  = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tl, sprintf('%s | %s', strrep(sc.name, '_', ' '), sc.desc), ...
          'FontWeight', 'bold');

    % [1,1] Source temperatures
    nexttile;
    plot(t_min, sc.T_0s, '-', 'LineWidth', 1.5); hold on;
    plot(t_min, sc.T_0r, '-', 'LineWidth', 1.5);
    yline(p.consumer.T_r_min, ':k', sprintf('T_r^{min} = %g', p.consumer.T_r_min));
    grid on;
    legend({'T_0^s commanded', 'T_0^r return'}, 'Location', 'best');
    ylabel('temperature [C]'); xlabel('time [min]');
    title('Source');

    % [1,2] Per-consumer supply temperature
    nexttile;
    hold on;
    for i = 1:n_user
        plot(t_min, sc.T_s_i(i, :), '-', 'Color', cmap(i, :), ...
             'LineWidth', 1.2, 'DisplayName', sprintf('C%d', i));
    end
    yline(p.consumer.T_r_min, ':k');
    grid on; legend('Location', 'best');
    ylabel('T_s^{(i)} [C]'); xlabel('time [min]');
    title('Supply temperature received per consumer');

    % [2,1] Demand vs delivered
    nexttile;
    hold on;
    for i = 1:n_user
        plot(t_min, sc.d_i(i, :)/1e3, '-',  'Color', cmap(i, :), ...
             'LineWidth', 1.2, 'DisplayName', sprintf('d^{(%d)}', i));
        plot(t_min, sc.c_i(i, :)/1e3, '--', 'Color', cmap(i, :), ...
             'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    grid on;
    legend('Location', 'best');
    ylabel('power [kW]'); xlabel('time [min]');
    title('Demand d^{(i)} (solid) vs delivered c^{(i)} (dashed)');

    % [2,2] Per-consumer flow + trunk
    nexttile;
    hold on;
    for i = 1:n_user
        plot(t_min, sc.q_users(i, :), '-', 'Color', cmap(i, :), ...
             'LineWidth', 1.2, 'DisplayName', sprintf('q^{(%d)}', i));
    end
    plot(t_min, sc.q_F0_F1, 'k--', 'LineWidth', 1.5, ...
         'DisplayName', 'q_{F_0 \rightarrow F_1}');
    grid on;
    legend('Location', 'best');
    ylabel('flow [kg/s]'); xlabel('time [min]');
    title('Per-consumer + trunk flow');

    fname = fullfile(outdir, sprintf('demo_%s', sc.name));
    exportgraphics(fig, [fname '.pdf'], 'ContentType', 'vector');
    exportgraphics(fig, [fname '.png'], 'Resolution', 200);
    close(fig);
    fprintf('  saved %s.{pdf,png}\n', fname);
end

%% Summary figure
% One bar/table figure capturing physical invariants per scenario.
n_scen = numel(scenarios);
max_cd  = zeros(1, n_scen);
min_Tr  = zeros(1, n_scen);
KCL_res = zeros(1, n_scen);
met_pct = zeros(1, n_scen);
labels  = strings(1, n_scen);
for k = 1:n_scen
    sc = scenarios(k);
    max_cd(k)  = max(sc.c_i(:) - sc.d_i(:));
    min_Tr(k)  = min(sc.T_r_i(:));
    n_last     = max(1, numel(sc.t)-4):numel(sc.t);
    KCL_res(k) = max(abs(sum(sc.q_users(:, n_last), 1) - sc.q_F0_F1(n_last)));
    met_pct(k) = 100 * sum(sc.c_i(:)) / max(sum(sc.d_i(:)), 1e-9);
    labels(k)  = strrep(sc.name, '_', ' ');
end

fig = figure('Visible', 'off', 'Position', [100 100 900 560]);
tl  = tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tl, 'Plant validation: physical invariants across scenarios', ...
      'FontWeight', 'bold');

nexttile;
bar(max_cd); set(gca, 'XTickLabel', labels, 'XTickLabelRotation', 30);
yline(0, '-k'); grid on; ylabel('max(c^{(i)} - d^{(i)}) [W]');
title('Substation: c \leq d');

nexttile;
bar(min_Tr); hold on;
yline(p.consumer.T_r_min, ':r', sprintf('T_r^{min} = %g', p.consumer.T_r_min));
set(gca, 'XTickLabel', labels, 'XTickLabelRotation', 30);
grid on; ylabel('min T_r^{(i)} [C]');
title('Substation: T_r \geq T_r^{min}');

nexttile;
bar(KCL_res); set(gca, 'XTickLabel', labels, 'XTickLabelRotation', 30);
set(gca, 'YScale', 'log'); grid on; ylabel('Kirchhoff residual [kg/s]');
title('Steady-state Kirchhoff: q_{F_0F_1} = \Sigma q^{(i)}');

nexttile;
bar(met_pct); set(gca, 'XTickLabel', labels, 'XTickLabelRotation', 30);
ylim([min(min(met_pct), 95)-1, 102]);
grid on; ylabel('total met [%]');
title('Total demand met (informational)');

fname = fullfile(outdir, 'demo_invariants');
exportgraphics(fig, [fname '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [fname '.png'], 'Resolution', 200);
close(fig);
fprintf('  saved %s.{pdf,png}\n', fname);
end
