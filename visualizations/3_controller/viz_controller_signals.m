% VIZ_CONTROLLER_SIGNALS  What the controller actually does, easy day vs hard day.
%
% For the same easy (nominal) and hard (mdot = 0.50) closed loops as the tracking
% figure, two figures: a hydraulic one (demand, commanded flow r_q, realized flow q,
% network total) and a thermal one (source supply/return temps with their [14,26] C
% band, consumer supply/return temps with the 15 C return floor, and the QP solve
% time). On the easy day the source temp barely moves and returns sit well above the
% floor; on the hard day the controller pushes the source supply to its 26 C cap,
% saturates flows, and lets the returns fall toward (never through) the floor, still
% meeting demand. T_ir is the realized return temp (the substation produces it); the
% controller actuates only the source supply temp and the flows.
%
% Runs live (see run_controller_scenario); no precomputed results are loaded.

clear; clc;
here = fileparts(mfilename('fullpath'));

cases = struct( ...
    'tag',   {'easy',              'hard'}, ...
    'mdot',  {1.00,                0.50}, ...
    'descr', {'nominal capacity',  'stressed capacity (mdot = 0.50)'});

for ci = 1:numel(cases)
    c = cases(ci);
    [res, ~, t_h, ei, p] = run_controller_scenario(c.mdot);
    n_user = ei.n_user;
    cons = compose('C%d', 1:n_user);

    % ---- FIGURE 1: hydraulic (demand -> commanded flow -> realized flow -> total) ----
    figure('Color', 'w', 'Position', [100 100 950 650]);

    subplot(2, 2, 1);
    plot(t_h, res.d_i'/1e3, 'LineWidth', 1.1);          % per-consumer demand [kW], the driver
    grid on; xlabel('time [h]'); ylabel('demand [kW]');
    legend(cons, 'Location', 'best'); title('Per-consumer demand d^{(i)}');

    subplot(2, 2, 2);
    plot(t_h, res.r_q(ei.user, :)', 'LineWidth', 1.1);  % commanded stub flow [kg/s]
    grid on; xlabel('time [h]'); ylabel('flow [kg/s]');
    legend(cons, 'Location', 'best'); title('Commanded flow r_q^{(i)}');

    subplot(2, 2, 3);
    plot(t_h, res.q_users', 'LineWidth', 1.1);          % realized stub flow [kg/s] (same C order)
    grid on; xlabel('time [h]'); ylabel('flow [kg/s]');
    legend(cons, 'Location', 'best'); title('Realized flow q^{(i)} (follows r_q within one step)');

    subplot(2, 2, 4);
    plot(t_h, res.Q_net, 'k', 'LineWidth', 1.5);        % total network flow [kg/s]
    grid on; xlabel('time [h]'); ylabel('flow [kg/s]');
    title('Network total flow Q_{net}');

    sgtitle(sprintf('Hydraulics, %s', c.descr));
    exportgraphics(gcf, fullfile(here, sprintf('viz_controller_hydraulic_%s.pdf', c.tag)), 'ContentType', 'vector');

    % ---- FIGURE 2: thermal + feasibility ----
    figure('Color', 'w', 'Position', [100 100 950 650]);

    subplot(2, 2, 1);
    % T_0s is the source supply-temperature COMMAND (the actuator), bounded to
    % Tin_nom +/- 6 = [14, 26] C. T_F0 is the temperature the source actually
    % produces: it can sit a little above the command at peak load because the warm
    % return water blends through the source, so the source is not an ideal
    % temperature source. T_0r is the mixed return; T_F0 - T_0r is the heat injected.
    plot(t_h, res.T_0s, 'k', 'LineWidth', 1.5); hold on;
    plot(t_h, res.T_F0, 'b', 'LineWidth', 1.2);
    plot(t_h, res.T_0r, 'Color', [.55 .55 .55], 'LineWidth', 1.3);
    yline(14, 'r--', 'LineWidth', 1.0); yline(26, 'r--', 'LineWidth', 1.0);
    grid on; xlabel('time [h]'); ylabel('[\circC]'); ylim([10 29]);
    legend('T_{0s} command', 'T_{F0} feed', 'T_{0r} return', 'Location', 'best');
    title('Source temperatures (command T_{0s} bounded 14-26 C)');

    subplot(2, 2, 2);
    % consumer supply temps; they follow the source feed T_F0 (delayed by transport),
    % so at peak load they can sit slightly above the 26 C source command bound too.
    plot(t_h, res.T_is', 'LineWidth', 1.1);
    grid on; xlabel('time [h]'); ylabel('[\circC]');
    legend(cons, 'Location', 'best'); title('Consumer supply temperature T_s^{(i)}');

    subplot(2, 2, 3);
    plot(t_h, res.T_ir', 'LineWidth', 1.1); hold on;    % realized consumer return temps
    yline(15, 'r--', '15 C return floor', 'LineWidth', 1.2);
    grid on; xlabel('time [h]'); ylabel('[\circC]');
    legend(cons, 'Location', 'best'); title('Consumer return temperature T_r^{(i)} (floor 15 C)');

    subplot(2, 2, 4);
    % res.solve_ms = wall-clock QP solve time per step; the step budget is Ts = 900 s
    stem(t_h, res.solve_ms, 'b', 'Marker', 'none', 'LineWidth', 1.0); hold on;
    yline(median(res.solve_ms), 'r--', 'LineWidth', 1.2);
    grid on; xlabel('time [h]'); ylabel('QP solve time [ms]');
    ylim([0, max(res.solve_ms) * 1.2]);
    title(sprintf('Solve time per step (median %.0f ms, step = %.0f s)', median(res.solve_ms), p.Ts));

    sgtitle(sprintf('Temperatures and feasibility, %s', c.descr));
    exportgraphics(gcf, fullfile(here, sprintf('viz_controller_thermal_%s.pdf', c.tag)), 'ContentType', 'vector');

    fprintf('Saved viz_controller_hydraulic_%s.pdf and viz_controller_thermal_%s.pdf\n', c.tag, c.tag);
end
