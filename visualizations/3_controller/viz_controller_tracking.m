% VIZ_CONTROLLER_TRACKING  Does the controller deliver the demand, easy day vs hard?
%
% I run the KPC controller in closed loop for a full day and, per consumer, lay the
% delivered heat (solid blue) under the demand (dotted black), with the hold-nominal
% delivery (dashed grey) as a no-controller reference. Done twice: nominal capacity
% (easy) and half capacity (hard, mdot = 0.50). On the easy day demand, KPC and hold
% all overlap at 100 %; on the hard day the hold reference falls short on the big
% consumers while KPC still hugs the demand at 100 %. The contrast is the point: same
% controller, much tighter network, still on demand. Each title gives the daily met %.
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
    [res, res_hold, t_h, ei, ~] = run_controller_scenario(c.mdot);
    n_user = ei.n_user;

    % per-consumer demand met over the day = total delivered / total demanded
    met = 100 * sum(res.c_i, 2) ./ max(sum(res.d_i, 2), 1e-9);

    figure('Color', 'w', 'Position', [100 100 1000 560]);
    for i = 1:n_user
        subplot(2, 3, i);
        plot(t_h, res.d_i(i, :)/1e3, 'k:', 'LineWidth', 1.6); hold on;          % demand
        plot(t_h, res_hold.c_i(i, :)/1e3, '--', 'Color', [.55 .55 .55], 'LineWidth', 1.2);  % hold (no control)
        plot(t_h, res.c_i(i, :)/1e3, 'b', 'LineWidth', 1.6);                    % KPC delivered
        grid on; xlabel('time [h]'); ylabel('heat [kW]'); xlim([0 24]);
        title(sprintf('C%d   (KPC met %.1f%%)', i, met(i)));
        if i == 1, legend('demand', 'hold', 'KPC', 'Location', 'best'); end
    end
    subplot(2, 3, 6);
    plot(t_h, sum(res.d_i, 1)/1e3, 'k:', 'LineWidth', 1.6); hold on;
    plot(t_h, sum(res_hold.c_i, 1)/1e3, '--', 'Color', [.55 .55 .55], 'LineWidth', 1.2);
    plot(t_h, sum(res.c_i, 1)/1e3, 'b', 'LineWidth', 1.6);
    grid on; xlabel('time [h]'); ylabel('heat [kW]'); xlim([0 24]);
    title(sprintf('whole network   (KPC met %.1f%%)', 100*sum(res.c_i(:))/sum(res.d_i(:))));

    sgtitle(sprintf('Demand tracking, %s', c.descr));
    out = fullfile(here, sprintf('viz_controller_tracking_%s.pdf', c.tag));
    exportgraphics(gcf, out, 'ContentType', 'vector');
    fprintf('Saved %s  (per-consumer KPC met: %s %%)\n', ...
            sprintf('viz_controller_tracking_%s.pdf', c.tag), mat2str(round(met', 1)));
end

% WHAT YOU SEE: on the easy day all three curves overlap (even hold meets demand);
% on the hard day the throttled network makes hold drop below demand on the big
% consumers while KPC stays on the dotted demand line at 100 %.
