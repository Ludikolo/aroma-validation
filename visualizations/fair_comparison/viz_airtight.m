% VIZ_AIRTIGHT  Renders the shipped comparison run itself.
%
% Three one-message figures straight from airtight_comparison_ls.mat (the run
% behind the thesis table; run run_comparison_airtight.m to regenerate it):
%   viz_airtight.pdf            comfort + speed head-to-head, binding consumer
%   viz_airtight_optimality.pdf per-step solver certificate: the convex QP is
%                               feasible and optimal on every step, NMPC is not
%   viz_airtight_npsweep.pdf    shared-horizon fairness: worst-met% vs Np for
%                               the three tuned controllers (from sweep_cheap_ls.mat),
%                               including the honest point at Np = 8 where the
%                               iterated baseline beats KPC

clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(fileparts(here))); startup();

A = load(fullfile(here, 'airtight_comparison_ls.mat'));
p = A.p; names = A.names;
R  = {A.res_kpc, A.res_klmpc, A.res_jlmpc, A.res_nmpc};
cc = {[0.10 0.30 0.80], [0.20 0.60 0.20], [0.85 0.50 0.10], [0.60 0.10 0.60]};
N  = size(A.res_kpc.c_i, 2); th = (0:N-1) * p.Ts / 3600;

% binding consumer = lowest demand-met across the four controllers
pcmet  = @(r) 100 * sum(r.c_i, 2) ./ max(sum(r.d_i, 2), 1e-9);
allmet = cell2mat(cellfun(@(r) pcmet(r), R, 'UniformOutput', false));
[~, bc] = min(min(allmet, [], 2));

%% comfort + speed
figure('Color', 'w', 'Position', [100 100 1000 700]);
subplot(2, 1, 1);
plot(th, A.res_kpc.d_i(bc, :)/1e3, 'k:', 'LineWidth', 2.4); hold on;
for j = 1:4
    plot(th, R{j}.c_i(bc, :)/1e3, 'Color', cc{j}, 'LineWidth', 1.3 + 0.9*(j==1));
end
grid on; xlabel('time [h]'); ylabel('delivered heat [kW]'); xlim([0 24]);
legend(['demand', names], 'Location', 'best');
title(sprintf('Delivered vs demand, binding consumer C%d (stressed day, mdot = %.2f)', bc, A.mdot_scale));

subplot(2, 1, 2);
for j = 1:4
    semilogy(th, max(R{j}.solve_ms, 1e-1), 'Color', cc{j}, 'LineWidth', 1.2 + 0.8*(j==1)); hold on;
end
yline(p.Ts * 1000, 'r--', 'Ts = 900 s budget', 'LineWidth', 1.0);
grid on; xlabel('time [h]'); ylabel('solve time [ms]'); xlim([0 24]);
legend(names, 'Location', 'best');
title('Solve time per step');
exportgraphics(gcf, fullfile(here, 'viz_airtight.pdf'), 'ContentType', 'vector');

%% per-step solver certificate: KPC vs NMPC
figure('Color', 'w', 'Position', [100 100 900 330]);
hold on;
spec = {1, 'o', [0.15 0.55 0.15]; 2, '^', [0.90 0.55 0.10]; -2, 'x', [0.80 0.10 0.10]};
rows = {A.res_kpc.exitflag, 2, 'KPC (convex QP)'; A.res_nmpc.exitflag, 1, 'NMPC (SQP)'};
for r = 1:2
    ef = rows{r, 1}; y = rows{r, 2};
    for s = 1:size(spec, 1)
        k = find(ef == spec{s, 1});
        if isempty(k), continue; end
        plot(th(k), y * ones(size(k)), spec{s, 2}, 'Color', spec{s, 3}, ...
             'MarkerSize', 5, 'LineWidth', 1.1, 'HandleVisibility', 'off');
    end
end
h1 = plot(nan, nan, 'o', 'Color', spec{1, 3}, 'LineWidth', 1.1);
h2 = plot(nan, nan, '^', 'Color', spec{2, 3}, 'LineWidth', 1.1);
h3 = plot(nan, nan, 'x', 'Color', spec{3, 3}, 'LineWidth', 1.1);
legend([h1 h2 h3], {'optimal (exitflag 1)', 'StepTolerance stall (2)', 'no feasible point (-2)'}, ...
       'Location', 'eastoutside');
ylim([0.5 2.5]); yticks([1 2]); yticklabels({'NMPC', 'KPC'});
xlim([0 24]); xlabel('time [h]'); grid on;
title(sprintf('Solver certificate per step: KPC %d/%d optimal, NMPC %d/%d', ...
      sum(A.res_kpc.exitflag == 1), N, sum(A.res_nmpc.exitflag == 1), N));
exportgraphics(gcf, fullfile(here, 'viz_airtight_optimality.pdf'), 'ContentType', 'vector');

%% shared-horizon fairness sweep
W = load(fullfile(here, 'sweep_cheap_ls.mat'));
nps = [W.swp.Np];
figure('Color', 'w', 'Position', [100 100 640 420]);
plot(nps, [W.swp.kpc],   '-o', 'LineWidth', 1.8, 'Color', cc{1}); hold on;
plot(nps, [W.swp.klmpc], '-s', 'LineWidth', 1.6, 'Color', cc{2});
plot(nps, [W.swp.jlmpc], '-^', 'LineWidth', 1.6, 'Color', cc{3});
grid on; xticks(nps);
xlabel('prediction horizon N_p'); ylabel('worst-consumer demand met [%]');
legend(names(1:3), 'Location', 'southeast');
title('Shared horizon N_p = 12 does not favour KPC');
exportgraphics(gcf, fullfile(here, 'viz_airtight_npsweep.pdf'), 'ContentType', 'vector');

%% headline
fprintf('\nFour-controller comparison (mdot = %.2f, Np = %d):\n', A.mdot_scale, A.Np);
fprintf('%-15s  worst-met%%  med-solve  optimal\n', 'controller');
for j = 1:4
    m = A.M.(matlab.lang.makeValidName(names{j}));
    fprintf('%-15s  %8.3f  %8.0f ms  %3d/%d\n', names{j}, m.worst, m.med_ms, ...
            sum(R{j}.exitflag == 1), m.n);
end
fprintf('Saved viz_airtight.pdf, viz_airtight_optimality.pdf, viz_airtight_npsweep.pdf\n');
