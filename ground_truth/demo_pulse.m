% DEMO_PULSE  Transport-delay demonstration on the AROMA plant.
%   Puts a +5 K block pulse of one hour on the source supply command
%   and watches the warm front travel through the supply network to
%   the five consumer taps. Per consumer the measured arrival time is
%   checked against the transport delay computed from the network
%   itself (pipe volume over volumetric flow, summed along the supply
%   path), the peak response has to attenuate with distance, and the
%   whole network has to settle back to its baseline afterwards.
%
%   See plant_validation.pdf, Section 3 (the pulse and delay-parity figures).
%
%   Timing. The terminal branch pipes are slow: F6 -> F7 carries only
%   0.12 kg/s through an 83 mm pipe, a 6.2 h transit at nominal flow.
%   C5 therefore sits about 8.7 h downstream of the source via its
%   faster feed (F0-F1-F2-F3-F4-F7-F8) and 9.8 h via the other one,
%   far beyond the 94 min a naive trunk-velocity estimate suggests.
%   The run is 25 h: 1 h pre-pulse baseline, 1 h pulse, then enough
%   tail for the pulse to clear C5 and wash out. The warmup is a full
%   24 h (generate_data uses 4 h) for the same reason: the slow
%   branch pipes need most of a day to flush the cold-start water, so
%   a shorter warmup would leave the far consumers off their diurnal
%   steady behaviour during the baseline hour.
%
%   Baseline. Demand follows the standard diurnal profiles, which are
%   exactly 24 h periodic. With a 25 h run the final hour covers the
%   same time of day as the pre-pulse baseline hour, so the
%   return-to-baseline check compares like with like and the diurnal
%   drift cancels out.
%
%   Amplitude. The source is not a Dirichlet boundary. The CSTR
%   mixing rule T_int = (1-theta1) T_in + theta1 T_out with
%   theta1 = 0.7 gives the outlet a steady gain of
%   kQ / (cp q + kQ) / theta1 ~ 1.09 on the command, so the commanded
%   5 K enters the network as roughly 5.4 K at F0. The consumer
%   amplitudes are therefore checked against the pulse as realised at
%   the source outlet, and that realised pulse is itself checked
%   against the commanded one.

clear; clc;
startup;
p = params();

pulse_amp = 5;           % [K] block pulse on the source command
t_on      = 1 * 3600;    % [s] pulse start, after 1 h of baseline
t_off     = 2 * 3600;    % [s] pulse end
T_sim     = 25 * 3600;   % [s] total run, see header
T_warm    = 24 * 3600;   % [s] full-day warmup, see header

fprintf('\n=== AROMA source-pulse demo ===\n');

fprintf('Building plant...\n');
[net, z0_cold] = build_plant(p);
ei     = edge_user_index(net);
n_user = ei.n_user;
F0_idx = find(strcmp({net.Nodes.name}, 'F0'));

% Warmup at nominal inputs. Ts = 900 s is fine here: the pipe ODE runs
% in continuous time either way and nothing changes between samples.
fprintf('Warming up plant (24 h at nominal)...\n');
res_wu = simulate_plant(net, z0_cold, p, @(t) 0, @(t) 1.0, T_warm);

% Pulse run at the fine sample rate so arrival detection has 60 s
% resolution. Nominal constant flows, standard diurnal demand, wall
% clock continues where the warmup stopped (midnight).
p_run = p;
p_run.Ts       = 60;
p_run.t_offset = T_warm;
u_fun = @(t) pulse_amp * (t >= t_on & t < t_off);
fprintf('Running the pulse experiment (+%g K on the command from %g to %g min, 25 h total)...\n', ...
    pulse_amp, t_on/60, t_off/60);
res = simulate_plant(net, res_wu.z_final, p_run, u_fun, @(t) 1.0, T_sim);

t     = res.t;
T_F0  = res.Tout(F0_idx, :);     % realised source outlet temperature
T_s_i = res.T_s_i;               % supply temperature at each consumer tap

%% Expected transport delay from the network itself
% Delay along a supply path = sum over its edges of pipe volume over
% volumetric flow, using the realised edge flows of this run. C5 is
% fed through both F4->F7 and F6->F7 (alpha_split = 0.5); the first
% arrival comes via the faster of its two paths, so take the min.
pre    = t < t_on;                       % pre-pulse baseline window
q_mean = mean(res.q_edges(:, pre), 2);   % realised edge flows [kg/s]

paths = { ...
    {'F0','F1','F2','C1'}, ...
    {'F0','F1','F2','F3','C2'}, ...
    {'F0','F1','F2','F3','F4','F5','C3'}, ...
    {'F0','F1','F6','C4'}, ...
    {'F0','F1','F2','F3','F4','F7','F8','C5'}};
tau_comp = zeros(n_user, 1);
for i = 1:n_user
    tau_comp(i) = supply_path_delay(net, paths{i}, q_mean);
end
tau_comp(5) = min(tau_comp(5), ...
    supply_path_delay(net, {'F0','F1','F6','F7','F8','C5'}, q_mean));

%% Measured arrival and peak response
% Arrival = first sample after pulse start where T_s^(i) crosses 10
% percent of the commanded amplitude above its pre-pulse baseline.
base_i  = mean(T_s_i(:, pre), 2);
base_F0 = mean(T_F0(pre));
post    = t >= t_on;
thr     = 0.1 * pulse_amp;

t_arr = zeros(n_user, 1);
A_i   = zeros(n_user, 1);
for i = 1:n_user
    k_hit = find(post & T_s_i(i, :) > base_i(i) + thr, 1);
    assert(~isempty(k_hit), 'pulse never crossed the 10%% threshold at %s', ei.consumers{i});
    t_arr(i) = t(k_hit) - t_on;
    A_i(i)   = max(T_s_i(i, post)) - base_i(i);
end
A_src = max(T_F0(post)) - base_F0;

fprintf('\n  Pulse: %.2f K commanded, %.2f K realised at the F0 outlet\n\n', pulse_amp, A_src);
fprintf('  consumer   computed delay [min]   measured arrival [min]   ratio   peak response [K]\n');
for i = 1:n_user
    fprintf('  %-8s   %20.1f   %22.1f   %5.2f   %17.2f\n', ...
        ei.consumers{i}, tau_comp(i)/60, t_arr(i)/60, t_arr(i)/tau_comp(i), A_i(i));
end

%% Save
outdir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end
pulse = struct('t', t, 'T_F0', T_F0, 'T_s_i', T_s_i, ...
    'T_0s_cmd', p.Tin_nom + res.u, 'consumers', {ei.consumers}, ...
    'pulse_amp', pulse_amp, 't_on', t_on, 't_off', t_off, 'Ts', p_run.Ts, ...
    'base_F0', base_F0, 'base_i', base_i, 'A_src', A_src, 'A_i', A_i, ...
    'tau_comp', tau_comp, 't_arr', t_arr);
save(fullfile(outdir, 'pulse.mat'), 'pulse');
fprintf('\nSaved %s/pulse.mat\n', outdir);

%% Plot: source outlet plus the five consumer taps, pulse window shaded
figdir = fullfile(fileparts(mfilename('fullpath')), 'figures');
if ~exist(figdir, 'dir'), mkdir(figdir); end

fig = figure('Visible', 'off', 'Position', [100 100 900 480]);
hold on;
t_h  = t / 3600;
cmap = lines(n_user);
yl   = [min([T_F0, T_s_i(:)']) - 0.3, max([T_F0, T_s_i(:)']) + 0.3];
patch([t_on t_off t_off t_on]/3600, [yl(1) yl(1) yl(2) yl(2)], ...
      [0.92 0.92 0.92], 'EdgeColor', 'none', 'DisplayName', 'pulse window');
plot(t_h, T_F0, 'k-', 'LineWidth', 1.5, 'DisplayName', 'F0 source outlet');
for i = 1:n_user
    plot(t_h, T_s_i(i, :), '-', 'Color', cmap(i, :), 'LineWidth', 1.2, ...
         'DisplayName', sprintf('C%d', i));
    k_arr = find(t >= t_on + t_arr(i), 1);
    plot(t_h(k_arr), T_s_i(i, k_arr), 'o', 'MarkerSize', 5, ...
         'MarkerFaceColor', cmap(i, :), 'MarkerEdgeColor', 'k', ...
         'HandleVisibility', 'off');
end
grid on; ylim(yl); xlim([0 T_sim/3600]);
xlabel('time [h]'); ylabel('supply temperature [C]');
title(sprintf('+%g K source pulse (1 h) travelling to the consumers', pulse_amp));
legend('Location', 'northeast');
fname = fullfile(figdir, 'pulse');
exportgraphics(fig, [fname '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [fname '.png'], 'Resolution', 200);
close(fig);
fprintf('  saved %s.{pdf,png}\n', fname);

%% Inline validation: four pulse checks
fprintf('\n=== Pulse checks ===\n');

% P1: the measured arrival matches the transport delay the network
% geometry and flows predict. The band is wide because the 10 percent
% threshold triggers on the diffused leading edge (upwind mixing over
% 10 m cells), which runs ahead of the plug-flow transit time.
for i = 1:n_user
    r = t_arr(i) / tau_comp(i);
    assert(r >= 0.7 && r <= 1.2, ...
        'P1 fail: %s arrival %.1f min vs computed %.1f min (ratio %.2f)', ...
        ei.consumers{i}, t_arr(i)/60, tau_comp(i)/60, r);
end
fprintf('P1 pass: measured arrival within [0.7, 1.2] x computed transport delay at every consumer\n');

% P2: every consumer sees the pulse, and nobody sees more than what
% the source actually injected. The realised source pulse itself must
% sit at the commanded amplitude times the ~1.09 CSTR outlet gain.
assert(A_src >= pulse_amp && A_src <= 1.2 * pulse_amp, ...
    'P2 fail: realised source pulse %.2f K outside [1.0, 1.2] x %.1f K', A_src, pulse_amp);
for i = 1:n_user
    assert(A_i(i) > 0, 'P2 fail: %s peak response %.2f K not positive', ei.consumers{i}, A_i(i));
    assert(A_i(i) <= A_src, ...
        'P2 fail: %s peak response %.2f K exceeds the realised source pulse %.2f K', ...
        ei.consumers{i}, A_i(i), A_src);
end
fprintf('P2 pass: every consumer responds and no response exceeds the realised source pulse (%.2f K)\n', A_src);

% P3: pipe mixing spreads the pulse out, so the peak response cannot
% grow with distance. Consumers close together may tie within 0.1 K.
[~, order] = sort(tau_comp);
A_sorted   = A_i(order);
for j = 2:n_user
    assert(A_sorted(j) <= A_sorted(j-1) + 0.1, ...
        'P3 fail: %s (%.2f K) exceeds nearer %s (%.2f K) by more than 0.1 K', ...
        ei.consumers{order(j)}, A_sorted(j), ei.consumers{order(j-1)}, A_sorted(j-1));
end
fprintf('P3 pass: peak response does not increase with computed path delay (attenuation with distance)\n');

% P4: the pulse washes out. The final hour of the run covers the same
% time of day as the baseline hour, so after the pulse has cleared the
% network both hours must look alike.
fin = t >= T_sim - 3600;
res_F0 = abs(mean(T_F0(fin)) - base_F0);
fprintf('  F0  end-vs-baseline residual = %.3f K\n', res_F0);
assert(res_F0 <= 0.2, 'P4 fail: F0 ends %.3f K away from baseline', res_F0);
for i = 1:n_user
    res_i = abs(mean(T_s_i(i, fin)) - base_i(i));
    fprintf('  %-3s end-vs-baseline residual = %.3f K\n', ei.consumers{i}, res_i);
    assert(res_i <= 0.2, 'P4 fail: %s ends %.3f K away from baseline', ei.consumers{i}, res_i);
end
fprintf('P4 pass: source outlet and every consumer supply back within 0.2 K of baseline in the final hour\n');

fprintf('\nAll four pulse checks verified. The pulse travels at the speed the pipe volumes and flows dictate.\n');


%% Local functions
function tau = supply_path_delay(net, names, q_edges)
% Transport delay along a chain of node names: sum over the path's
% edges of pipe volume over volumetric flow (plug-flow transit time).
tau = 0;
for s = 1:numel(names)-1
    e   = find_edge(net, names{s}, names{s+1});
    V   = net.Edges(e).Acs * net.Edges(e).L;   % pipe volume [m^3]
    F   = q_edges(e) / net.rho;                % volumetric flow [m^3/s]
    tau = tau + V / F;
end
end

function e = find_edge(net, from_name, to_name)
for k = 1:numel(net.Edges)
    if strcmp(net.Nodes(net.Edges(k).from).name, from_name) && ...
       strcmp(net.Nodes(net.Edges(k).to).name,   to_name)
        e = k;
        return;
    end
end
error('find_edge: %s -> %s not found', from_name, to_name);
end
