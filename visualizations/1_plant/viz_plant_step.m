% VIZ_PLANT_STEP  Supply-temperature step response.
%
% I step the source supply temperature up by +6 C and watch the warm front arrive at
% each consumer. The point is the transport physics: each consumer reacts after a
% delay set by its distance from the source, and the rise is smoothed by the pipe
% thermal mass. Transit time is distance over flow speed; at nominal flow the F0->F1
% trunk is about 25 min and the furthest consumer (C5, ~1.9 km) about 90 min, and this
% run uses 2x flow so the delays you read off are roughly half those. What matters is
% the ordering: later consumers lag the nearer ones, the signature of 1D pipe
% transport, which means the plant moves heat through the network correctly.
%
% Runs live against the AROMA plant; nothing is loaded from a saved .mat.

clear; clc;

% --- put the aroma-build code on the path (same setup the demos use) ---
here = fileparts(mfilename('fullpath'));   % .../visualizations/1_plant
root = fileparts(fileparts(here));         % .../aroma-build (repo root, two levels up)
addpath(root); startup();

p = params();          % all physical constants (cp, Tin_nom, pipe + source params, ...)
p.Ts = 60;             % read-out cadence [s]; the plant ODE is continuous, 60 s just gives smooth curves

% --- build the plant (cold steady state at Tin_nom) ---
[net, z0_cold] = build_plant(p);
ei     = edge_user_index(net);   % index helper: which edges/nodes are the consumers
n_user = ei.n_user;              % number of consumers (5 in AROMA)
mdotE  = net.mdotEdges(:);       % nominal per-edge mass flow [kg/s] = the design operating point
F0_idx = find(strcmp({net.Nodes.name}, 'F0'));   % source supply node, for the actual source temperature

% Isolate the pure transport response: set the consumer demand to zero so the
% supply temperature is driven only by the source step and pipe transport, not
% by load-driven variation. This is the clean way to read a step response.
% (The normal-operation script runs WITH the full demand.)
for i = 1:n_user
    ci = find(strcmp({net.Nodes.name}, ei.consumers{i}));
    net.Nodes(ci).params.Q_peak = 0;
end

% --- step settings ---
t_start = 12.5 * 3600;   % wall-clock start [s] (12:30); only fixes the demand phase
T_warm  = 6   * 3600;    % warmup [s]; long enough that even the slow far branches are
                         % fully settled before the step, so the response is a clean rise
T_sim   = 8   * 3600;    % how long we watch the response [s] (long enough for the far consumers)
step_K  = 6;             % supply-temperature step [C]: T_0s goes Tin_nom -> Tin_nom + 6

% Hold the flow constant so the ONLY change is the supply step. We run a bit
% above nominal (2x) so transport dominates over pipe heat loss on the
% low-flow far branches, and every consumer's response is clearly visible
% within the window (at nominal flow the furthest branch is transport-slow).
flow_factor = 2.0;
r_q_fun = @(t) flow_factor * mdotE;

% --- 1) warm up at nominal supply (offset u = 0) so we step from steady state ---
p_wu = p; p_wu.t_offset = t_start; p_wu.r_q_fun = r_q_fun;
res_wu = simulate_plant(net, z0_cold, p_wu, @(t) 0, @(t) 1.0, T_warm);

% --- 2) apply the +step_K supply step and simulate the response ---
%   u_fun returns the supply offset above Tin_nom, so @(t) step_K is a +6 C step.
net_run = net; net_run.q0 = res_wu.q_edges(:, end);   % carry the flow state across the warmup
p_run = p; p_run.t_offset = t_start + T_warm; p_run.r_q_fun = r_q_fun;
res = simulate_plant(net_run, res_wu.z_final, p_run, @(t) step_K, @(t) 1.0, T_sim);

t_h      = res.t / 3600;             % time since the step [h]
T_0s_set = p.Tin_nom + res.u;        % the supply-temperature setpoint (the +6 step command) [C]
T_F0     = res.Tout(F0_idx, :);      % the actual source supply temperature the plant produces [C]

% --- figure: the step at the source + each consumer's supply temperature ---
%   setpoint  = the +6 command; source = what the CSTR source actually produces;
%   consumers = what arrives after pipe transport (always <= the source).
figure('Color', 'w');
plot(t_h, T_0s_set, ':', 'Color', [.55 .55 .55], 'LineWidth', 1.3); hold on;   % setpoint (+6)
plot(t_h, T_F0, 'k', 'LineWidth', 1.8);                                         % actual source supply
plot(t_h, res.T_s_i', 'LineWidth', 1.2);                                        % each consumer
grid on; xlabel('time since step [h]'); ylabel('supply temperature [\circC]');
leg = [{'setpoint (+6)', 'source supply'}, arrayfun(@(i) sprintf('C%d', i), 1:n_user, 'UniformOutput', false)];
legend(leg, 'Location', 'southeast');
title(sprintf('Supply-temperature step (+%d C): it reaches each consumer after its transport delay', step_K));

% WHAT YOU SEE: the source steps up; each consumer's supply temperature rises
% later and more gradually the further it sits from the source (the two nearest,
% C1 and C4, rise first and overlap; the far C3/C5 last), and never exceeds the
% source. That delay ordering is the 1D pipe-transport physics.

exportgraphics(gcf, fullfile(here, 'viz_plant_step.pdf'), 'ContentType', 'vector');
fprintf('Saved viz_plant_step.pdf\n');
