function [res, res_hold, t_h, ei, p] = run_controller_scenario(mdot_scale)
% RUN_CONTROLLER_SCENARIO  Run one 24 h KPC closed loop (and a hold-nominal
% reference) at a given network capacity, for the controller figures.
%
%   [res, res_hold, t_h, ei, p] = run_controller_scenario(mdot_scale)
%
% mdot_scale = 1.0 is the nominal (easy) case. mdot_scale < 1 throttles every
% pipe, so the same demand has to be met through a tighter network (the hard
% case). The controller (V_seq predictor + locked tune, shipped in this repo) and
% the warmup are identical in both cases, so the only thing that changes between
% easy and hard is the capacity.
%
% res      = the KPC closed-loop trajectory (one column per Ts step)
% res_hold = hold-nominal on the same plant + warmup (apply nominal flow, source
%            at nominal): the "no controller" reference
% t_h      = time axis [h]

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(root); startup();

p = params();
pred_dir = fullfile(root, 'predictor_open_loop', 'results');
ctrl_dir = fullfile(root, 'controller_closed_loop', 'results');

% build the controller: V_seq predictor maps + locked tune
F  = load(fullfile(pred_dir, 'vseq_fits_full.mat'));
Hb = load(fullfile(ctrl_dir, 'best_horizons.mat'));
Ab = load(fullfile(ctrl_dir, 'best_alpha.mat'));
Bb = load(fullfile(ctrl_dir, 'best_tune.mat'));
pred = truncate_pred_to_Np(build_kpc_v2_matrices(F.fits), Hb.Np_best);
tune = Bb.best_tune;
tune.dT_0s_max = 7.5; tune.dr_q_max = 4.5; tune.dT_ir_max = 15.0;   % actuator rate caps
tune.use_kirchhoff = true; tune.use_Tr_le_Ts = true; tune.use_mixing = true; tune.rho_slack_mix = 1;
tune.Nc = Hb.Nc_best; tune.alpha_energy = Ab.alpha_star;

% plant at the requested capacity (scale every pipe by mdot_scale)
T_warm = 30 * 60; T_sim = 24 * 3600;
[net, z0] = build_plant(p);
net.mdotEdges = mdot_scale * net.mdotEdges;
net.flows     = mdot_scale * net.flows;
net.q0        = mdot_scale * net.q0;

ei = edge_user_index(net); n_user = ei.n_user;
F0 = find(strcmp({net.Nodes.name}, 'F0'));
R0 = find(strcmp({net.Nodes.name}, 'R0'));
con = arrayfun(@(i) find(strcmp({net.Nodes.name}, ei.consumers{i})), 1:n_user);

% warmup, then the KPC closed loop
p_wu = p; p_wu.t_offset = 0;
res_wu = simulate_plant(net, z0, p_wu, @(t) 0, @(t) 1.0, T_warm);
res = kpc_step_loop(net, res_wu, p, 0, T_warm, T_sim, pred, tune, ei, F0, R0, con);

% hold-nominal on the same plant + warmup: apply nominal flow, source at nominal
p_hold = p; p_hold.t_offset = T_warm;
res_hold = simulate_plant(net, res_wu.z_final, p_hold, @(t) 0, @(t) 1.0, T_sim);
res_hold.d_i = res_hold.d_i(:, 2:end);   % drop the warmup-boundary sample so it
res_hold.c_i = res_hold.c_i(:, 2:end);   % lines up with the receding-horizon window

t_h = (0:size(res.c_i, 2) - 1) * p.Ts / 3600;
end
