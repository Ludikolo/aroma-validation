% GENERATE_OPERATIONAL  A few smooth operational days to augment the PRBS
% training set. Nominal flow, the diurnal demand profile at staggered start
% times. PRBS gives the broad excitation for identifiability; these calibrate
% the predictor on the smooth operating regime PRBS under-represents, which
% removes the peak under-prediction at the applied control step.

clear; clc;
startup;
p = params();

outdir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
[net, z0] = build_plant(p);
mdotE = net.mdotEdges(:);
R0 = find(strcmp({net.Nodes.name}, 'R0'));
F0 = find(strcmp({net.Nodes.name}, 'F0'));

phases_h = 0:3:21;                     % 8 days, every 3 h over the daily cycle
for j = 1:numel(phases_h)
    pr = p; pr.t_offset = phases_h(j) * 3600; pr.r_q_fun = @(t) mdotE;
    res = simulate_plant(net, z0, pr, @(t) 0, @(t) 1.0, p.data.traj_dur_s + 2 * 3600);

    traj = struct();
    traj.t        = res.t;
    traj.t_offset = phases_h(j) * 3600;
    traj.T_0s     = p.Tin_nom + res.u;
    traj.r_q      = res.r_q;
    traj.T_ir     = res.T_r_i;
    traj.T_is     = res.T_s_i;
    traj.T_0r     = res.Tout(R0, :);
    traj.T_F0     = res.Tout(F0, :);
    traj.q_users  = res.q_users;
    traj.q_edges  = res.q_edges;
    traj.Tout     = res.Tout;
    traj.d        = res.d_i;
    traj.theta    = (res.T_s_i - res.T_r_i) .* res.q_users;

    save(fullfile(outdir, sprintf('op_%02d.mat', j)), 'traj');
end
fprintf('Wrote %d operational days op_01..op_%02d.mat\n', numel(phases_h), numel(phases_h));
