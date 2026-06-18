% GENERATE_DATA  Generate the PRBS train/val/test trajectories used to fit
% and evaluate the predictor. Each trajectory is a 24 h PRBS-excited run on
% the AROMA plant at the sample time Ts set in params.m. The set sizes
% (n_train / n_val / n_test) and trajectory duration come from p.data.* in
% params.m. Trajectories are saved as individual .mat files in results/.
%
% THEORY index:
%   line 31: data split

clear; clc;
startup;
p = params();

[net, z0_cold] = build_plant(p);

% warm up the network so each trajectory starts near thermal equilibrium
fprintf('Warming up plant (4h at nominal)...\n');
r_wu = simulate_plant(net, z0_cold, p, @(t) 0, @(t) 1.0, 4*3600);
z0_warm = r_wu.z_final;

outdir = fullfile(fileparts(mfilename('fullpath')), '..', 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end

T_sim     = p.data.traj_dur_s;
seed_base = 1000;       % keep this seed space well apart from the 100s-based one

% Spread t_offsets evenly across 24 h so the predictor sees every
% demand regime (morning ramp, midday peak, evening peak, night).
% Deterministic per-index for reproducibility.
T_day = 24 * 3600;
% THEORY (data split): 20 train / 2 val / 3 test, staggered day offsets so the demand phase differs per trajectory
offsets_train = mod((0:p.data.n_train-1) * T_day / p.data.n_train, T_day);
offsets_val   = mod((0:p.data.n_val-1)   * T_day / p.data.n_val   + 0.5*T_day/p.data.n_train, T_day);
offsets_test  = mod((0:p.data.n_test-1)  * T_day / p.data.n_test  + 0.7*T_day/p.data.n_train, T_day);

fprintf('\nTraining set (%d trajectories, %.1f h each, t_offset spread over 24 h)...\n', ...
    p.data.n_train, T_sim/3600);
for j = 1:p.data.n_train
    seed = seed_base + j;
    p_traj = p;
    p_traj.t_offset = offsets_train(j);
    fprintf('  train_%02d  seed=%d  t_offset=%5.2f h ... ', j, seed, offsets_train(j)/3600);
    tic;
    traj = generate_traj(net, z0_warm, p_traj, T_sim, seed);
    traj.t_offset = offsets_train(j);
    save(fullfile(outdir, sprintf('train_%02d.mat', j)), 'traj');
    fprintf('done (%.1f s)\n', toc);
end

fprintf('\nValidation set (%d trajectories)...\n', p.data.n_val);
for j = 1:p.data.n_val
    seed = seed_base + 1000 + j;
    p_traj = p;
    p_traj.t_offset = offsets_val(j);
    fprintf('  val_%02d  seed=%d  t_offset=%5.2f h ... ', j, seed, offsets_val(j)/3600);
    tic;
    traj = generate_traj(net, z0_warm, p_traj, T_sim, seed);
    traj.t_offset = offsets_val(j);
    save(fullfile(outdir, sprintf('val_%02d.mat', j)), 'traj');
    fprintf('done (%.1f s)\n', toc);
end

fprintf('\nTest set (%d trajectories)...\n', p.data.n_test);
for j = 1:p.data.n_test
    seed = seed_base + 2000 + j;
    p_traj = p;
    p_traj.t_offset = offsets_test(j);
    fprintf('  test_%02d  seed=%d  t_offset=%5.2f h ... ', j, seed, offsets_test(j)/3600);
    tic;
    traj = generate_traj(net, z0_warm, p_traj, T_sim, seed);
    traj.t_offset = offsets_test(j);
    save(fullfile(outdir, sprintf('test_%02d.mat', j)), 'traj');
    fprintf('done (%.1f s)\n', toc);
end

fprintf('\nAll trajectories saved to %s\n', outdir);
