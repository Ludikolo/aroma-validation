function res = lmpc_step_loop(net, res_wu, p, scenario_start, T_warm, T_sim, ...
                               fit_iter, Np, tune, ei, F0_idx, R0_idx, con_idx, lift_fn)
% LMPC_STEP_LOOP  Closed-loop driver for the iterated one-step Koopman LMPC.
%
% This mirrors kpc_step_loop EXACTLY: same warmup fill, same lift to z_0,
% same demand forecast, same one-Ts plant step via simulate_plant, same
% recorded outputs. That is what makes the comparison fair: identical
% plant, scenario, warmup and apply. The only differences are that we
% rebuild the predictor every step (F and G are fixed by A,B, but the
% demand roll-out offsets and the energy-term operating point change
% with the trajectory) and call lmpc_solve instead of kpc_v2_solve.
%
% Instead of a prebuilt `pred` we take `fit_iter` (the iterated A,B) plus
% `Np` (the same prediction horizon KPC uses). The 14th arg `lift_fn` is
% the same lift handle (traj, p) -> [Z, info, meta] KPC uses; defaults to
% candidate_library so the lifted state matches the production controller.

if nargin < 14 || isempty(lift_fn)
    lift_fn = @(t, pp) candidate_library(t, pp);
end

Ts      = p.Ts;
N_cl    = round(T_sim / Ts);
N_wu    = numel(res_wu.t);
n_edges = numel(net.Edges);
n_user  = ei.n_user;

N_tot = N_wu + N_cl;
traj = struct();
traj.t_offset = scenario_start;
traj.t        = zeros(1, N_tot);
traj.Tout     = zeros(numel(net.Nodes), N_tot);
traj.q_edges  = zeros(n_edges, N_tot);
traj.q_users  = zeros(n_user, N_tot);
traj.r_q      = zeros(n_edges, N_tot);
traj.T_0s     = zeros(1, N_tot);
traj.T_0r     = zeros(1, N_tot);
traj.T_F0     = zeros(1, N_tot);
traj.T_is     = zeros(n_user, N_tot);
traj.T_ir     = zeros(n_user, N_tot);
traj.d        = zeros(n_user, N_tot);
traj.c        = zeros(n_user, N_tot);
traj.theta    = zeros(n_user, N_tot);
traj.edge_idx = ei;

idx_wu = 1:N_wu;
traj.t(idx_wu)        = res_wu.t;
traj.Tout(:, idx_wu)  = res_wu.Tout;
traj.q_edges(:, idx_wu) = res_wu.q_edges;
traj.q_users(:, idx_wu) = res_wu.q_users;
traj.r_q(:, idx_wu)   = repmat(net.mdotEdges(:), 1, N_wu);
traj.T_0s(idx_wu)     = p.Tin_nom + res_wu.u;
traj.T_0r(idx_wu)     = res_wu.Tout(R0_idx, :);
traj.T_F0(idx_wu)     = res_wu.Tout(F0_idx, :);
traj.T_is(:, idx_wu)  = res_wu.T_s_i;
traj.T_ir(:, idx_wu)  = res_wu.T_r_i;
traj.d(:, idx_wu)     = res_wu.d_i;
traj.c(:, idx_wu)     = res_wu.c_i;
traj.theta(:, idx_wu) = (res_wu.T_s_i - res_wu.T_r_i) .* res_wu.q_users;

z_state = res_wu.z_final;
u_prev  = [traj.T_0s(N_wu); traj.r_q(:, N_wu); traj.T_ir(:, N_wu)];

solve_ms = zeros(1, N_cl);
exitflag = zeros(1, N_cl);

for k = 1:N_cl
    cur_idx = N_wu + k - 1;

    traj_to_now = trim_traj_local(traj, cur_idx);
    % optional controller-side measurement noise (same hook as KPC)
    if isfield(tune, 'noise_cfg') && ~isempty(tune.noise_cfg)
        nc = tune.noise_cfg;
        rs = RandStream('mt19937ar', 'Seed', nc.seed + cur_idx);
        sT = nc.sigma_T;
        sQ = nc.sigma_q_frac;
        traj_to_now.T_0s    = traj_to_now.T_0s    + sT * randn(rs, size(traj_to_now.T_0s));
        traj_to_now.T_0r    = traj_to_now.T_0r    + sT * randn(rs, size(traj_to_now.T_0r));
        traj_to_now.T_F0    = traj_to_now.T_F0    + sT * randn(rs, size(traj_to_now.T_F0));
        traj_to_now.T_is    = traj_to_now.T_is    + sT * randn(rs, size(traj_to_now.T_is));
        traj_to_now.T_ir    = traj_to_now.T_ir    + sT * randn(rs, size(traj_to_now.T_ir));
        traj_to_now.q_users = traj_to_now.q_users .* (1 + sQ * randn(rs, size(traj_to_now.q_users)));
    end
    [Z_hist, ~, ~] = lift_fn(traj_to_now, p);
    z_0 = Z_hist(:, end);

    t_now = traj.t(cur_idx);
    d_seq = zeros(n_user, Np);
    for h = 1:Np
        for i = 1:n_user
            d_seq(i, h) = compute_prosumer_Q(scenario_start + t_now + h*Ts, ...
                                             net.Nodes(con_idx(i)).params);
        end
    end

    % optional demand-forecast noise (same hook as KPC). Plant integrates
    % the true profile; only the controller sees a perturbed d_seq.
    if isfield(tune, 'demand_noise_cfg') && ~isempty(tune.demand_noise_cfg)
        dn = tune.demand_noise_cfg;
        rs_d = RandStream('mt19937ar', 'Seed', dn.seed + cur_idx);
        type = 'gaussian';
        if isfield(dn, 'type') && ~isempty(dn.type)
            type = lower(dn.type);
        end
        switch type
            case 'gaussian'
                d_seq = d_seq .* (1 + dn.sigma_frac * randn(rs_d, size(d_seq)));
            case 'bias'
                d_seq = d_seq * (1 + dn.bias_frac);
            case 'step'
                if cur_idx >= dn.step_idx
                    d_seq = d_seq * (1 + dn.step_frac);
                end
            case 'ar1'
                if ~exist('ar1_state', 'var') || isempty(ar1_state)
                    ar1_state = zeros(size(d_seq));
                end
                innov = dn.sigma_frac * sqrt(1 - dn.phi^2) * randn(rs_d, size(d_seq));
                ar1_state = dn.phi * ar1_state + innov;
                d_seq = d_seq .* (1 + ar1_state);
            otherwise
                error('lmpc_step_loop: unknown demand_noise type "%s"', type);
        end
        d_seq = max(0, d_seq);
    end

    % warm-start: shift-by-one of last U_seq, same as KPC
    if isfield(tune, 'use_warm_start') && tune.use_warm_start && k > 1
        prev_U = info.U_seq;
        shifted = [prev_U(:, 2:end), prev_U(:, end)];
        tune.warm_start_U = shifted(:);
    end

    % rebuild the iterated predictor at the current operating point, with
    % the known demand forecast folded in (the one place this differs from
    % KPC, which uses a fixed pred and folds demand via G_d at solve time)
    pred = build_lmpc_pred(fit_iter, Np, z_0, u_prev, p, n_edges, n_user, d_seq);

    [u_opt, info] = lmpc_solve(z_0, d_seq, u_prev, pred, p, net, tune);
    solve_ms(k) = info.solve_ms;
    exitflag(k) = info.exitflag;

    % apply only the first input for one Ts (receding horizon) - same as KPC
    T_0s_apply = u_opt(1);
    r_q_apply  = u_opt(2:1+n_edges);
    p_step = p;
    p_step.r_q_fun = @(t) r_q_apply;
    p_step.t_offset = scenario_start + t_now;
    u_fun = @(t) T_0s_apply - p.Tin_nom;
    w_fun = @(t) 1.0;
    net_step = net;
    net_step.q0 = traj.q_edges(:, cur_idx);
    res_step = simulate_plant(net_step, z_state, p_step, u_fun, w_fun, Ts);

    new_idx = cur_idx + 1;
    traj.t(new_idx)        = t_now + Ts;
    traj.Tout(:, new_idx)  = res_step.Tout(:, end);
    traj.q_edges(:, new_idx) = res_step.q_edges(:, end);
    traj.q_users(:, new_idx) = res_step.q_users(:, end);
    traj.r_q(:, new_idx)   = res_step.r_q(:, end);
    traj.T_0s(new_idx)     = p.Tin_nom + res_step.u(end);
    traj.T_0r(new_idx)     = res_step.Tout(R0_idx, end);
    traj.T_F0(new_idx)     = res_step.Tout(F0_idx, end);
    traj.T_is(:, new_idx)  = res_step.T_s_i(:, end);
    traj.T_ir(:, new_idx)  = res_step.T_r_i(:, end);
    traj.d(:, new_idx)     = res_step.d_i(:, end);
    traj.c(:, new_idx)     = res_step.c_i(:, end);
    traj.theta(:, new_idx) = (res_step.T_s_i(:, end) - res_step.T_r_i(:, end)) .* res_step.q_users(:, end);

    z_state = res_step.z_final;
    u_prev  = u_opt;
end

cl_idx = N_wu + (1:N_cl);
res.d_i      = traj.d(:, cl_idx);
res.c_i      = traj.c(:, cl_idx);
res.r_q      = traj.r_q(:, cl_idx);
res.T_0s     = traj.T_0s(cl_idx);
res.T_0r     = traj.T_0r(cl_idx);
res.T_F0     = traj.T_F0(cl_idx);
res.T_is     = traj.T_is(:, cl_idx);
res.T_ir     = traj.T_ir(:, cl_idx);
res.theta    = traj.theta(:, cl_idx);
res.q_users  = traj.q_users(:, cl_idx);
res.Q_net    = sum(res.q_users, 1);
res.solve_ms = solve_ms;
res.exitflag = exitflag;
end


function tt = trim_traj_local(traj, n)
flds = fieldnames(traj);
tt = struct();
for f = 1:numel(flds)
    nm = flds{f};
    v = traj.(nm);
    if ischar(v) || isstruct(v) || isscalar(v)
        tt.(nm) = v;
        continue;
    end
    if isrow(v)
        tt.(nm) = v(1:n);
    elseif ismatrix(v) && size(v, 2) >= n
        tt.(nm) = v(:, 1:n);
    else
        tt.(nm) = v;
    end
end
end
