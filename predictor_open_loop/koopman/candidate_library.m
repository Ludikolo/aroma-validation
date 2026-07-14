function [Z, names, meta] = candidate_library(traj, p)
% CANDIDATE_LIBRARY  Lifting psi(xi) for the Koopman predictor.
%
% Three feature blocks plus a constant:
%   A thermal    T_F0, T_s^i (5), T_r^i (5), T_0r                 12
%   B hydraulic  q^i (5), Q_net                                    6
%   C theta      theta^i = (T_s^i - T_r^i) q^i, per consumer       5
% Optional source-delay block adds T_F0(t-d) and delayed bilinears.
%
% Including theta^i in the lift is what makes c^i = cp theta^i
% linear in z, which is the trick that lets the QP run.
%
% The T_r^i features are the measured returns at time k (plant states, the
% substations set them); the same channel also appears in the input stack U
% as the planned T_ir sequence, which is never applied to the plant.
%
% THEORY index:
%   line 61: theta
%   line 73: delay embedding
%   line 102: delayed bilinear
%   line 183: exergy block

n_user = size(traj.q_users, 1);
N      = numel(traj.t);
% The per-consumer loops are parametric in n_user (5 in the production
% fits, more on the synthetic plants used for the scaling experiments).
feats = {};
nm    = {};
idx   = struct();

% Block A: thermal (12). State-only lift: the source-supply command T_0s is an
% input, so it enters the predictor through the B[u;d] regressor, not the lift.
feats{end+1} = traj.T_F0(:)';            nm{end+1} = 'T_F0';
idx.T_F0 = numel(feats);

idx.T_is = zeros(1, n_user);
for i = 1:n_user
    feats{end+1} = traj.T_is(i, :);      nm{end+1} = sprintf('T_s_C%d', i);
    idx.T_is(i) = numel(feats);
end

idx.T_ir = zeros(1, n_user);
for i = 1:n_user
    feats{end+1} = traj.T_ir(i, :);      nm{end+1} = sprintf('T_r_C%d', i);
    idx.T_ir(i) = numel(feats);
end

feats{end+1} = traj.T_0r(:)';            nm{end+1} = 'T_0r';
idx.T_0r = numel(feats);

% Block B: hydraulic (6)
idx.q = zeros(1, n_user);
for i = 1:n_user
    feats{end+1} = traj.q_users(i, :);   nm{end+1} = sprintf('q_C%d', i);
    idx.q(i) = numel(feats);
end

feats{end+1} = sum(traj.q_users, 1);     nm{end+1} = 'Q_net';
idx.Q_net = numel(feats);

% THEORY (theta): theta = (Ts - Tr)*q per consumer; c = cp*theta is a linear read-out of z, key for the convex QP
% Block C: theta (5)
idx.theta = zeros(1, n_user);
for i = 1:n_user
    feats{end+1} = traj.theta(i, :);     nm{end+1} = sprintf('theta_C%d', i);
    idx.theta(i) = numel(feats);
end

% constant
feats{end+1} = ones(1, N);               nm{end+1} = 'one';
idx.one = numel(feats);

% THEORY (delay embedding): transport delay source to consumer is 17-50 samples, long horizons need past T_F0 (Takens)
% optional delay block (Takens-style upstream embedding)
% Carries the source-delay-embedding insight. The transport
% delay F_0 -> consumer is 17-50 samples depending on flow, so without
% past T_F0 the linear lift cannot reproduce theta^i at long horizons.
max_lag = 0;
if isfield(p, 'o2') && isfield(p.o2, 'use_delays') && p.o2.use_delays
    lags = p.o2.delay_lags;
    max_lag = max(max_lag, max(lags));

    idx.T_F0_lag = zeros(1, numel(lags));
    for L = 1:numel(lags)
        d = lags(L);
        feats{end+1} = lag_signal(traj.T_F0(:)', d);
        nm{end+1}    = sprintf('T_F0_lag%d', d);
        idx.T_F0_lag(L) = numel(feats);
    end

    if isfield(p.o2, 'use_flow_disc') && p.o2.use_flow_disc
        Q_net = sum(traj.q_users, 1);
        idx.T_F0_lag_x_Qnet = zeros(1, numel(lags));
        for L = 1:numel(lags)
            d = lags(L);
            feats{end+1} = lag_signal(traj.T_F0(:)', d) .* Q_net;
            nm{end+1}    = sprintf('T_F0_lag%d_Qnet', d);
            idx.T_F0_lag_x_Qnet(L) = numel(feats);
        end
    end

    % THEORY (delayed bilinear): T_F0(t-l)*q couples what was injected upstream with how fast it travels
    % delayed bilinear T_F0(t-d) * q^i(t): what was injected at the
    % source d steps ago, modulated by current user-side flow.
    if isfield(p.o2, 'use_delayed_bilinear') && p.o2.use_delayed_bilinear
        idx.T_F0_lag_x_q = zeros(n_user, numel(lags));
        for i = 1:n_user
            for L = 1:numel(lags)
                d = lags(L);
                feats{end+1} = lag_signal(traj.T_F0(:)', d) .* traj.q_users(i, :);
                nm{end+1}    = sprintf('T_F0_lag%d_q_C%d', d, i);
                idx.T_F0_lag_x_q(i, L) = numel(feats);
            end
        end
    end
end

% optional substation-aware extras (gated)
% These break through the linear-lift ceiling on the substation
% clipping nonlinearity c = min(d, c_p q (T_s - T_r_min)). Without
% them mean test R^2 at h=1 caps at 0.74 across consumers; with them
% the lift can represent both regimes (saturated vs demand-limited).
% Literature: Peitz & Klus 2019, Mauroy/Susuki/Mezic 2020 (saturation
% observables); Mauroy & Mezic 2016, Kaiser/Kutz/Brunton 2018 (diurnal
% trig basis for periodically-forced systems). Gated by
% p.o2.use_substation_extras = true so legacy fits keep loading.
use_extras = isfield(p, 'o2') && isfield(p.o2, 'use_substation_extras') ...
             && p.o2.use_substation_extras;
if use_extras
    cp_w  = p.cp;
    Trmin = p.consumer.T_r_min;

    % saturation observables: (T_s^(i) - T_r_min)_+ * q^(i) * c_p  (W)
    capacity = cp_w * traj.q_users .* (traj.T_is - Trmin);
    sat_W    = max(0, capacity);
    idx.sat = zeros(1, n_user);
    for i = 1:n_user
        feats{end+1} = sat_W(i, :);
        nm{end+1}    = sprintf('sat_W_C%d', i);
        idx.sat(i)   = numel(feats);
    end

    % smooth regime indicator s_i = sigmoid((d - capacity)/scale)
    scale = 5000;     % W
    regime = 1 ./ (1 + exp(-(traj.d - capacity) / scale));
    idx.regime = zeros(1, n_user);
    for i = 1:n_user
        feats{end+1} = regime(i, :);
        nm{end+1}    = sprintf('regime_C%d', i);
        idx.regime(i) = numel(feats);
    end

    % regime cross-products
    reg_x_d   = regime .* traj.d;            % anticipated c when saturating
    nreg_x_d  = (1 - regime) .* traj.d;       % anticipated c when demand-limited
    reg_x_cap = regime .* capacity;
    for i = 1:n_user
        feats{end+1} = reg_x_d(i, :);   nm{end+1} = sprintf('reg_d_C%d', i);
        feats{end+1} = nreg_x_d(i, :);  nm{end+1} = sprintf('nreg_d_C%d', i);
        feats{end+1} = reg_x_cap(i, :); nm{end+1} = sprintf('reg_cap_C%d', i);
    end

    % diurnal trig basis (Kaiser/Kutz/Brunton 2018)
    t_abs = traj.t;
    if isfield(traj, 't_offset'), t_abs = t_abs + traj.t_offset; end
    omega = 2 * pi / 86400;
    sin_t = sin(omega * t_abs);
    cos_t = cos(omega * t_abs);
    feats{end+1} = sin_t;  nm{end+1} = 'sin_diurnal';
    feats{end+1} = cos_t;  nm{end+1} = 'cos_diurnal';
    for i = 1:n_user
        feats{end+1} = sin_t .* traj.q_users(i, :); nm{end+1} = sprintf('sin_q_C%d', i);
        feats{end+1} = cos_t .* traj.q_users(i, :); nm{end+1} = sprintf('cos_q_C%d', i);
    end
end

% optional exergy block (gated for exergy production)
% 11 bilinears appended AFTER the base 36-feature lift, filling rows 37..47
% in this fixed order (Tamb = p.Text):
%   q^i*(T_s^i - Tamb) (5), q^i*(T_r^i - Tamb) (5), Q_net*(T_0r - Tamb) (1).
% OFF by default -> base 36-feature lift. Enabled only when exergy is the
% adopted production lift.
% THEORY (exergy block): bilinears q*(T - Tamb) on top of the base 36, gives the production lift n_z = 47
use_exergy = isfield(p, 'o2') && isfield(p.o2, 'use_exergy') && p.o2.use_exergy;
if use_exergy
    T_amb = p.Text;
    Qnet  = sum(traj.q_users, 1);
    T0r   = traj.T_0r(:)';
    idx.exergy_s = zeros(1, n_user);
    for i = 1:n_user
        feats{end+1} = traj.q_users(i, :) .* (traj.T_is(i, :) - T_amb);
        nm{end+1}    = sprintf('exergy_s_C%d', i);
        idx.exergy_s(i) = numel(feats);
    end
    idx.exergy_r = zeros(1, n_user);
    for i = 1:n_user
        feats{end+1} = traj.q_users(i, :) .* (traj.T_ir(i, :) - T_amb);
        nm{end+1}    = sprintf('exergy_r_C%d', i);
        idx.exergy_r(i) = numel(feats);
    end
    feats{end+1} = Qnet .* (T0r - T_amb);  nm{end+1} = 'exergy_net_r';
    idx.exergy_net_r = numel(feats);
end

% pack into matrix
Z = vertcat(feats{:});                   % (n_z x N)

names = nm;
meta.n_z       = size(Z, 1);
meta.idx       = idx;
meta.n_user    = n_user;
meta.N         = N;
meta.max_lag   = max_lag;
meta.valid_start = max_lag + 1;          % skip warm-up padded samples in fits

end


function s = lag_signal(s, d)
% shift s by d steps; pad the head with the first valid value
N = numel(s);
if d <= 0
    return;
end
s = [repmat(s(1), 1, d), s(1:max(0, N - d))];
end
