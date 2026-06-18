function sig = excitation(type, t, opts)
% EXCITATION  Generate excitation signals for data collection.
%   sig = excitation('prbs', t, struct('amp', 6, 'dwell', 900, 'seed', 1))
%   sig = excitation('multisine', t, struct('amp', 5, 'n_harm', 5, 'seed', 1))
%   sig = excitation('steps', t, struct('levels', [-6 -3 0 3 6], 'dwell', 3600))
%   sig = excitation('ramp', t, struct('lo', 0.5, 'hi', 1.5))
%   sig = excitation('diurnal', t, struct('offset', 1.0, 'amp', 0.3, 'period', 8*3600))
%   sig = excitation('constant', t, struct('value', 1.0))
%
%   Inputs:
%     type - signal type string
%     t    - time vector [s]
%     opts - struct with signal-specific options
%
%   Output:
%     sig  - signal vector, same length as t

if nargin < 3, opts = struct(); end

switch lower(type)
    case 'prbs',      sig = gen_prbs(t, opts);
    case 'multisine', sig = gen_multisine(t, opts);
    case 'steps',     sig = gen_steps(t, opts);
    case 'ramp',      sig = gen_ramp(t, opts);
    case 'diurnal',   sig = gen_diurnal(t, opts);
    case 'constant',  sig = gen_constant(t, opts);
    otherwise,        error('excitation: unknown type "%s"', type);
end
end

%% ---- PRBS ----
function sig = gen_prbs(t, o)
amp   = getf(o, 'amp', 6);
dwell = getf(o, 'dwell', 900);
seed  = getf(o, 'seed', 1);
lo    = getf(o, 'lo', []);
hi    = getf(o, 'hi', []);

rng(seed, 'twister');
n_sw = ceil((t(end) - t(1)) / dwell) + 1;
bits = 2 * (randi(2, 1, n_sw) - 1.5);  % +/-1

sig = zeros(size(t));
for k = 1:numel(t)
    idx = min(floor((t(k) - t(1)) / dwell) + 1, n_sw);
    sig(k) = bits(idx);
end

if ~isempty(lo) && ~isempty(hi)
    sig = lo + (sig + 1)/2 * (hi - lo);
else
    sig = sig * amp;
end
end

%% ---- MULTISINE ----
function sig = gen_multisine(t, o)
amp    = getf(o, 'amp', 5);
n_harm = getf(o, 'n_harm', 5);
seed   = getf(o, 'seed', 1);

rng(seed, 'twister');
T_tot  = t(end) - t(1);
phases = 2*pi * rand(1, n_harm);

sig = zeros(size(t));
for h = 1:n_harm
    sig = sig + sin(2*pi * h / T_tot * t + phases(h));
end
sig = sig / max(abs(sig)) * amp;
end

%% ---- STEPS ----
function sig = gen_steps(t, o)
levels = getf(o, 'levels', [-6 -3 0 3 6]);
dwell  = getf(o, 'dwell', 3600);

n_rep = ceil((t(end) - t(1)) / (dwell * numel(levels)));
seq = repmat(levels(:)', 1, n_rep);

sig = zeros(size(t));
for k = 1:numel(t)
    idx = min(floor((t(k) - t(1)) / dwell) + 1, numel(seq));
    sig(k) = seq(idx);
end
end

%% ---- RAMP (triangle) ----
function sig = gen_ramp(t, o)
lo = getf(o, 'lo', 0.5);
hi = getf(o, 'hi', 1.5);
T_half = (t(end) - t(1)) / 2;

sig = zeros(size(t));
for k = 1:numel(t)
    tau = t(k) - t(1);
    if tau <= T_half
        sig(k) = lo + (hi - lo) * tau / T_half;
    else
        sig(k) = hi - (hi - lo) * (tau - T_half) / T_half;
    end
end
end

%% ---- DIURNAL ----
function sig = gen_diurnal(t, o)
offset = getf(o, 'offset', 1.0);
amp    = getf(o, 'amp', 0.3);
period = getf(o, 'period', 8*3600);
sig = offset + amp * sin(2*pi * t / period);
end

%% ---- CONSTANT ----
function sig = gen_constant(t, o)
sig = getf(o, 'value', 0) * ones(size(t));
end

%% ---- Helper ----
function v = getf(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
