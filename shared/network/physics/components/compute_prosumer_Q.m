function Q = compute_prosumer_Q(t, prm)
% COMPUTE_PROSUMER_Q   Heat demand [W] for one prosumer at time t.
%
% Supported profile types (prm.profile):
%   'flat'         - constant Q_peak
%   'diurnal'      - legacy single sin with peak at 15:00
%   'residential'  - bimodal (morning + evening peaks), IEA Annex 60 style
%   'commercial'   - unimodal with business-hours window (~13:00 peak)
%   'shop'         - unimodal retail window (~17:00 peak, 10:00-20:00 active)
%
% Shape references: IEA Annex 60 neighbourhood benchmark; ORNL commercial
% and residential hourly load profiles.

Q_peak = prm.Q_peak;

if ~isfield(prm,'profile') || isempty(prm.profile)
    Q = Q_peak;
    return;
end

% optional time offset for closed-loop runners that step the plant in
% relative-time chunks but want the demand profile to track absolute
% wall-clock time. simulate_plant injects prm.t_offset before each call.
if isfield(prm, 't_offset') && ~isempty(prm.t_offset)
    t = t + prm.t_offset;
end

t_day = mod(t/3600, 24);   % hour of day in [0, 24)

switch lower(prm.profile)
    case 'flat'
        Q = Q_peak;

    case 'diurnal'
        % Legacy profile: sin with phase shift so peak lands at 15:00.
        Q = Q_peak * (0.25 + 0.75 * 0.5 * (1 + sin(2*pi*(t/86400 - 0.375))));

    case 'residential'
        % Domestic pattern: morning (~7:30) and evening (~19:00) peaks,
        % with a baseload around 20 percent (always-on hot water, standing losses).
        g_m = exp(-0.5 * ((t_day - 7.5) / 1.5)^2);
        g_e = exp(-0.5 * ((t_day - 19.0) / 2.0)^2);
        Q = Q_peak * (0.20 + 0.40 * g_m + 0.40 * g_e);

    case 'commercial'
        % Office/commercial pattern: single workhours peak (~13:00), gated by
        % sigmoid open/close at 7:00 and 18:00, baseload 10 percent at night.
        g       = exp(-0.5 * ((t_day - 13.0) / 3.0)^2);
        s_open  = 1 / (1 + exp(-2 * (t_day - 7.0)));
        s_close = 1 / (1 + exp(-2 * (18.0 - t_day)));
        Q = Q_peak * (0.10 + 0.90 * g * s_open * s_close);

    case 'shop'
        % Retail pattern: afternoon peak (~17:00), active 10:00-20:00,
        % very small baseload outside opening hours.
        g       = exp(-0.5 * ((t_day - 17.0) / 2.5)^2);
        s_open  = 1 / (1 + exp(-2 * (t_day - 10.0)));
        s_close = 1 / (1 + exp(-2 * (20.0 - t_day)));
        Q = Q_peak * (0.05 + 0.95 * g * s_open * s_close);

    otherwise
        error('compute_prosumer_Q: unknown profile type "%s"', prm.profile);
end
end
