function [Tout, c_i, d_i] = prosumer_house(t, Tsup, mdot, cp, prm, mode)
% PROSUMER_HOUSE  Substation: heat extraction with T_r >= T_r_min and
% c <= d. Falls back to legacy c = d when prm.T_r_min is not set.
% Energy balance: c = mdot * cp * (Tsup - Tout) in heating.
% 
%
% THEORY index:
%   line 27: substation saturation



% pick local mode (per-consumer override is allowed)
if isfield(prm,'mode') && ~isempty(prm.mode)
    local_mode = lower(prm.mode);
else
    local_mode = lower(mode);
end

d_i = compute_prosumer_Q(t, prm);

% when prm.T_r_min is set, the substation clips c_i; otherwise it lets c_i = d_i pass through
has_floor = isfield(prm, 'T_r_min') && ~isempty(prm.T_r_min);

switch local_mode
    case 'heating'
        if has_floor
            % THEORY (substation saturation): c = min(d, mdot*cp*(Tsup - Tr_min)); low flow caps delivery, the hard regime
            c_max = max(mdot * cp * (Tsup - prm.T_r_min), 0);
            c_i   = min(max(d_i, 0), c_max);
        else
            c_i = d_i;
        end
        Tout = Tsup - c_i / (mdot * cp);

    case 'cooling'
        % cooling mode does not use T_r_min, so c_i = d_i passes through
        c_i  = d_i;
        Tout = Tsup + c_i / (mdot * cp);

    otherwise
        error('prosumer_house: unknown mode "%s". Use cooling/heating.', local_mode);
end
end
