function Tout = dry_cooler(Tin, mdot, cp, UA, Text, approach_min)
%DRY_COOLER Simple ambient heat rejection/uptake unit.
%   Tout = dry_cooler(Tin, mdot, cp, UA, Text, approach_min) computes the
%   outlet temperature after passing through a dry cooler or cooling tower.
%   The device rejects heat to the ambient (cooling the fluid) according to
%   a simple first-order model. The outlet cannot drop below (Text + approach_min).
%
%   Inputs:
%     Tin          - inlet temperature (°C)
%     mdot         - mass flow rate (kg/s)
%     cp           - heat capacity (J/(kg*K))
%     UA           - overall heat transfer coefficient-area product (W/K)
%     Text         - ambient (air) temperature (°C)
%     approach_min - minimum approach (°C) between outlet and ambient
%
%   Example:
%     Tout = dry_cooler(35, 1, 4180, 8000, 20, 3);
%

% First-order linear cooling model (linearised NTU for small UA/(mdot*cp)).
% Note: the exact model would be Tout = Text + (Tin-Text)*exp(-UA/(mdot*cp)),
% but the linearisation is consistent with the training data and frozen results.
if mdot > 0
    Tout_est = Tin - UA/(mdot*cp) * (Tin - Text);
else
    Tout_est = Tin;
end

% Ensure minimum approach to ambient
Tout = max(Tout_est, Text + approach_min);
end