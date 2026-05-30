function [dTtank, Tout] = tes_cold(t, Ttank, u, p)
%TES_COLD  Well-mixed cold TES energy balance.
%   [dTtank, Tout] = tes_cold(t, Ttank, u, p)
%
%   Inputs:
%     t     - time [s] (unused, kept for ODE interface)
%     Ttank - tank temperature [°C]
%     u     - struct: u.mdot [kg/s], u.Tin [°C]
%     p     - struct: p.C [J/K], p.UA [W/K], p.cp [J/(kg*K)], p.Text [°C]
%
%   Dynamics:  dTtank = (mdot*cp*(Tin - Ttank) - UA*(Ttank - Text)) / C
%   Outlet temperature equals tank temperature (well-mixed assumption).
%

% Unpack inputs
mdot = u.mdot;
Tin  = u.Tin;

C  = p.C;
UA = p.UA;
cp = p.cp;
Text = p.Text;

% Energy balance on the mixed tank
dTtank = (mdot * cp * (Tin - Ttank) - UA * (Ttank - Text)) / C;

% Mixed outlet temperature equals tank temperature
Tout = Ttank;
end