function dT = pipe_rhs_bidir(~, T, p, Tin)
%PIPE_RHS_BIDIR  1D pipe semi-discrete RHS: upwind advection + ambient loss.
%   dT = pipe_rhs_bidir(t, T, p, Tin) returns dT/dt for a pipe discretized
%   into p.n cells of length p.dx. A constant velocity p.v advects heat
%   along the pipe. A linear loss term -alpha*(T - Text) models ambient
%   heat exchange (paper-style ambient shift).
%
% Inputs
%   T    : [n x 1] cell temperatures (absolute, °C)
%   p    : struct with fields
%            .n     number of cells (integer >=1)
%            .dx    cell length [m] (>0)
%            .v     fluid velocity [m/s] (sign gives flow direction)
%            .alpha ambient loss coeff [1/s]
%            .Text  ambient/external temperature [°C]
%   Tin  : scalar inlet temperature [°C] at the UPWIND boundary
%
% Notes
% - Upwind first-order advection:
%     v >= 0: inlet at cell 1 (ghost = Tin), flux uses T(i)-T(i-1)
%     v <  0: inlet at cell n (ghost = Tin), flux uses T(i)-T(i+1)
% - Loss term is applied as -alpha*(T - Text) per cell.
% - This routine is deliberately simple and close to the paper's semi-discrete form.
%
% Example
%   p = struct('n',10,'dx',1,'v',1,'alpha',0.01,'Text',20);
%   T = 20*ones(p.n,1);
%   dT = pipe_rhs_bidir(0, T, p, 10);

    T = T(:);
    n = numel(T);

    v  = p.v;
    dx = max(p.dx, eps);                 % avoid division by zero
    a  = p.alpha;
    Ta = p.Text;

    dT = zeros(n,1);

    if v >= 0
        % Upwind from the left (inlet at cell 1)
        % Cell 1 uses ghost value Tin; others use previous cell.
        adv1      = -(v/dx) * (T(1)   - Tin);
        loss1     = -a * (T(1) - Ta);
        dT(1)     = adv1 + loss1;

        if n > 1
            adv_i     = -(v/dx) * (T(2:n) - T(1:n-1));
            loss_i    = -a * (T(2:n) - Ta);
            dT(2:n)   = adv_i + loss_i;
        end
    else
        % Upwind from the right (inlet at cell n)
        % Cell n uses ghost value Tin; others use next cell.
        advn      = -(v/dx) * (T(n)   - Tin);
        lossn     = -a * (T(n) - Ta);
        dT(n)     = advn + lossn;

        if n > 1
            adv_i     = -(v/dx) * (T(1:n-1) - T(2:n));
            loss_i    = -a * (T(1:n-1) - Ta);
            dT(1:n-1) = adv_i + loss_i;
        end
    end
end


