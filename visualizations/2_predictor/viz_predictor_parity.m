% VIZ_PREDICTOR_PARITY  Predicted vs actual delivered heat, on the diagonal.
%
% Four panels for a normal day, every consumer: the forecast at 1, 5, 10 and 16
% steps ahead (15 min .. 4 h). The 1-step is the move the controller applies; the
% other three are genuine multi-step forecasts. On the diagonal = accurate, and
% it stays there at every horizon. The demand forecast over the horizon is part
% of the predictor input (as in the paper).

clear; clc;
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(root); startup();

% all plant/predictor constants live in params(); p.cp = 4180 J/(kg.K) does the
% theta -> heat conversion below, p.Ts = 900 s is the 15-min sample (so step h = h*15 min)
p = params(); p.Ts = 900;
% the trained V_seq predictor: one (cp, dp) pair per (horizon h, consumer i), fit
% offline by fit_vseq. cp weights the lift z_k, dp weights the input/demand sequence
F = load(fullfile(root, 'predictor_open_loop', 'results', 'vseq_fits_full.mat')); fits = F.fits;
% the fit carries the demand forecast as a regressor block, so dp has a demand tail
% (demand-in-V). flag it so the rollout below adds that block
has_dV = isfield(fits, 'includes_d_in_V') && fits.includes_d_in_V;

% build the network and drive it as a normal operational day: hold the edge flows at
% their nominal values mdotEdges with the built-in diurnal demand. this is the
% deployment regime we forecast on, not the PRBS data the predictor was trained on
[net, z0] = build_plant(p); mdotE = net.mdotEdges(:);
R0 = find(strcmp({net.Nodes.name}, 'R0')); F0 = find(strcmp({net.Nodes.name}, 'F0'));
% nominal flows for every edge, no source-temp offset (handled below), 26 h of sim
pr = p; pr.t_offset = 0; pr.r_q_fun = @(t) mdotE;
res = simulate_plant(net, z0, pr, @(t) 0, @(t) 1.0, 26 * 3600);
tr.t = res.t; tr.t_offset = 0; tr.T_0s = p.Tin_nom + res.u; tr.r_q = res.r_q; tr.T_ir = res.T_r_i;
% theta^i = (T_s^i - T_r^i) * q^i [K.kg/s], the per-consumer delivered-heat feature;
% multiply by p.cp later to get heat in W
tr.theta = (res.T_s_i - res.T_r_i) .* res.q_users; tr.T_is = res.T_s_i; tr.T_0r = res.Tout(R0, :);
tr.T_F0 = res.Tout(F0, :); tr.q_users = res.q_users; tr.q_edges = res.q_edges; tr.Tout = res.Tout; tr.d = res.d_i;
% Z = the 49-feature Koopman lift (37 base + 12 exergy bilinears), one column per step.
% U = the control input stacked as [T_0s; r_q; T_ir]: 1 source supply temperature,
% 29 edge flow setpoints, 5 consumer return-temp setpoints, so n_u = 35.
% D = per-consumer heat demand [W], the known forecast that feeds the dp demand tail
[Z, ~, meta] = candidate_library(tr, p); U = [tr.T_0s(:)'; tr.r_q; tr.T_ir]; D = tr.d;
% n_user = 5 (meta.idx.theta points at the 5 theta features); n_u = 35 from the stack above
n_user = numel(meta.idx.theta); n_u = size(U, 1);
r2 = @(a, b) 1 - sum((a - b).^2) / max(sum((b - mean(b)).^2), eps);
col = lines(n_user); lims = [0, 26];

% one panel per horizon: 1 step (= the controller's applied move) plus 5, 10, 16
% steps = 15 min, 75 min, 150 min, 240 min ahead. 16 is the top of the trained sweep
horizons = [1 5 10 16];
figure('Color', 'w', 'Position', [100 100 900 840]);
for pidx = 1:numel(horizons)
    h = horizons(pidx);
    hi = find(fits.horizons == h);     % row of fits.cp/dp that holds this horizon's fit
    % skip the lift warm-up window, and stop h steps early so target k+h stays in range
    N = size(Z, 2); kk = max(1, meta.valid_start) : N - h;
    subplot(2, 2, pidx); hold on;
    for i = 1:n_user
        % actual delivered heat: p.cp [J/(kg.K)] turns theta [K.kg/s] into W
        act = p.cp * tr.theta(i, kk + h); prd = zeros(size(kk));
        for n = 1:numel(kk)
            % direct h-step forecast of theta: cp on the lift at k, plus the control part
            % of dp on the h applied inputs u_k..u_{k+h-1}. dp's first n_u*h entries are
            % the control coefficients, the rest are the demand tail
            k = kk(n); base = fits.cp{hi,i}(:)' * Z(:,k) + fits.dp{hi,i}(1:n_u*h)' * reshape(U(:,k:k+h-1),[],1);
            % demand-in-V: add the demand tail acting on the forecast d_{k+1}..d_{k+h}
            if has_dV, base = base + fits.dp{hi,i}(n_u*h+1:end)' * reshape(D(:,k+1:k+h),[],1); end
            prd(n) = p.cp * base;   % theta -> heat [W], same conversion as act
        end
        plot(act/1e3, prd/1e3, '.', 'Color', col(i,:), 'DisplayName', sprintf('C%d (R^2 %.2f)', i, r2(prd, act)));
    end
    plot(lims, lims, 'k--', 'HandleVisibility', 'off');
    xlim(lims); ylim(lims); axis square; grid on;
    xlabel('actual delivered heat [kW]'); ylabel('predicted delivered heat [kW]');
    legend('Location', 'southeast');
    if h * 15 < 120, tlab = sprintf('%d min', h * 15); else, tlab = sprintf('%g h', h * 15 / 60); end
    title(sprintf('%d-step forecast (%s ahead)', h, tlab));
end
sgtitle('Delivered-heat forecast parity, 1 to 16 steps ahead (15 min to 4 h)');
exportgraphics(gcf, fullfile(here, 'viz_predictor_parity.pdf'), 'ContentType', 'vector');
fprintf('Saved viz_predictor_parity.pdf\n');
