%==========================================================================
% FileName: ex05_generalization
% Author: Yi CAI
% Description:
%   Pure angular approximation study for the hybrid velocity representation
%   on S^1 (2DT topology). Compares five fixed representations across seven
%   test cases that arise as angular profiles of free-transport solutions
%   at a fixed observation point (x0, y0) at time t0:
%
%       f(theta) = g(x0 - cos(theta)*t0, y0 - sin(theta)*t0, theta)
%
%   Test cases:
%       1. Smooth transport pattern
%       2-4. Gaussian beams (sigma = 0.30, 0.15, 0.08)
%       5-7. Delta beams (epsilon = 0.20, 0.10, 0.05)
%
%   Representations (MdeMat element API):
%       1. Pure modal  (L2SphereElement.modal)
%       2-4. Hybrid (3/4-1/4, 1/2-1/2, 1/4-3/4) (SumSpace)
%       5. Pure nodal  (L2SphereElement.nodal)
%
%   Metrics: relative L2 error, relative L-inf error, TV penalty.
%==========================================================================

clear; clc; close all;

%% ============================================================
%  Reference grid on [0, 2pi)
%% ============================================================

Nt_ref = 4096;
theta  = linspace(0, 2*pi, Nt_ref + 1);
theta(end) = [];
dtheta = 2*pi / Nt_ref;

% Fixed parameters in the test functions
x0 = 0.25;
y0 = -0.15;
t0 = 1.0;

%% ============================================================
%  Total velocity DOFs
%% ============================================================

D_list = [9, 17, 25, 33, 49, 65];
nD     = numel(D_list);

%% ============================================================
%  Fixed representation configurations
%% ============================================================

config_names = { ...
    'pure_modal', ...
    'hybrid_3_4_1_4', ...
    'hybrid_1_2_1_2', ...
    'hybrid_1_4_3_4', ...
    'pure_nodal'};

config_labels = { ...
    'Pure modal', ...
    'Hybrid (3/4, 1/4)', ...
    'Hybrid (1/2, 1/2)', ...
    'Hybrid (1/4, 3/4)', ...
    'Pure nodal'};

nConfigs = numel(config_names);

%% ============================================================
%  Representative case library
%% ============================================================

case_names = { ...
    'smooth_transport', ...
    'gaussian_beam_sigma_030', ...
    'gaussian_beam_sigma_015', ...
    'gaussian_beam_sigma_008', ...
    'delta_beam_eps_020', ...
    'delta_beam_eps_010', ...
    'delta_beam_eps_005'};

case_labels = { ...
    'Smooth transport', ...
    'Gaussian beam (\sigma=0.30)', ...
    'Gaussian beam (\sigma=0.15)', ...
    'Gaussian beam (\sigma=0.08)', ...
    'Delta beam (\epsilon=0.20)', ...
    'Delta beam (\epsilon=0.10)', ...
    'Delta beam (\epsilon=0.05)'};

nCases = numel(case_names);

%% ============================================================
%  Storage
%% ============================================================

errL2   = zeros(nCases, nConfigs, nD);
errLinf = zeros(nCases, nConfigs, nD);
penTV   = zeros(nCases, nConfigs, nD);

Nm_store = zeros(nConfigs, nD);
Mr_store = zeros(nConfigs, nD);

%% ============================================================
%  Precompute configuration parameters
%% ============================================================

fprintf('\n============================================================\n');
fprintf('Fixed representation configurations\n');
fprintf('============================================================\n');

for iconfig = 1:nConfigs
    fprintf('\n%s\n', config_labels{iconfig});
    for id = 1:nD
        D = D_list(id);
        [Nm_store(iconfig, id), Mr_store(iconfig, id)] = ...
            choose_representation(D, config_names{iconfig});

        modal_dof = 2 * Nm_store(iconfig, id) + 1;
        nodal_dof = Mr_store(iconfig, id);

        fprintf('  D = %2d -> Nm = %2d, modal DOF = %2d, Mr = %2d\n', ...
            D, Nm_store(iconfig, id), modal_dof, nodal_dof);
    end
end

%% ============================================================
%  Main experiment loop
%% ============================================================

for icase = 1:nCases

    case_name = case_names{icase};
    f_func    = @(th) target_function(th, case_name, x0, y0, t0);
    f         = f_func(theta);
    normL2    = sqrt(sum(abs(f).^2) * dtheta);

    fprintf('\n============================================================\n');
    fprintf('Case %d / %d: %s\n', icase, nCases, case_name);
    fprintf('============================================================\n');

    for iconfig = 1:nConfigs
        fprintf('%s\n', config_labels{iconfig});

        for id = 1:nD
            D  = D_list(id);
            Nm = Nm_store(iconfig, id);
            Mr = Mr_store(iconfig, id);

            fh = approximate_with_mdemat(f_func, theta, ...
                config_names{iconfig}, Nm, Mr);

            errL2(icase, iconfig, id)   = rel_L2_error(f, fh, dtheta, normL2);
            errLinf(icase, iconfig, id) = rel_Linf_error(f, fh);
            penTV(icase, iconfig, id)   = oscillation_penalty_tv(f, fh);

            fprintf('  D = %2d -> L2 = %.3e, Linf = %.3e, TV = %.3e\n', ...
                D, ...
                errL2(icase, iconfig, id), ...
                errLinf(icase, iconfig, id), ...
                penTV(icase, iconfig, id));
        end
    end
end

%% ============================================================
%  Aggregate metrics
%% ============================================================

meanL2    = squeeze(mean(errL2, 1));
worstL2   = squeeze(max(errL2, [], 1));

meanLinf  = squeeze(mean(errLinf, 1));
worstLinf = squeeze(max(errLinf, [], 1));

meanTV    = squeeze(mean(penTV, 1));
worstTV   = squeeze(max(penTV, [], 1));

%% ============================================================
%  Print summary
%% ============================================================

fprintf('\n============================================================\n');
fprintf('Summary by metric\n');
fprintf('============================================================\n');

for iconfig = 1:nConfigs
    fprintf('\n%s\n', config_labels{iconfig});
    for id = 1:nD
        fprintf(['  D = %2d -> ', ...
                 'mean L2 = %.3e, worst L2 = %.3e, ', ...
                 'mean Linf = %.3e, worst Linf = %.3e, ', ...
                 'mean TV = %.3e, worst TV = %.3e\n'], ...
            D_list(id), ...
            meanL2(iconfig, id),   worstL2(iconfig, id), ...
            meanLinf(iconfig, id), worstLinf(iconfig, id), ...
            meanTV(iconfig, id),   worstTV(iconfig, id));
    end
end

%% ============================================================
%  Plot group 1: average metric vs DOF
%% ============================================================

figure;
for iconfig = 1:nConfigs
    semilogy(D_list, meanL2(iconfig, :), '-o', ...
        'LineWidth', 1.5, 'MarkerSize', 8); hold on;
end
grid on;
xlabel('Total velocity DOF');
ylabel('Mean relative L^2 error');
legend(config_labels, 'Location', 'southwest');
title('Average L^2 error vs DOF (MdeMat API)');

figure;
for iconfig = 1:nConfigs
    semilogy(D_list, meanLinf(iconfig, :), '-o', ...
        'LineWidth', 1.5, 'MarkerSize', 8); hold on;
end
grid on;
xlabel('Total velocity DOF');
ylabel('Mean relative L^\infty error');
legend(config_labels, 'Location', 'southwest');
title('Average L^\infty error vs DOF (MdeMat API)');

figure;
for iconfig = 1:nConfigs
    semilogy(D_list, meanTV(iconfig, :), '-o', ...
        'LineWidth', 1.5, 'MarkerSize', 8); hold on;
end
grid on;
xlabel('Total velocity DOF');
ylabel('Mean TV penalty');
legend(config_labels, 'Location', 'southwest');
title('Average TV penalty vs DOF (MdeMat API)');

%% ============================================================
%  Plot group 2: worst-case metric vs DOF
%% ============================================================

figure;
for iconfig = 1:nConfigs
    semilogy(D_list, worstL2(iconfig, :), '-o', ...
        'LineWidth', 1.5, 'MarkerSize', 8); hold on;
end
grid on;
xlabel('Total velocity DOF');
ylabel('Worst-case relative L^2 error');
legend(config_labels, 'Location', 'southwest');
title('Worst-case L^2 error vs DOF (MdeMat API)');

figure;
for iconfig = 1:nConfigs
    semilogy(D_list, worstLinf(iconfig, :), '-o', ...
        'LineWidth', 1.5, 'MarkerSize', 8); hold on;
end
grid on;
xlabel('Total velocity DOF');
ylabel('Worst-case relative L^\infty error');
legend(config_labels, 'Location', 'southwest');
title('Worst-case L^\infty error vs DOF (MdeMat API)');

figure;
for iconfig = 1:nConfigs
    semilogy(D_list, worstTV(iconfig, :), '-o', ...
        'LineWidth', 1.5, 'MarkerSize', 8); hold on;
end
grid on;
xlabel('Total velocity DOF');
ylabel('Worst-case TV penalty');
legend(config_labels, 'Location', 'southwest');
title('Worst-case TV penalty vs DOF (MdeMat API)');

%% ============================================================
%  Plot group 3: cross-case bars at representative DOF
%% ============================================================

D_show   = 25;
id_show  = find(D_list == D_show, 1);

if isempty(id_show)
    error('D_show must be contained in D_list.');
end

figure;
vals = squeeze(errL2(:, :, id_show));
bar(vals);
grid on;
set(gca, 'XTick', 1:nCases, 'XTickLabel', case_labels);
xtickangle(20);
ylabel('Relative L^2 error');
legend(config_labels, 'Location', 'northoutside', 'Orientation', 'horizontal');
title(sprintf('Cross-case L^2 error at DOF = %d (MdeMat API)', D_show));

figure;
vals = squeeze(errLinf(:, :, id_show));
bar(vals);
grid on;
set(gca, 'XTick', 1:nCases, 'XTickLabel', case_labels);
xtickangle(20);
ylabel('Relative L^\infty error');
legend(config_labels, 'Location', 'northoutside', 'Orientation', 'horizontal');
title(sprintf('Cross-case L^\\infty error at DOF = %d (MdeMat API)', D_show));

figure;
vals = squeeze(penTV(:, :, id_show));
bar(vals);
grid on;
set(gca, 'XTick', 1:nCases, 'XTickLabel', case_labels);
xtickangle(20);
ylabel('TV penalty');
legend(config_labels, 'Location', 'northoutside', 'Orientation', 'horizontal');
title(sprintf('Cross-case TV penalty at DOF = %d (MdeMat API)', D_show));

%% ============================================================
%  Plot group 4: fitting curves for each case at D_show
%% ============================================================

for icase = 1:nCases
    case_name = case_names{icase};
    f_func    = @(th) target_function(th, case_name, x0, y0, t0);
    f         = f_func(theta);

    figure;
    plot(theta, f, 'k', 'LineWidth', 2); hold on;

    for iconfig = 1:nConfigs
        Nm = Nm_store(iconfig, id_show);
        Mr = Mr_store(iconfig, id_show);

        fh = approximate_with_mdemat(f_func, theta, ...
            config_names{iconfig}, Nm, Mr);

        plot(theta, fh, 'LineWidth', 1.2);
    end

    grid on;
    xlabel('\theta');
    ylabel('f(\theta)');
    legend(['Exact', config_labels], 'Location', 'best');
    title(sprintf('Function fitting: %s, DOF = %d (MdeMat API)', ...
        case_labels{icase}, D_show));
end

%% ============================================================
%%  local functions
%% ============================================================

function [Nm, Mr] = choose_representation(D, config_name)

    switch config_name
        case 'pure_modal'
            modal_dof = D;
            Mr        = 0;
            Nm        = (modal_dof - 1) / 2;

        case 'pure_nodal'
            Nm = 0;
            Mr = D;

        case 'hybrid_1_4_3_4'
            target_modal_dof = D / 4;
            modal_dof        = nearest_odd_dof(D, target_modal_dof);
            Nm               = (modal_dof - 1) / 2;
            Mr               = D - modal_dof;

        case 'hybrid_1_2_1_2'
            target_modal_dof = D / 2;
            modal_dof        = nearest_odd_dof(D, target_modal_dof);
            Nm               = (modal_dof - 1) / 2;
            Mr               = D - modal_dof;

        case 'hybrid_3_4_1_4'
            target_modal_dof = 3 * D / 4;
            modal_dof        = nearest_odd_dof(D, target_modal_dof);
            Nm               = (modal_dof - 1) / 2;
            Mr               = D - modal_dof;

        otherwise
            error('Unknown configuration: %s', config_name);
    end
end

function modal_dof = nearest_odd_dof(D, target_modal_dof)
    candidates = 1:2:(D-2);
    [~, idx]   = min(abs(candidates - target_modal_dof));
    modal_dof  = candidates(idx);
end

function fh = approximate_with_mdemat(f_func, theta_ref, config_name, Nm, Mr)
% APPROXIMATE_WITH_MDEMAT Approximate f using the MdeMat element API.
%
%   f_func:    function handle f(theta) -> row vector
%   theta_ref: reference evaluation grid (1 x Nt_ref)
%   config_name, Nm, Mr: configuration parameters from choose_representation
%
%   Returns fh as (1 x Nt_ref) row vector.
%
%   NOTE: To change the decomposition method, edit fitHybrid below.

    nDims = 2;
    vRed  = 'topology';
    V_ref = [cos(theta_ref); sin(theta_ref)];   % (2 x Nt_ref)

    switch config_name

        case 'pure_modal'
            nu   = Nm + 1;
            elem = approx.element.L2SphereElement.modal( ...
                nDims, nu, reduction=vRed);
            C    = fitModal(elem, f_func);
            fh   = elem.eval(V_ref, C);

        case 'pure_nodal'
            elem = approx.element.L2SphereElement.nodal( ...
                nDims, Mr, reduction=vRed);
            C    = fitNodal(elem, f_func);
            fh   = elem.eval(V_ref, C);

        otherwise
            nu      = Nm + 1;
            macElem = approx.element.L2SphereElement.modal( ...
                nDims, nu, reduction=vRed);
            micElem = approx.element.L2SphereElement.nodal( ...
                nDims, Mr, reduction=vRed);
            [C_mac, C_mic] = fitHybrid(macElem, micElem, f_func);
            fh = macElem.eval(V_ref, C_mac) + micElem.eval(V_ref, C_mic);
    end
end

function C = fitModal(elem, fTheta)
% FITMODAL L2-project fTheta onto the modal element.
    V     = elem.Volume.Nodes;
    w     = elem.Volume.Weights;
    theta = mod(atan2(V(2,:), V(1,:)), 2*pi);
    f     = fTheta(theta);
    B     = elem.Approximator.Basis.eval(V);
    rhs   = B * (w(:) .* f(:));
    C     = elem.Approximator.Mass \ rhs;
end

function C = fitNodal(elem, fTheta)
% FITNODAL Fit fTheta at the nodal element's RBF nodes.
    V     = elem.Volume.Nodes;
    theta = mod(atan2(V(2,:), V(1,:)), 2*pi);
    f     = fTheta(theta);
    C     = elem.fit(f(:));
end

function [C_mac, C_mic] = fitHybrid(macElem, micElem, fTheta)
% FITHYBRID Hybrid decomposition: modal L2 projection + nodal residual fit.
%
%   NOTE: To modify the decomposition, edit this function.
%   Step 1: L2-project f onto the modal space -> macro coefficients C_mac
%   Step 2: Fit the residual (f - u_mac) at nodal points -> micro coefficients C_mic

    C_mac = fitModal(macElem, fTheta);

    V_n   = micElem.Volume.Nodes;
    th_n  = mod(atan2(V_n(2,:), V_n(1,:)), 2*pi);
    f_n   = fTheta(th_n);
    u_n   = macElem.eval(V_n, C_mac);
    g_n   = f_n - u_n;
    C_mic = micElem.fit(g_n(:));
end

function f = target_function(theta, case_name, x, y, t)

    switch case_name
        case 'smooth_transport'
            f = 1 + cos(pi*(x - cos(theta)*t)) .* cos(pi*(y - sin(theta)*t));

        case 'gaussian_beam_sigma_030'
            sigma = 0.30;
            r2    = (x - cos(theta)*t).^2 + (y - sin(theta)*t).^2;
            f     = exp(-r2/(2*sigma^2)) / (2*pi*sigma^2);

        case 'gaussian_beam_sigma_015'
            sigma = 0.15;
            r2    = (x - cos(theta)*t).^2 + (y - sin(theta)*t).^2;
            f     = exp(-r2/(2*sigma^2)) / (2*pi*sigma^2);

        case 'gaussian_beam_sigma_008'
            sigma = 0.08;
            r2    = (x - cos(theta)*t).^2 + (y - sin(theta)*t).^2;
            f     = exp(-r2/(2*sigma^2)) / (2*pi*sigma^2);

        case 'delta_beam_eps_020'
            f = periodic_gaussian(theta, pi/6, 0.20);

        case 'delta_beam_eps_010'
            f = periodic_gaussian(theta, pi/6, 0.10);

        case 'delta_beam_eps_005'
            f = periodic_gaussian(theta, pi/6, 0.05);

        otherwise
            error('Unknown case name: %s', case_name);
    end
end

function d = wrap_angle(x)
    d = mod(x + pi, 2*pi) - pi;
end

function f = periodic_gaussian(theta, theta0, eps)
    d      = wrap_angle(theta - theta0);
    f      = exp(-d.^2 / (2*eps^2));
    dth    = 2*pi / numel(theta);
    f      = f / (sum(f) * dth);
end

function err = rel_L2_error(f, fh, dtheta, normL2)
    err = sqrt(sum(abs(f - fh).^2) * dtheta) / normL2;
end

function err = rel_Linf_error(f, fh)
    err = max(abs(f - fh)) / (max(abs(f)) + 1e-14);
end

function tv = total_variation(u)
    tv = sum(abs(diff(u)));
end

function pen = oscillation_penalty_tv(f, fh)
    tv_f  = total_variation(f);
    tv_fh = total_variation(fh);
    pen   = max(tv_fh - tv_f, 0) / (tv_f + 1e-14);
end
