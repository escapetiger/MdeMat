%==========================================================================
% FileName: ex00_func_approx
% Author: Yi CAI
% Description:
%   Angular function approximation comparison on S^1 (2DT topology
%   reduction). Evaluates the representation quality of pure modal (P_N),
%   pure nodal (S_N), and hybrid (P_N + S_N) velocity bases for four
%   representative test functions: smooth, single beam, two beams, mixed.
%
%   Run this script before any radiation scheme examples to understand
%   what angular approximation can and cannot capture.
%==========================================================================

clc; clear; close all;

%% CONFIGURATION

nDims  = 2;          % physical dimensions of velocity space
reduc  = 'topology'; % dimension reduction type (2D sphere → S^1)
kappa  = 15;         % beam concentration (larger = narrower beam)
N_ref  = 1024;       % fine reference grid on S^1
outDir = fullfile(fileparts(mfilename('fullpath')), 'ex00_func_approx');
if ~isfolder(outDir), mkdir(outDir); end

%% REFERENCE GRID

theta_ref = (0:N_ref-1) * (2*pi / N_ref);   % 1 × N_ref, equispaced in [0, 2*pi)
V_ref     = [cos(theta_ref); sin(theta_ref)]; % 2 × N_ref, for element eval

%% TEST FUNCTIONS ON S^1  (all defined as f(theta) with theta in [0, 2*pi))

E   = 1 / (2*pi);   % isotropic equilibrium density: int f dtheta = 1 → E = 1/(2*pi)
th0 = pi / 6;        % beam 1 peak angle
th1 = 5*pi / 6;      % beam 2 peak angle

% Von Mises kernel: exp(kap*(cos(theta - theta0) - 1)) / Z
%   Z = 2*pi * besseli(0,kap,1)  [scaled Bessel avoids numerical overflow]
%   Integral over [0, 2*pi): 1   (normalized angular distribution)
vmNorm   = @(kap)       2*pi * besseli(0, kap, 1);
vonMises = @(th, th0, kap) exp(kap * (cos(th - th0) - 1)) / vmNorm(kap);

% Test 1: Smooth (Fourier polynomial up to degree 3 = cos/sin with k <= 3)
%   int f dtheta = 1;  few Fourier modes → modal converges exponentially
fSmooth = @(th) E * (1 + 0.5*cos(th) + 0.3*cos(2*th) + 0.2*sin(th));

% Test 2: Single narrow Gaussian beam at theta = pi/6
%   int f dtheta = 1;  modal suffers Gibbs oscillations for large kappa
fBeam = @(th) vonMises(th, th0, kappa);

% Test 3: Two narrow Gaussian beams at theta = pi/6 and 5*pi/6 (equal weight)
%   int f dtheta = 1;  needs both angular resolution and smoothness handling
fTwoBeam = @(th) 0.5 * vonMises(th, th0, kappa) + 0.5 * vonMises(th, th1, kappa);

% Test 4: Smooth background + localized beam — the key benchmark for hybrid
%   Modal captures the smooth part; nodal resolves the beam residual
fMixed = @(th) E + 0.3*E*cos(th) + 0.2*E*cos(2*th) + 0.5 * vonMises(th, th0, kappa);

testFuns  = {fSmooth, fBeam, fTwoBeam, fMixed};
testNames = {'Smooth', 'Single Beam', 'Two Beams', 'Mixed'};
nTests    = numel(testFuns);

%% DOF CONFIGURATIONS

% Modal: nu = Fourier order, NDofs = 2*nu - 1
%   nu=1 → 1 DOF (constant), nu=2 → 3 DOFs, nu=k → 2k-1 DOFs
nuList     = [2, 3, 4, 5, 6, 7, 8, 9, 10, 12];
nDofsModal = 2*nuList - 1;

% Nodal: nv equispaced RBF nodes on S^1, NDofs = nv
nvList     = [4, 6, 8, 12, 16, 20, 24, 32];
nDofsNodal = nvList;

% Hybrid: (nu_macro, nv_micro) pairs
%   Total DOFs = (2*nu - 1) + nv_micro
%   Modal part captures smooth low-frequency content
%   Nodal part captures residual high-frequency / localized features
hybridCfg = [2,  4;   %  3 +  4 = 7  total
             2,  8;   %  3 +  8 = 11 total
             3,  6;   %  5 +  6 = 11 total
             3, 10;   %  5 + 10 = 15 total
             4,  8;   %  7 +  8 = 15 total
             5,  8;   %  9 +  8 = 17 total
             5, 12;   %  9 + 12 = 21 total
             6, 12;   % 11 + 12 = 23 total
             8,  8;   % 15 +  8 = 23 total
             8, 16];  % 15 + 16 = 31 total
nDofsHybrid = (2*hybridCfg(:,1) - 1) + hybridCfg(:,2);

%% COMPUTE APPROXIMATION ERRORS

errModal  = zeros(nTests, length(nuList));
errNodal  = zeros(nTests, length(nvList));
errHybrid = zeros(nTests, size(hybridCfg, 1));

for it = 1:nTests
    fTest   = testFuns{it};
    f_exact = fTest(theta_ref);                         % 1 × N_ref
    f_norm  = l2norm(f_exact, N_ref);

    fprintf('Test %d/%d: %s\n', it, nTests, testNames{it});

    % Modal sweep
    for im = 1:length(nuList)
        nu   = nuList(im);
        elem = approx.element.L2SphereElement.modal(nDims, nu, reduction=reduc);
        C    = fitModal(elem, fTest);
        f_h  = elem.eval(V_ref, C);
        errModal(it, im) = l2norm(f_h - f_exact, N_ref) / f_norm;
    end

    % Nodal sweep
    for in_ = 1:length(nvList)
        nv   = nvList(in_);
        elem = approx.element.L2SphereElement.nodal(nDims, nv, reduction=reduc);
        C    = fitNodal(elem, fTest);
        f_h  = elem.eval(V_ref, C);
        errNodal(it, in_) = l2norm(f_h - f_exact, N_ref) / f_norm;
    end

    % Hybrid sweep
    for ih = 1:size(hybridCfg, 1)
        nu_m   = hybridCfg(ih, 1);
        nv_n   = hybridCfg(ih, 2);
        macElem = approx.element.L2SphereElement.modal(nDims, nu_m, reduction=reduc);
        micElem = approx.element.L2SphereElement.nodal(nDims, nv_n, reduction=reduc);
        [C_mac, C_mic] = fitHybrid(macElem, micElem, fTest);
        f_h = macElem.eval(V_ref, C_mac) + micElem.eval(V_ref, C_mic);
        errHybrid(it, ih) = l2norm(f_h - f_exact, N_ref) / f_norm;
    end
end

%% PROFILE PLOTS (fixed DOF count ~15)
%
%   Purpose: visualise approximation quality at the same DOF budget.
%   Modal:  nu = 8  → 15 DOFs  (Fourier modes 0 … 7)
%   Nodal:  nv = 16 → 16 DOFs  (16 equispaced RBF nodes)
%   Hybrid: nu = 4, nv = 8 → 7 + 8 = 15 DOFs

nu_prof = 8;   % modal order for profile plot
nv_prof = 16;  % nodal count for profile plot
nu_h    = 4;   % hybrid macro order (7 DOFs)
nv_h    = 8;   % hybrid micro count (8 DOFs)

strategy = physics.visual.Strategy1d();
colors   = strategy.ColorOrder;

for it = 1:nTests
    fTest   = testFuns{it};
    f_exact = fTest(theta_ref);

    % Modal approximation
    elemM = approx.element.L2SphereElement.modal(nDims, nu_prof, reduction=reduc);
    f_m   = elemM.eval(V_ref, fitModal(elemM, fTest));

    % Nodal approximation
    elemN = approx.element.L2SphereElement.nodal(nDims, nv_prof, reduction=reduc);
    f_n   = elemN.eval(V_ref, fitNodal(elemN, fTest));

    % Hybrid approximation
    elemMac = approx.element.L2SphereElement.modal(nDims, nu_h, reduction=reduc);
    elemMic = approx.element.L2SphereElement.nodal(nDims, nv_h, reduction=reduc);
    [C_mac, C_mic] = fitHybrid(elemMac, elemMic, fTest);
    f_h = elemMac.eval(V_ref, C_mac) + elemMic.eval(V_ref, C_mic);

    fig = figure(it);
    clf(fig);
    ax = axes(fig);
    hold(ax, 'on');

    plot(ax, theta_ref, f_exact, '-',  'Color', colors(1,:), 'LineWidth', 2.0, ...
        'DisplayName', 'Exact');
    plot(ax, theta_ref, f_m,     '--', 'Color', colors(2,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Modal P_{%d}  (%d DOF)', nu_prof-1, 2*nu_prof-1));
    plot(ax, theta_ref, f_n,     ':',  'Color', colors(3,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Nodal S_{%d}  (%d DOF)', nv_prof, nv_prof));
    plot(ax, theta_ref, f_h,     '-.', 'Color', colors(4,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Hybrid P_{%d}+S_{%d}  (%d DOF)', ...
                               nu_h-1, nv_h, 2*nu_h-1+nv_h));

    xlim(ax, [0, 2*pi]);
    xticks(ax, 0:pi/4:2*pi);
    xticklabels(ax, {'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'});
    xl = xlabel(ax, '\theta');
    yl = ylabel(ax, 'f(\theta)');
    legend(ax, 'show', 'FontSize', strategy.DefaultLegendFontSize, 'Location', 'best');
    set(ax, 'FontSize', strategy.DefaultFontSize);
    set(xl, 'FontSize', strategy.DefaultAxisLabelFontSize);
    set(yl, 'FontSize', strategy.DefaultAxisLabelFontSize);
    box(ax, 'on');

    fname = fullfile(outDir, sprintf('profile_%s.pdf', lower(strrep(testNames{it}, ' ', '_'))));
    exportgraphics(ax, fname, 'ContentType', 'vector');
    fprintf('  Saved: %s\n', fname);
end

%% CONVERGENCE PLOTS (relative L2 error vs total angular DOFs)

for it = 1:nTests
    fig = figure(nTests + it);
    clf(fig);
    ax = axes(fig);
    hold(ax, 'on');

    semilogy(ax, nDofsModal,  errModal(it,:),  'o-', ...
             'Color', colors(2,:), 'LineWidth', 1.5, 'MarkerSize', 7, ...
             'DisplayName', 'Modal P_N');
    semilogy(ax, nDofsNodal,  errNodal(it,:),  's-', ...
             'Color', colors(3,:), 'LineWidth', 1.5, 'MarkerSize', 7, ...
             'DisplayName', 'Nodal S_N');
    semilogy(ax, nDofsHybrid, errHybrid(it,:), '^-', ...
             'Color', colors(4,:), 'LineWidth', 1.5, 'MarkerSize', 7, ...
             'DisplayName', 'Hybrid P_N+S_N');

    xl = xlabel(ax, 'Total angular DOFs');
    yl = ylabel(ax, 'Relative L^2 error');
    legend(ax, 'show', 'FontSize', strategy.DefaultLegendFontSize, 'Location', 'best');
    grid(ax, 'on');
    set(ax, 'FontSize', strategy.DefaultFontSize, 'YScale', 'log');
    set(xl, 'FontSize', strategy.DefaultAxisLabelFontSize);
    set(yl, 'FontSize', strategy.DefaultAxisLabelFontSize);
    box(ax, 'on');

    fname = fullfile(outDir, sprintf('convergence_%s.pdf', lower(strrep(testNames{it}, ' ', '_'))));
    exportgraphics(ax, fname, 'ContentType', 'vector');
    fprintf('  Saved: %s\n', fname);
end

%% SAVE NUMERICAL RESULTS

results.testNames  = testNames;
results.nuList     = nuList;
results.nDofsModal = nDofsModal;
results.nvList     = nvList;
results.nDofsNodal = nDofsNodal;
results.hybridCfg  = hybridCfg;
results.nDofsHybrid = nDofsHybrid;
results.errModal   = errModal;
results.errNodal   = errNodal;
results.errHybrid  = errHybrid;
results.kappa      = kappa;
save(fullfile(outDir, 'results.mat'), 'results');
fprintf('\nResults saved to %s\n', outDir);

%% LOCAL FUNCTIONS

function C = fitModal(elem, fTheta)
% FITMODAL L2-project fTheta(theta) onto modal (Fourier) basis.
%
%   fTheta is a function handle f(theta) with theta in [0, 2*pi).
%   Element quadrature nodes (unit vectors) are converted to angles
%   via theta = atan2(v(2), v(1)), mapped to [0, 2*pi).
%
%   C = M^{-1} * B * W * f   (Galerkin L2 projection)

    V      = elem.Volume.Nodes;              % 2 × Nq  (unit vectors)
    w      = elem.Volume.Weights;            % 1 × Nq  (weights, includes 1/E)
    theta  = mod(atan2(V(2,:), V(1,:)), 2*pi); % 1 × Nq  (angles in [0, 2*pi))
    f      = fTheta(theta);                  % 1 × Nq
    B      = elem.Approximator.Basis.eval(V); % NDofs × Nq
    rhs    = B * (w(:) .* f(:));             % NDofs × 1
    C      = elem.Approximator.Mass \ rhs;   % NDofs × 1
end

function C = fitNodal(elem, fTheta)
% FITNODAL Fit fTheta(theta) at nodal (RBF) points via interpolation.
%
%   fTheta is a function handle f(theta) with theta in [0, 2*pi).
%   Solves the RBF system D*C = f at the equispaced nodal angles.

    V      = elem.Volume.Nodes;              % 2 × Nn  (unit vectors)
    theta  = mod(atan2(V(2,:), V(1,:)), 2*pi); % 1 × Nn
    f      = fTheta(theta);                  % 1 × Nn
    C      = elem.fit(f(:));                 % Nn × 1
end

function [C_mac, C_mic] = fitHybrid(macElem, micElem, fTheta)
% FITHYBRID Decompose fTheta into modal macro part u and nodal residual g.
%
%   fTheta is a function handle f(theta) with theta in [0, 2*pi).
%
%   Step 1: L2-project f onto the modal (macro) space → C_mac
%   Step 2: Evaluate residual r = f - u at the nodal (micro) angles
%   Step 3: Fit the residual with the nodal RBF basis → C_mic
%
%   Reconstruction:  f(theta) ≈ u(theta) + g(theta)
%
%   NOTE: This function defines the hybrid coupling strategy.
%   To modify the decomposition, edit this function and the reconstruction
%   in the calling code (macElem.eval + micElem.eval).

    % Step 1: modal projection of f(theta)
    C_mac = fitModal(macElem, fTheta);

    % Step 2–3: nodal fit of residual at micro nodal angles
    V_n   = micElem.Volume.Nodes;                  % 2 × Nn
    th_n  = mod(atan2(V_n(2,:), V_n(1,:)), 2*pi);  % 1 × Nn
    f_n   = fTheta(th_n);                          % 1 × Nn  (exact)
    u_n   = macElem.eval(V_n, C_mac);              % 1 × Nn  (modal part)
    g_n   = f_n - u_n;                             % 1 × Nn  (residual)
    C_mic = micElem.fit(g_n(:));                   % Nn × 1
end

function n = l2norm(f, N)
% L2NORM Approximate L2 norm on [0, 2*pi) using trapezoidal quadrature.
%
%   n = sqrt( (2*pi/N) * sum(f.^2) )

    n = sqrt(sum(f.^2) * (2*pi / N));
end
