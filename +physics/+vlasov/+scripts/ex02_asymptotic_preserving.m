%==========================================================================
% FileName: ex02_asymptotic_preserving
% Author: Yi CAI
% Date: 27/01/2026
% Description:
%   Asymptotic preserving test.
%==========================================================================

clc, clear, close all;

%% EXECUTION
% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = physics.scheme.ConfigParser.parseFile(common);

% Checkpoint information
config.id = 'ex02i';
config.ckptDir = fullfile(fileparts(mfilename('fullpath')), config.id);
config.ckptPrefix = [upper(config.schemeName), num2str(config.xBasisOrder)];
if config.useFilter
    config.ckptPrefix = ['F', config.ckptPrefix];
end
switch config.id
    case 'ex02h'
        config.ckptPrefix = strjoin({'LDH', config.ckptPrefix}, '_');
    case 'ex02i'
        config.ckptPrefix = strjoin({'LDI', config.ckptPrefix}, '_');
end
config.tau0 = 10^5;
config.vBBox = [-3.2, 3.2];
config.ckptFigureName = 'f';

override = 0;
strategy1d = physics.visual.Strategy1d();
strategy2d = physics.visual.Strategy2d();
epsilonLevel = [0,-1,-2,-3,-4,-5,-6];
epsilon = 10.^epsilonLevel;
peLegs = cell(1, length(epsilon));
k0 = inf;
kmax = length(epsilon);
for k = 1:kmax
    config.epsilon = epsilon(k);
    config.tFinal = 100;
    config.ckptTimeStamps = [];
    config.ckptPostfix = sprintf('epsilon%.0e', epsilon(k));
    fileName = sprintf('%s_state_%s.mat', config.ckptPrefix, config.ckptPostfix);
    fileName = fullfile(config.ckptDir, fileName);

    if exist(fileName, 'file') && ~override
        state = physics.vlasov.HermiteState.load(fileName);
    else
        [config, scheme, state] = run(config);
        state.save(fileName);
    end

    % Save distribution plot
    figure(2);
    ax = gca;
    xRef = 0;
    V = linspace(config.vBBox(1), config.vBBox(2), 512);
    X = state.XDisc.Mesh.collocate(xRef);
    F = state.distribution(xRef, V);
    imagesc(ax, X, V, F.', 'Interpolation', 'bilinear');
    axis(ax, 'xy', 'equal', 'tight');
    xl = xlabel(ax, 'x');
    yl = ylabel(ax, 'v');
    set(xl, 'FontSize', strategy2d.DefaultAxisLabelFontSize);
    set(yl, 'FontSize', strategy2d.DefaultAxisLabelFontSize);
    colormap(ax, strategy2d.ColorMap);
    colorbar(ax);
    fileName = sprintf('%s_f_%s.pdf', config.ckptPrefix, config.ckptPostfix);
    fileName = fullfile(config.ckptDir, fileName);
    exportgraphics(ax, fileName, 'ContentType', 'vector');

    % Save potential energy plot
    figure(3);
    ax = gca;
    hold(ax, 'on');
    t = state.History.time;
    potentialEnergy = max(1e-35, state.History.potentialEnergy);
    % mask = state.History.time <= config.tFinal;
    % t = state.History.time(mask) / epsilon(k);
    % potentialEnergy = max(1e-35, state.History.potentialEnergy(mask));
    if k <= k0
        style = strategy1d.getDefaultLineStyle(k+1, k);
    else
        style = strategy1d.getDefaultScatterStyle(k+1, k-k0);
    end
    line = plot(t, potentialEnergy, style{:});
    set(line, 'LineWidth', 3, 'MarkerSize', 8);
    set(ax, 'YScale', 'log');
    hold(ax, 'off');
    peLegs{k} = sprintf('$\\varepsilon=10^{%d}$', epsilonLevel(k));
    leg = legend(ax, peLegs(1:k));
    set(leg, 'Location', strategy1d.DefaultLegendPosition);
    set(leg, 'FontSize', strategy1d.DefaultLegendFontSize);
    set(leg, 'Interpreter', 'latex');
    xl = xlabel('$t$');
    % xl = xlabel('$\varepsilon t$');
    set(xl, 'Interpreter','latex', 'FontSize', 16);
    ax.XAxis.Exponent = 0;
    ax.XAxis.TickLabelFormat = '%.0f';
    xlim([0, config.tFinal]);
    % xlim([0, config.tFinal/epsilon(k)]);
    ylim([1e-32, 1e2]);
    yticks(logspace(-35, 0, 8));
    yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-35, 0, 8), 'Un', 0));
    title('Potential energy', 'FontSize', 16);
    fileName = sprintf('%s_potential_energy.pdf', config.ckptPrefix);
    % fileName = sprintf('%s_rescaled_potential_energy.pdf', config.ckptPrefix);
    fileName = fullfile(config.ckptDir, fileName);
    exportgraphics(ax, fileName, 'ContentType', 'vector');

    % Save relative entropy plot
    figure(4);
    ax = gca;
    hold(ax, 'on');
    t = state.History.time;
    L2FDistance = max(1e-35, state.History.L2FDistance);
    % mask = state.History.time <= config.tFinal;
    % t = state.History.time(mask) / epsilon(k);
    % L2FDistance = max(1e-35, state.History.L2FDistance(mask));
    if k <= k0
        style = strategy1d.getDefaultLineStyle(k+1, k);
    else
        style = strategy1d.getDefaultScatterStyle(k+1, k-k0);
    end
    line = plot(t, L2FDistance, style{:});
    set(line, 'LineWidth', 3, 'MarkerSize', 8);
    peLegs{k} = sprintf('$\\varepsilon=10^{%d}$', epsilonLevel(k));
    leg = legend(ax, peLegs(1:k));
    set(leg, 'Location', strategy1d.DefaultLegendPosition);
    set(leg, 'FontSize', strategy1d.DefaultLegendFontSize);
    set(leg, 'Interpreter', 'latex');
    set(ax, 'YScale', 'log');
    hold(ax, 'off');
    xlim([0, config.tFinal]);
    % xlim([0, config.tFinal/epsilon(k)]);
    xl = xlabel('$t$');
    % xl = xlabel('$\varepsilon t$');
    set(xl, 'Interpreter','latex', 'FontSize', 16);
    ax.XAxis.Exponent = 0;
    ax.XAxis.TickLabelFormat = '%.0f';
    ylim([5e-12, 2e-1]);
    yticks(logspace(-16, 0, 17));
    yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-16, 0, 17), 'Un', 0));
    ti = title('$\|f-f_\infty\|_{L^2(f_\infty^{-1})}$');
    set(ti, 'Interpreter','latex', 'FontSize', 16);
    fileName = sprintf('%s_relative_entropy.pdf', config.ckptPrefix);
    % fileName = sprintf('%s_rescaled_relative_entropy.pdf', config.ckptPrefix);
    fileName = fullfile(config.ckptDir, fileName);
    exportgraphics(ax, fileName, 'ContentType', 'vector');

    % Save free energy plot
    figure(5);
    ax = gca;
    hold(ax, 'on');
    t = state.History.time;
    freeEnergy = max(1e-35, L2FDistance.^2 / 2 + potentialEnergy);
    % mask = state.History.time <= config.tFinal;
    % t = state.History.time(mask) / epsilon(k);
    % freeEnergy = max(1e-35, freeEnergy(mask));
    if k <= k0
        style = strategy1d.getDefaultLineStyle(k+1, k);
    else
        style = strategy1d.getDefaultScatterStyle(k+1, k-k0);
    end
    line = plot(t, freeEnergy, style{:});
    set(line, 'LineWidth', 3, 'MarkerSize', 8);
    peLegs{k} = sprintf('$\\varepsilon=10^{%d}$', epsilonLevel(k));
    leg = legend(ax, peLegs(1:k));
    set(leg, 'Location', strategy1d.DefaultLegendPosition);
    set(leg, 'FontSize', strategy1d.DefaultLegendFontSize);
    set(leg, 'Interpreter', 'latex');
    set(ax, 'YScale', 'log');
    hold(ax, 'off');
    xlim([0, config.tFinal]);
    % xlim([0, config.tFinal/epsilon(k)]);
    xl = xlabel('$t$');
    % xl = xlabel('$\varepsilon t$');
    set(xl, 'Interpreter','latex', 'FontSize', 16);
    ax.XAxis.Exponent = 0;
    ax.XAxis.TickLabelFormat = '%.0f';
    ylim([1e-27, 1e2]);
    yticks(logspace(-35, 0, 8));
    yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-35, 0, 8), 'Un', 0));
    ti = title('$\mathcal{E}_h$');
    set(ti, 'Interpreter','latex', 'FontSize', 16);
    fileName = sprintf('%s_free_energy.pdf', config.ckptPrefix);
    % fileName = sprintf('%s_rescaled_free_energy.pdf', config.ckptPrefix);
    fileName = fullfile(config.ckptDir, fileName);
    exportgraphics(ax, fileName, 'ContentType', 'vector');
end



%% SIMULATION

function [config, scheme, state] = run(config)

config.nDims = 1; % Dimension
switch config.id
    case 'ex02h'
        config.L = 4*pi; % Length of domain
        config.kappa = 2*pi / config.L; % Wave length
        config.xBBox = [0, config.L]; % Spatial bounding box
    case 'ex02i'
        config.L = 6; % Length of domain
        config.kappa = pi / config.L; % Wave length
        config.xBBox = [-config.L, config.L]; % Spatial bounding box
end
config.T0 = 1; % Temperature parameter
config.delta = 0.01; % Perturbation
config.nx = 128; % Grid resolution
config.nh = 640; % Number of Hermite modes - 1
config.tau = config.tau0 * config.epsilon;
config.phiInf = phiInf(config); % Electric potential equilibrium
config.EInf = EInf(config); % Electric field equilibrium
config.cInf = cInf(config); % Normalization constant
config.rhoInf = rhoInf(config); % Density equilibrium
config.sqrtRhoInf = sqrtRhoInf(config); % Square root of density equilibrium
config.dSqrtRhoInf = dSqrtRhoInf(config); % Derivative of square root of density equilibrium
config.ic = DInit(config); % Initial hermite coefficients
config.bc = []; % Boundary condition
config.exact = []; % Exact solution
config.verbose = 1; % Verbose flag

% Construct and configure scheme
scheme = physics.vlasov.ApDghsScheme(config = config);

% Configure visualizer
visualizer = scheme.getConfig('visualizer');
visualizer.setDensity(1);
visualizer.setTimeline(scheme.getConfig('tFinal'), 1);
visualizer.setComponents(struct('D', [], 'Psi', [], 'Q', []));
visualizer.addDataset(scheme.getConfig('sId'), 'numeric');
visualizer.addPlotter('profile', '1d');

vElement = approx.element.L2OrthotopeElement.hermite(1, config.nh, T = config.T0);
vDisc = approx.space.SpectralSpace(vElement);

switch scheme.getConfig('xBasisType')
    case 'modal'
        xElement = approx.element.BH1OrthotopeElement.modal(config.nDims, ...
            scheme.getConfig('xBasisOrder'), ...
            pattern=scheme.getConfig('xBasisPattern'));
    case 'nodal'
        xElement = approx.element.BH1OrthotopeElement.nodal(config.nDims, ...
            scheme.getConfig('xBasisOrder'));
end
xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx/2);
xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
state = physics.vlasov.ApHermiteState(xDisc, vDisc);
state = scheme.run(state);
end

%% HELPER FUNCTIONS

function h = DInit(config)
h = @(x, varargin) DInitImpl(config, x, varargin{:});
end

function D = DInitImpl(config, x, varargin)
if nargin >= 1, nh = varargin{1}; end
kappa = config.kappa;
delta = config.delta;
rhoInf = rhoInfImpl(config, x);
D0 = (rhoInf + delta * cos(kappa*x)) ./ sqrt(rhoInf);
D = [D0, zeros(1, size(x, 2)*(nh-1))];
end

function h = phiInf(config)
h = @(x) phiInfImpl(config, x);
end

function phi = phiInfImpl(config, x)
switch config.id
    case 'ex02h'
        phi = zeros(1, size(x, 2));
    case 'ex02i'
        kappa = config.kappa;
        phi = 0.2*sin(kappa*x);
end
end

function h = rhoInf(config)
h = @(x) rhoInfImpl(config, x);
end

function rho = rhoInfImpl(config, x)
cInf = config.cInf;
T0 = config.T0;
rho = cInf * exp(-phiInfImpl(config, x) / T0);
end

function h = sqrtRhoInf(config)
h = @(x) sqrtRhoInfImpl(config, x);
end

function rho = sqrtRhoInfImpl(config, x)
rho = sqrt(rhoInfImpl(config, x));
end

function h = EInf(config)
% Automatic symbolic differentiation of sqrt(rhoInf)
syms xs;
phiSym = -phiInfImpl(config, xs);
ESym = diff(phiSym, xs);

% Check if derivative is a constant (independent of xs)
if ~has(ESym, xs)
    % Derivative is constant - convert to numeric value
    constVal = double(ESym);
    % Return function that returns constant with same size as input
    h = @(x) constVal * ones(1, size(x, 2));
else
    % Derivative depends on x - convert symbolic expression to function handle
    h = matlabFunction(ESym, 'Vars', xs);
end
end

function h = dSqrtRhoInf(config)
% Automatic symbolic differentiation of sqrt(rhoInf)
syms xs;
cInf = config.cInf;
T0 = config.T0;
phiSym = phiInfImpl(config, xs);
sqrtRhoSym = sqrt(cInf) * exp(-phiSym / (2*T0));
dSqrtRhoSym = diff(sqrtRhoSym, xs);

% Check if derivative is a constant (independent of xs)
if ~has(dSqrtRhoSym, xs)
    % Derivative is constant - convert to numeric value
    constVal = double(dSqrtRhoSym);
    % Return function that returns constant with same size as input
    h = @(x) constVal * ones(1, size(x, 2));
else
    % Derivative depends on x - convert symbolic expression to function handle
    h = matlabFunction(dSqrtRhoSym, 'Vars', xs);
end
end

function c = cInf(config)
L = config.L;
T0 = config.T0;
phiInf = config.phiInf;
c = 2 * L / integral(@(x) exp(-phiInf(x)/T0), -L, L);
end