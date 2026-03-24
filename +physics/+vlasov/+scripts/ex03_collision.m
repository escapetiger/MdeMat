%==========================================================================
% FileName: ex03_collision
% Author: Yi CAI
% Date: 27/01/2026
% Description:
%   Collision comparison.
%==========================================================================

clc, clear, close all;

%% EXECUTION
% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = physics.scheme.ConfigParser.parseFile(common);

% Checkpoint information
config.id = 'ex03ldh';
config.ckptDir = fullfile(fileparts(mfilename('fullpath')), config.id);
config.ckptPrefix = [upper(config.schemeName), num2str(config.xBasisOrder)];
if config.useFilter
    config.ckptPrefix = ['F', config.ckptPrefix];
end
switch config.id
    case 'ex03ldh'
        config.ckptPrefix = strjoin({'LDH', config.ckptPrefix}, '_');
        config.tFinal = 50;
        config.ckptTimeStamps = [4, 8, 16, 32, 40, 50];
    case 'ex03ldi'
        config.ckptPrefix = strjoin({'LDI', config.ckptPrefix}, '_');
        config.tFinal = 70;
        config.ckptTimeStamps = [4, 8, 16, 30, 50, 70];
        % config.tFinal = 400;
        % config.ckptTimeStamps = [];
    case 'ex03tsh'
        config.ckptPrefix = strjoin({'TSH', config.ckptPrefix}, '_');
        config.tFinal = 50;
        config.ckptTimeStamps = [4, 8, 16, 30, 50];
    case 'ex03tsi'
        config.ckptPrefix = strjoin({'TSI', config.ckptPrefix}, '_');
        config.tFinal = 60;
        config.ckptTimeStamps = [8, 16, 30, 60];
end
config.epsilon = 1;
config.vBBox = [-3.2, 3.2];
config.ckptFigureName = 'f';
config.verbose = 2;

override = 1;
strategy1d = physics.visual.Strategy1d();
strategy2d = physics.visual.Strategy2d();
tau0Level = [inf];
tau0 = 10.^tau0Level;
peLegs = cell(1, length(tau0));
k0 = inf;
kmax = length(tau0);
for k = 1:kmax
    config.tau0 = tau0(k);
    config.ckptPostfix = sprintf('tau0%.0e', tau0(k));
    % fileName = sprintf('%s_state_%s.mat', config.ckptPrefix, config.ckptPostfix);
    tckpt = 40;
    fileName = sprintf('%s_state_%s_t%dp0.mat', config.ckptPrefix, config.ckptPostfix, tckpt);
    fileName = fullfile(config.ckptDir, fileName);

    if exist(fileName, 'file') && ~override
        state = physics.vlasov.ApHermiteState.load(fileName);
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
    if config.verbose == 2
        imagesc(ax, X, V, F.', 'Interpolation', 'bilinear');
    else
        FInf = state.equilibrium(xRef, V);
        imagesc(ax, X, V, (F-FInf).', 'Interpolation', 'bilinear');
    end
    axis(ax, 'xy', 'equal', 'tight');
    xl = xlabel(ax, 'x');
    yl = ylabel(ax, 'v');
    set(ax, 'FontSize', strategy2d.DefaultFontSize);
    set(xl, 'FontSize', strategy2d.DefaultAxisLabelFontSize);
    set(yl, 'FontSize', strategy2d.DefaultAxisLabelFontSize);
    colormap(ax, strategy2d.ColorMap);
    colorbar(ax);
    ti = title(sprintf('t = %d', tckpt));
    set(ti, 'FontSize', strategy2d.DefaultTitleFontSize);
    caxis([-0.05, 0.65]);
    % fileName = sprintf('%s_f_%s_t%dp0.pdf', config.ckptPrefix, config.ckptPostfix, tckpt);
    fileName = sprintf('%s_f_%s.pdf', config.ckptPrefix, config.ckptPostfix);
    fileName = fullfile(config.ckptDir, fileName);
    exportgraphics(ax, fileName, 'ContentType', 'vector');

    % % Save free energy plot
    % figure(3);
    % ax = gca;
    % hold(ax, 'on');
    % t = state.History.time;
    % potentialEnergy = sqrt(2*state.History.potentialEnergy);
    % freeEnergy = max(1e-35, potentialEnergy);
    % if k <= k0
    %     style = strategy1d.getDefaultLineStyle(k+1, k);
    % else
    %     style = strategy1d.getDefaultScatterStyle(k+1, k-k0);
    % end
    % line = plot(t, freeEnergy, style{:});
    % set(line, 'LineWidth', 3, 'MarkerSize', 8);
    % peLegs{k} = sprintf('$\\tau_0=10^{%d}$', tau0Level(k));
    % leg = legend(ax, peLegs(1:k));
    % set(leg, 'Location', strategy1d.DefaultLegendPosition);
    % set(leg, 'FontSize', strategy1d.DefaultLegendFontSize);
    % set(leg, 'Interpreter', 'latex');
    % set(ax, 'YScale', 'log');
    % hold(ax, 'off');
    % xlim([0, config.tFinal]);
    % xl = xlabel('t');
    % ax.XAxis.Exponent = 0;
    % ax.XAxis.TickLabelFormat = '%.0f';
    % ylim([5e-8, 5e2]);
    % yticks(logspace(-18, 2, 21));
    % yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-18, 2, 21), 'Un', 0));
    % set(ax, 'FontSize', strategy2d.DefaultFontSize);
    % set(xl, 'FontSize', strategy2d.DefaultAxisLabelFontSize);
    % ti = title('Potential energy');
    % set(ti, 'FontSize', strategy2d.DefaultTitleFontSize);
    % fileName = sprintf('%s_potential_energy.pdf', config.ckptPrefix);
    % fileName = fullfile(config.ckptDir, fileName);
    % exportgraphics(ax, fileName, 'ContentType', 'vector');
    % 
    % % Save relative entropy plot
    % figure(4);
    % ax = gca;
    % hold(ax, 'on');
    % t = state.History.time;
    % L2FDistance = max(1e-35, state.History.L2FDistance);
    % if k <= k0
    %     style = strategy1d.getDefaultLineStyle(k+1, k);
    % else
    %     style = strategy1d.getDefaultScatterStyle(k+1, k-k0);
    % end
    % line = plot(t, L2FDistance, style{:});
    % set(line, 'LineWidth', 3, 'MarkerSize', 8);
    % peLegs{k} = sprintf('$\\tau_0=10^{%d}$', tau0Level(k));
    % leg = legend(ax, peLegs(1:k));
    % set(leg, 'Location', strategy1d.DefaultLegendPosition);
    % set(leg, 'FontSize', strategy1d.DefaultLegendFontSize);
    % set(leg, 'Interpreter', 'latex');
    % set(ax, 'YScale', 'log');
    % hold(ax, 'off');
    % xlim([0, config.tFinal]);
    % xl = xlabel('t');
    % ax.XAxis.Exponent = 0;
    % ax.XAxis.TickLabelFormat = '%.0f';
    % ylim([5e-5, 5e2]);
    % yticks(logspace(-16, 4, 21));
    % yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-16, 4, 21), 'Un', 0));
    % set(ax, 'FontSize', strategy2d.DefaultFontSize);
    % set(xl, 'FontSize', strategy2d.DefaultAxisLabelFontSize);
    % % ti = title('$\|f-f_\infty\|_{L^2(f_\infty^{-1})}$');
    % % set(ti, 'Interpreter','latex', 'FontSize', 16);
    % ti = title('Distance to kinetic equilibrium');
    % set(ti, 'FontSize', strategy2d.DefaultTitleFontSize);
    % fileName = sprintf('%s_relative_entropy.pdf', config.ckptPrefix);
    % fileName = fullfile(config.ckptDir, fileName);
    % exportgraphics(ax, fileName, 'ContentType', 'vector');
    % 
    % % Save trend to macro equilibrium plot
    % figure(5);
    % ax = gca;
    % hold(ax, 'on');
    % t = state.History.time;
    % L2RhoDistance = max(1e-35, state.History.L2RhoDistance);
    % if k <= k0
    %     style = strategy1d.getDefaultLineStyle(k+1, k);
    % else
    %     style = strategy1d.getDefaultScatterStyle(k+1, k-k0);
    % end
    % line = plot(t, L2RhoDistance, style{:});
    % set(line, 'LineWidth', 3, 'MarkerSize', 8);
    % peLegs{k} = sprintf('$\\tau_0=10^{%d}$', tau0Level(k));
    % leg = legend(ax, peLegs(1:k));
    % set(leg, 'Location', strategy1d.DefaultLegendPosition);
    % set(leg, 'FontSize', strategy1d.DefaultLegendFontSize);
    % set(leg, 'Interpreter', 'latex');
    % set(ax, 'YScale', 'log');
    % hold(ax, 'off');
    % xlim([0, config.tFinal]);
    % xl = xlabel('t');
    % ax.XAxis.Exponent = 0;
    % ax.XAxis.TickLabelFormat = '%.0f';
    % ylim([5e-8, 5e2]);
    % yticks(logspace(-16, 4, 21));
    % yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-16, 4, 21), 'Un', 0));
    % set(ax, 'FontSize', strategy2d.DefaultFontSize);
    % set(xl, 'FontSize', strategy2d.DefaultAxisLabelFontSize);
    % ti = title('Distance to macroscopic equilibrium');
    % set(ti, 'FontSize', strategy2d.DefaultTitleFontSize);
    % fileName = sprintf('%s_macro_relative_entropy.pdf', config.ckptPrefix);
    % fileName = fullfile(config.ckptDir, fileName);
    % exportgraphics(ax, fileName, 'ContentType', 'vector');
    % 
    % % Save trend to local equilibrium plot
    % figure(6);
    % ax = gca;
    % hold(ax, 'on');
    % t = state.History.time;
    % L2Entropy = max(1e-35, state.History.L2Entropy);
    % if k <= k0
    %     style = strategy1d.getDefaultLineStyle(k+1, k);
    % else
    %     style = strategy1d.getDefaultScatterStyle(k+1, k-k0);
    % end
    % line = plot(t, L2Entropy, style{:});
    % set(line, 'LineWidth', 3, 'MarkerSize', 8);
    % peLegs{k} = sprintf('$\\tau_0=10^{%d}$', tau0Level(k));
    % leg = legend(ax, peLegs(1:k));
    % set(leg, 'Location', strategy1d.DefaultLegendPosition);
    % set(leg, 'FontSize', strategy1d.DefaultLegendFontSize);
    % set(leg, 'Interpreter', 'latex');
    % set(ax, 'YScale', 'log');
    % hold(ax, 'off');
    % xlim([0, config.tFinal]);
    % xl = xlabel('t');
    % ax.XAxis.Exponent = 0;
    % ax.XAxis.TickLabelFormat = '%.0f';
    % ylim([5e-5, 5e2]);
    % yticks(logspace(-16, 4, 21));
    % yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-16, 4, 21), 'Un', 0));
    % set(ax, 'FontSize', strategy2d.DefaultFontSize);
    % set(xl, 'FontSize', strategy2d.DefaultAxisLabelFontSize);
    % ti = title('Distance to local equilibrium');
    % set(ti, 'FontSize', strategy2d.DefaultTitleFontSize);
    % fileName = sprintf('%s_micro_relative_entropy.pdf', config.ckptPrefix);
    % fileName = fullfile(config.ckptDir, fileName);
    % exportgraphics(ax, fileName, 'ContentType', 'vector');
end


% Save conservation plot
figure(8);
ax = gca;
hold(ax, 'on');
style1 = strategy1d.getDefaultLineStyle(1, 1);
style2 = strategy1d.getDefaultLineStyle(2, 2);
style3 = strategy1d.getDefaultLineStyle(3, 3);
t = state.History.time;
mass = state.History.mass;
momentum = state.History.momentum;
totalEnergy = state.History.totalEnergy;
mass = mass - mass(1);
momentum = momentum - momentum(1);
totalEnergy = totalEnergy - totalEnergy(1);
line = plot(t, mass, style1{:});
set(line, 'LineWidth', 2);
line = plot(t, momentum, style2{:});
set(line, 'LineWidth', 2);
line = plot(t, totalEnergy, style3{:});
set(line, 'LineWidth', 2);
distLegs = cell(1, 3);
distLegs{1} = 'mass';
distLegs{2} = 'momentum';
distLegs{3} = 'total energy';
leg = legend(ax, distLegs);
set(leg, 'Location', strategy1d.DefaultLegendPosition);
set(leg, 'FontSize', strategy1d.DefaultLegendFontSize);
set(leg, 'Interpreter', 'latex');
xlim([0, config.tFinal]);
xl = xlabel('$t$');
set(xl, 'Interpreter', 'latex');
ylim([-1e-5, 1e-5]);
yticks(-1e-5:5e-6:1e-5);
yticklabels(arrayfun(@(x) sprintf('%.0e', x), -1e-5:5e-6:1e-5, 'Un', 0));
title('Conservation');
fileName = sprintf('%s_conservation.pdf', config.ckptPrefix);
fileName = fullfile(config.ckptDir, fileName);
exportgraphics(ax, fileName, 'ContentType', 'vector');

%% SIMULATION

function [config, scheme, state] = run(config)

config.nDims = 1; % Dimension
switch config.id
    case 'ex03ldh'
        config.L = 4*pi; % Length of domain
        config.kappa = 2*pi / config.L; % Wave length
        config.xBBox = [0, config.L]; % Spatial bounding box
        config.delta = 0.01; % Perturbation
    case 'ex03ldi'
        config.L = 6; % Length of domain
        config.kappa = pi / config.L; % Wave length
        config.xBBox = [-config.L, config.L]; % Spatial bounding box
        config.delta = 0.01; % Perturbation
    case 'ex03tsh'
        config.L = 4*pi; % Length of domain
        config.kappa = 2*pi / config.L; % Wave length
        config.xBBox = [0, config.L]; % Spatial bounding box
        config.delta = 0.01; % Perturbation
        config.alpha1 = 1 / 6;
        config.alpha2 = 5 / 6;
    case 'ex03tsi'
        config.L = 6; % Length of domain
        config.kappa = pi / config.L; % Wave length
        config.xBBox = [-config.L, config.L]; % Spatial bounding box
        config.T0 = 1; % Temperature parameter
        config.delta = 0.01; % Perturbation
        config.alpha1 = 1 / 6;
        config.alpha2 = 5 / 6;
end
config.T0 = 1; % Temperature parameter
config.nx = 128; % Grid resolution
config.nh = 50; % Number of Hermite modes - 1
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
switch config.id
    case {'ex03ldh', 'ex03ldi'}
        D0 = (rhoInf + delta * cos(kappa*x)) ./ sqrt(rhoInf);
        D = [D0, zeros(1, size(x, 2)*(nh-1))];
    case 'ex03tsh'
        alpha1 = config.alpha1;
        alpha2 = config.alpha2;
        tmp = (cos(2*kappa*x) + cos(3*kappa*x)) * alpha2 + cos(kappa*x);
        D0 = (alpha1 + alpha2) * (1 + delta * tmp) ;
        D2 = sqrt(2) * alpha2 * (1 + delta * tmp);
        if nh > 2
            D = [D0, zeros(1, size(x, 2)), D2, zeros(1, size(x,2)*(nh - 3))];
        else
            error('nh must be larger than 2.');
        end
    case 'ex03tsi'
        alpha1 = config.alpha1;
        alpha2 = config.alpha2;
        tmp = cos(kappa*x);
        D0 = (alpha1 + alpha2) * (1 + delta * tmp) ./ sqrt(rhoInf);
        D2 = sqrt(2) * alpha2 * (1 + delta * tmp) ./ sqrt(rhoInf);
        if nh > 2
            D = [D0, zeros(1, size(x, 2)), D2, zeros(1, size(x,2)*(nh-3))];
        else
            error('nh must be larger than 2.');
        end
end
end

function h = phiInf(config)
h = @(x) phiInfImpl(config, x);
end

function phi = phiInfImpl(config, x)
switch config.id
    case {'ex03ldh', 'ex03tsh'}
        phi = zeros(1, size(x, 2));
    case 'ex03ldi'
        kappa = config.kappa;
        phi = 0.2*sin(kappa*x);
    case 'ex03tsi'
        kappa = config.kappa;
        phi = 0.1*(1-cos(kappa*x));
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