%==========================================================================
% FileName: ex01h_convergence
% Author: Yi CAI
% Date: 13/09/2025
% Description:
%   Homogeneous Landau damping convergence test.
%==========================================================================

clc, clear, close all;

%% EXECUTION
% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = physics.scheme.ConfigParser.parseFile(common);

% Checkpoint information
config.ckptDir = fullfile(fileparts(mfilename('fullpath')), 'ex05h');
config.ckptPrefix = [upper(config.schemeName), num2str(config.xBasisOrder)];
if config.useFilter
    config.ckptPrefix = ['F', config.ckptPrefix];
end
config.ckptPrefix = strjoin({'LDH', config.ckptPrefix}, '_');
config.tau0 = 1;
config.ckptFigureName = 'f';

override = 1;
strategy1d = physics.visual.Strategy1d();
strategy2d = physics.visual.Strategy2d();
epsilonLevel = [-1, -2];
epsilon = 10.^epsilonLevel;
peLegs = cell(1, length(epsilon));
for k = 1:length(epsilon)
    config.epsilon = epsilon(k);
    config.tFinal = 100 * epsilon(k);
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
    V = linspace(-6, 6, 512);
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
    style = strategy1d.getDefaultLineStyle(k+1, k+1);
    t = state.History.time / epsilon(k);
    potentialEnergy = max(1e-16, state.History.potentialEnergy);
    line = plot(t, potentialEnergy, style{:});
    set(line, 'LineWidth', 2);
    set(ax, 'YScale', 'log');
    hold(ax, 'off');
    peLegs{k} = sprintf('$\\varepsilon$=%.0e', epsilon(k));
    leg = legend(ax, peLegs(1:k));
    set(leg, 'Location', strategy1d.DefaultLegendPosition);
    set(leg, 'FontSize', strategy1d.DefaultLegendFontSize);
    set(leg, 'Interpreter', 'latex');
    xl = xlabel('$\varepsilon t$');
    xlim([0, config.tFinal/epsilon(k)]);
    set(xl, 'Interpreter', 'latex');
    ylim([1e-16, 1e0]);
    yticks(logspace(-16, 0, 9));
    yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-16, 0, 9), 'Un', 0));
    title('Potential energy');
    fileName = sprintf('%s_potential_energy.pdf', config.ckptPrefix);
    fileName = fullfile(config.ckptDir, fileName);
    exportgraphics(ax, fileName, 'ContentType', 'vector');

    % Save trend to equilibrium plot
    figure(4);
    ax = gca;
    hold(ax, 'on');
    style = strategy1d.getDefaultLineStyle(k+1, k+1);
    t = state.History.time / epsilon(k);
    L2FDistance = max(1e-25, state.History.L2FDistance);
    line = plot(t, L2FDistance, style{:});
    set(line, 'LineWidth', 2);
    peLegs{k} = sprintf('$\\varepsilon$=%.0e', epsilon(k));
    leg = legend(ax, peLegs(1:k));
    set(leg, 'Location', strategy1d.DefaultLegendPosition);
    set(leg, 'FontSize', strategy1d.DefaultLegendFontSize);
    set(leg, 'Interpreter', 'latex');
    set(ax, 'YScale', 'log');
    hold(ax, 'off');
    xlim([0, config.tFinal/epsilon(k)]);
    xl = xlabel('$\varepsilon t$');
    set(xl, 'Interpreter', 'latex');
    ylim([1e-16, 1e0]);
    yticks(logspace(-16, 0, 9));
    yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-16, 0, 9), 'Un', 0));
    title('Trend to equilibrium');
    fileName = sprintf('%s_trend2eq.pdf', config.ckptPrefix);
    fileName = fullfile(config.ckptDir, fileName);
    exportgraphics(ax, fileName, 'ContentType', 'vector');
end



%% SIMULATION

function [config, scheme, state] = run(config)

config.nDims = 1; % Dimension
config.kappa = 0.5; % Wave length
config.L = 2 * pi / config.kappa; % Length of domain
config.xBBox = [0, config.L]; % Spatial bounding box
config.vBBox = [-6, 6]; % Velocity bounding box
config.T0 = 1; % Temperature parameter
config.delta = 0.05; % Perturbation
config.nx = 64; % Grid resolution
config.nh = 128; % Number of Hermite modes - 1
config.ic = DInit(config.kappa, config.delta); % Initial hermite coefficients
config.bc = []; % Boundary condition
config.exact = []; % Exact solution
config.verbose = 1; % Verbose flag

% Construct and configure scheme
switch config.schemeName
    case 'dghs'
        scheme = physics.vlasov.DghsScheme(config = config);
    case 'cdghs'
        scheme = physics.vlasov.CdghsScheme(config = config);
end

% Configure visualizer
visualizer = scheme.getConfig('visualizer');
visualizer.setDensity(1);
visualizer.setTimeline(scheme.getConfig('tFinal'), 1);
visualizer.setComponents(struct('D', [], 'P', [], 'E', []));
visualizer.addDataset(scheme.getConfig('sId'), 'numeric');
visualizer.addPlotter('profile', '1d');

% Configure analyzer
analyzer = scheme.getConfig('analyzer');
analyzer.setNLevels(5);
analyzer.setDensity(4^config.nDims);
analyzer.setComponents(struct('D', 1:config.nh, 'P', 1, 'E', 1));

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
state = physics.vlasov.HermiteState(xDisc, vDisc);
state = scheme.run(state);
end

%% HELPER FUNCTIONS

function h = DInit(kappa, delta)
h = @(x, varargin) DInitImpl(kappa, delta, x, varargin{:});
end

function D = DInitImpl(kappa, delta, x, varargin)
if nargin >= 1, nh = varargin{1}; end
D0 = (1 + delta * cos(kappa*x));
D = [D0, zeros(1, size(x, 2)*(nh-1))];
end