%==========================================================================
% FileName: ex05_weak_inhomogeneous_landau_damping
% Author: Yi CAI
% Date: 13/09/2025
% Description:
%   Weak inhomogeneous Landau damping with Vlasov solver.
%==========================================================================

clc, clear, close all;

%% EXECUTION
dataPath = fullfile(fileparts(mfilename('fullpath')), 'ex05');

strategy1d = physics.visual.Strategy1d();
legends = cell(1, 4);
for k = 0:3
    epsilon = 10^(-k);
    [scheme, state] = run(epsilon, 10^5, 100);
    figure(3);
    ax = gca;
    hold(ax, 'on');
    style = strategy1d.getDefaultLineStyle(k+1, k+1);
    t = state.history.t;
    potential = max(1e-25, state.history.potential);
    line = plot(t, potential, style{:});
    set(line, 'LineWidth', 2);
    set(ax, 'YScale', 'log');
    hold(ax, 'off');
    legends{k+1} = sprintf('$\\varepsilon=10^%d$', k);
end

%% PLOTTING
leg = legend(ax, legends);
set(leg, 'Location', strategy1d.DefaultLegendPosition);
set(leg, 'FontSize', strategy1d.DefaultLegendFontSize);
set(leg, 'Interpreter', 'latex');
ax = gca;
delete(ax.Children(2:end))
fileName = sprintf('figure_LDHSM%d.pdf', scheme.getConfig('xBasisOrder'));
exportgraphics(ax, fullfile(dataPath, fileName), 'ContentType', 'vector');

%% SIMULATION

function [scheme, state] = run(epsilon, tau0, tFinal)

% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');

% Experiment-specific configuration struct
config = struct();
config.schemeName = 'dghs';
config.nDims = 1; % Dimension
config.L = 6; % Length
config.xBBox = [-config.L, config.L]; % Spatial bounding box
config.tFinal = tFinal; % Final time
config.k = pi / config.L; % Wave length
config.T0 = 1; % Temperature parameter
config.delta = 0.01; % Perturbation
config.phiInf = phiInf(config.k); % Potential steady state
config.EInf = EInf(config.k); % Electrical field steady state
config.cInf = 2 * config.L / integral(@(x) exp(-config.phiInf(x)), -config.L, config.L); % Normalization parameter
config.rhoInf = rhoInf(config.k, config.cInf); % Density steady state
config.sqrtRhoInf = @(x) sqrt(config.rhoInf(x));
config.invRhoInf = @(x) 1 ./ config.rhoInf(x);
config.invSqrtRhoInf = @(x) 1 ./ sqrt(config.rhoInf(x));
config.nx = repmat(64, 1, config.nDims); % Grid resolution
config.nh = 80; % Number of Hermite modes - 1
config.ic.D = DInit(config.k, config.delta, config.cInf, config.nh); % Initial hermite coefficients
config.ic.omega = omegaInit(config.k, config.delta, config.cInf); % Initial potential
config.bc = []; % Boundary condition
config.exact = [];
config.tau0 = tau0; %
config.epsilon = epsilon; % Scaling parameter
config.tau = config.tau0;
config.verbose = 1; % Verbose flag
config.eId = sprintf('AP%d', config.nDims); % Experiment ID

% Construct and configure scheme
switch config.schemeName
    case 'dghs'
        scheme = physics.vlasov.DghsScheme(file = common, config = config);
end

% Configure visualizer
visualizer = scheme.getConfig('visualizer');
visualizer.setDensity(ones(1, config.nDims));
visualizer.setTimeline(scheme.getConfig('tFinal'), 10);
visualizer.setComponents(struct('D', [], 'omega', []));
visualizer.addDataset(scheme.getConfig('sId'), 2);
visualizer.addPlotter('profile', '1d');

xElement = approx.element.BH1OrthotopeElement(config.nDims, ...
    config.xBasisOrder, config.xBasisType, config.xBasisPattern);
xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx);
xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
state = physics.vlasov.HermiteState(xDisc, config.nh);
state = scheme.run(state);
end

%% HELPER FUNCTIONS

function h = fInit(k, delta)
h = @(x, v) fInitImpl(x, v, k, delta);
end

function f = fInitImpl(x, v, k, delta)
M = exp(-abs(v).^2/(2)) / sqrt(2*pi);
f = (rhoInf(x) + delta * cos(k*x)) .* M;
end

function h = DInit(k, delta, cInf, nh)
h = @(x) DInitImpl(x, k, delta, cInf, nh);
end

function D = DInitImpl(x, k, delta, cInf, nh)
rhoInf = rhoInfImpl(x, k, cInf);
D0 = (rhoInf + delta * cos(k*x)) ./ rhoInf;
D = [D0, zeros(1, size(x, 2)*nh)];
end

function h = omegaInit(k, delta, cInf)
h = @(x) omegaInitImpl(x, k, delta, cInf);
end

function omega = omegaInitImpl(x, k, delta, cInf)
rhoInf = rhoInfImpl(x, k, cInf);
omega = delta / k^2 * cos(k*x) .* sqrt(rhoInf);
end

function h = phiInf(k)
h = @(x) phiInfImpl(x, k);
end

function phi = phiInfImpl(x, k)
phi = 0.2*sin(k*x);
end

function h = EInf(k)
h = @(x) EInfImpl(x, k);
end

function E = EInfImpl(x, k)
E = -0.2*k*cos(k*x);
end

function h = rhoInf(k, cInf)
h = @(x) rhoInfImpl(x, k, cInf);
end

function rho = rhoInfImpl(x, k, cInf)
rho = cInf * exp(-phiInfImpl(x, k));
end