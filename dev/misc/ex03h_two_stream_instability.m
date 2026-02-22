%==========================================================================
% FileName: ex03h_two_stream_instability
% Author: Yi CAI
% Date: 13/09/2025
% Description:
%   Two-stream instability with Hermite-based Vlasov solver.
%==========================================================================

clc, clear, close all;

%% EXECUTION
dataPath = fullfile(fileparts(mfilename('fullpath')), 'ex03h');

strategy1d = physics.visual.Strategy1d();
tau0 = 10.^(2:4);
nh = 2;
for j = 1:length(tau0)
    stateName = sprintf('state_%d', j);
    figureName = sprintf('f_%d', j);
    [scheme, state] = run(1, tau0(j), 60, nh, dataPath, stateName, figureName);

    % Save distribution plot
    figure(1);
    ax = gca;
    fileName = sprintf('TS%d_%s.pdf', scheme.getConfig('xBasisOrder'), figureName);
    exportgraphics(ax, fullfile(dataPath, fileName), 'ContentType', 'vector');

    % Save state
    fileName = sprintf('TS%d_%s.mat', ...
        scheme.getConfig('xBasisOrder'), stateName);
    state.save(fullfile(dataPath, fileName));

    % Plot time series
    figure(2);
    ax = gca;
    hold(ax, 'on');
    style = strategy1d.getDefaultLineStyle(j, j);
    t = state.History.t;
    potential = max(1e-25, state.History.potential);
    line = plot(t, potential, style{:});
    set(line, 'LineWidth', 2);
    set(ax, 'YScale', 'log');
    hold(ax, 'off');
end

legends = cell(1, length(tau0));
for j = 1:length(tau0)
    legends{j-1} = sprintf('$\\tau_0=10^{%d}$', log10(tau0(j)));
end
leg = legend(ax, legends);
set(leg, 'Location', strategy1d.DefaultLegendPosition);
set(leg, 'FontSize', strategy1d.DefaultLegendFontSize);
set(leg, 'Interpreter', 'latex');
xl = xlabel('$t$');
yl = ylabel('$\|\partial_x \Psi\|_2$');
set(xl, 'Interpreter', 'latex');
set(yl, 'Interpreter', 'latex');
ylim([1e-4, 1e4]);
yticks(10.^(-4:2:4));
yticklabels(arrayfun(@(x) sprintf('10^{%d}', x), -4:2:4, 'Un', 0));
title('Potential energy');
ax = gca;
fileName = sprintf('TS%d_potential.pdf', scheme.getConfig('xBasisOrder'));
exportgraphics(ax, fullfile(dataPath, fileName), 'ContentType', 'vector');

%% SIMULATION

function [scheme, state] = run(epsilon, tau0, tFinal, nh, dataPath, stateName, figureName)

% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');

% Experiment-specific configuration struct
config = struct();
config.schemeName = 'dghs';
config.nDims = 1; % Dimension
config.k = 0.5; % Wave length
config.L = 2*pi / config.k; % Length of domain
config.xBBox = [0, config.L]; % Spatial bounding box
config.vBBox = [-6, 6]; % Velocity bounding box
config.tFinal = tFinal; % Final time
config.T0 = 1; % Temperature parameter
config.alpha1 = 1 / 6;
config.alpha2 = 5 / 6;
config.delta = 0.01; % Perturbation
config.phiInf = phiInf(config.k); % Potential steady state
config.EInf = EInf(config.k); % Electrical field steady state
config.cInf = config.L / integral(@(x) exp(-config.phiInf(x)), 0, config.L); % Normalization parameter
config.rhoInf = rhoInf(config.k, config.cInf); % Density steady state
config.nx = repmat(64, 1, config.nDims); % Grid resolution
config.nh = nh; % Number of Hermite modes - 1
config.ic.D = DInit(config.k, config.alpha1, config.alpha2, config.delta, config.cInf, config.nh); % Initial hermite coefficients
config.ic.omega = omegaInit(config.k, config.alpha1, config.alpha2, config.delta, config.cInf); % Initial potential
config.bc = []; % Boundary condition
config.exact = []; % Exact solution
config.tau0 = tau0; % Collision parameter
config.epsilon = epsilon; % Scaling parameter
config.verbose = 1; % Verbose flag
config.eId = sprintf('AP%d', config.nDims); % Experiment ID
config.ckptDir = dataPath;
config.ckptCounts = [8, 16, 32, 50] ./ config.dt;
config.ckptStateName = stateName;
config.ckptFigureName = figureName;

% Construct and configure scheme
switch config.schemeName
    case 'dghs'
        scheme = physics.vlasov.DghsScheme(file = common, config = config);
end

% Configure visualizer
visualizer = scheme.getConfig('visualizer');
visualizer.setDensity(1);
visualizer.setTimeline(scheme.getConfig('tFinal'), 1);
visualizer.setComponents(struct('D', [], 'omega', []));
visualizer.addDataset(scheme.getConfig('sId'), 'numeric');
visualizer.addPlotter('profile', '1d');

vElement = approx.element.C0OrthotopeElement(config.vBBox, ...
    config.nh+1, config.vBasisType, lambda = sqrt(config.T0));
vDisc = approx.space.SpectralSpace(vElement);

xElement = approx.element.BH1OrthotopeElement(config.nDims, ...
    config.xBasisOrder, config.xBasisType, config.xBasisPattern);
xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx);
xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
state = physics.vlasov.HermiteState(xDisc, vDisc, config.T0, config.rhoInf);
state = scheme.run(state);
end

%% HELPER FUNCTIONS

function h = DInit(k, alpha1, alpha2, delta, cInf, nh)
h = @(x) DInitImpl(k, alpha1, alpha2, delta, cInf, nh, x);
end

function D = DInitImpl(k, alpha1, alpha2, delta, cInf, nh, x)
rhoInf = rhoInfImpl(k, cInf, x);
O = zeros(1, size(x, 2));
if nh >= 0
    D0 = (alpha1 + alpha2) * (rhoInf + delta * cos(k*x)) ./ sqrt(rhoInf);
end

if nh >= 2
    D2 = sqrt(2) * alpha2 * (rhoInf + delta * cos(k*x)) ./ sqrt(rhoInf);
end

if nh < 2
    D = [D0, repmat(O, 1, nh)];
else
    D = [D0, O, D2, repmat(O, 1, nh-2)];
end
end

function h = omegaInit(k, alpha1, alpha2, delta, cInf)
h = @(x) omegaInitImpl(k, alpha1, alpha2, delta, cInf, x);
end

function omega = omegaInitImpl(k, alpha1, alpha2, delta, cInf, x)
rhoInf = rhoInfImpl(k, cInf, x);
omega = (alpha1 + alpha2) * delta / k^2 * cos(k*x) .* sqrt(rhoInf);
end

function h = phiInf(k)
h = @(x) phiInfImpl(k, x);
end

function phi = phiInfImpl(~, x)
phi = zeros(size(x));
end

function h = EInf(k)
h = @(x) EInfImpl(k, x);
end

function E = EInfImpl(~, x)
E = zeros(size(x));
end

function h = rhoInf(k, cInf)
h = @(x) rhoInfImpl(k, cInf, x);
end

function rho = rhoInfImpl(k, cInf, x)
rho = cInf * exp(-phiInfImpl(k, x));
end