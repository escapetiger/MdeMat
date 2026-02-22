%==========================================================================
% FileName: ex04_asymptotic_preserving
% Author: Yi CAI
% Date: 13/09/2025
% Description:
%   Asymptotic preserving Vlasov scheme test.
%==========================================================================

clc, clear, close all;

%% EXECUTION
dataPath = fullfile(fileparts(mfilename('fullpath')), 'ex04');

% tFinal = [0.2, 0.8, 1.6];
% for j = 3
%     for k = 1
%         [config, state] = simulate(10^(-k), tFinal(j));
%         fileName = sprintf('%s_%d.mat', config.SID, k);
%         state.save(fullfile(dataPath, fileName));
%     end
%
%     ax = gca;
%     delete(ax.Children(2:end))
%     fileName = sprintf('figure_LDHS%d_%d.pdf', config.xBasisOrder, j);
%     exportgraphics(ax, fullfile(dataPath, fileName), 'ContentType', 'vector');
% end

%% PLOTTING

strategy1d = physics.visual.Strategy1d();
for k = 0:3
    epsilon = 10^(-k);
    [scheme, state] = run(epsilon, 1, 40 * epsilon);
    figure(3);
    ax = gca;
    hold(ax, 'on');
    style = strategy1d.getDefaultLineStyle(k+1, k+1);
    t = state.History.t / epsilon;
    potential = max(1e-25, state.History.potential);
    line = plot(t, potential, style{:});
    set(line, 'LineWidth', 2);
    set(ax, 'YScale', 'log');
    hold(ax, 'off');
end

legends = cell(1, 4);
for k = 0:3
    legends{k+1} = sprintf('$\\varepsilon=10^{-%d}$', k);
end
leg = legend(ax, legends);
set(leg, 'Location', strategy1d.DefaultLegendPosition);
set(leg, 'FontSize', strategy1d.DefaultLegendFontSize);
set(leg, 'Interpreter', 'latex');
xl = xlabel('$\varepsilon t$');
yl = ylabel('$\|\partial_x \Psi\|_2$');
set(xl, 'Interpreter', 'latex');
set(yl, 'Interpreter', 'latex');
title('Potential energy in log scale');
ax = gca;
fileName = sprintf('figure_LDHWAP%d_scaled.pdf', scheme.getConfig('xBasisOrder'));
exportgraphics(ax, fullfile(dataPath, fileName), 'ContentType', 'vector');

%% Original timeline
% strategy1d = physics.visual.Strategy1d();
% for k = 0:3
%     epsilon = 10^(-k);
%     [config, state] = simulate(epsilon, 1, 40);
%     figure(3);
%     ax = gca;
%     hold(ax, 'on');
%     style = strategy1d.getDefaultLineStyle(k+1, k+1);
%     t = state.History.t;
%     potential = max(1e-25, state.History.potential);
%     line = plot(t, potential, style{:});
%     set(line, 'LineWidth', 2);
%     set(ax, 'YScale', 'log');
%     hold(ax, 'off');
% end
% 
% legends = cell(1, 4);
% for k = 0:3
%     legends{k+1} = sprintf('$\\varepsilon=10^{-%d}$', k);
% end
% leg = legend(ax, legends);
% set(leg, 'Location', strategy1d.DEFAULT_LEGEND_POSITION);
% set(leg, 'FontSize', strategy1d.DEFAULT_LEGEND_FONT_SIZE);
% set(leg, 'Interpreter', 'latex');
% xl = xlabel('$\varepsilon t$');
% yl = ylabel('$\|\partial_x \Psi\|_2$');
% set(xl, 'Interpreter', 'latex');
% set(yl, 'Interpreter', 'latex');
% title('Potential energy in log scale');
% ax = gca;
% fileName = sprintf('figure_LDHWAP%d_original.pdf', config.xBasisOrder);
% exportgraphics(ax, fullfile(dataPath, fileName), 'ContentType', 'vector');

%% SIMULATION

function [scheme, state] = run(epsilon, tau0, tFinal)

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
config.delta = 0.01; % Perturbation
config.phiInf = phiInf(config.k); % Potential steady state
config.EInf = EInf(config.k); % Electrical field steady state
config.cInf = config.L / integral(@(x) exp(-config.phiInf(x)), 0, config.L); % Normalization parameter
config.rhoInf = rhoInf(config.k, config.cInf); % Density steady state
config.nx = repmat(64, 1, config.nDims); % Grid resolution
config.nh = 80; % Number of Hermite modes - 1
config.ic.D = DInit(config.k, config.delta, config.cInf, config.nh); % Initial hermite coefficients
config.ic.omega = omegaInit(config.k, config.delta, config.cInf); % Initial potential
config.bc = []; % Boundary condition
config.exact = []; % Exact solution
config.tau0 = tau0; % Collision parameter
config.epsilon = epsilon; % Scaling parameter
config.verbose = 1; % Verbose flag

config.eId = sprintf('AP%d', config.nDims); % Experiment ID

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
    config.nh + 1, config.vBasisType, lambda=sqrt(config.T0));
vDisc = approx.space.SpectralSpace(vElement);

xElement = approx.element.BH1OrthotopeElement(config.nDims, ...
    config.xBasisOrder, config.xBasisType, config.xBasisPattern);
xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx);
xDisc = approx.space.FiniteElementSpace(xMesh, xElement);

state = physics.vlasov.HermiteState(xDisc, vDisc, config.T0, config.rhoInf);
state = scheme.run(state);
end

%% HELPER FUNCTIONS

function h = DInit(k, delta, cInf, nh)
h = @(x) DInitImpl(k, delta, cInf, nh, x);
end

function D = DInitImpl(k, delta, cInf, nh, x)
rhoInf = rhoInfImpl(k, cInf, x);
D0 = (rhoInf + delta * cos(k*x)) ./ rhoInf;
D = [D0, zeros(1, size(x, 2)*nh)];
end

function h = omegaInit(k, delta, cInf)
h = @(x) omegaInitImpl(k, delta, cInf, x);
end

function omega = omegaInitImpl(k, delta, cInf, x)
rhoInf = rhoInfImpl(k, cInf, x);
omega = delta / k^2 * cos(k*x) .* sqrt(rhoInf);
end

function h = phiInf(k)
h = @(x) phiInfImpl(k, x);
end

function phi = phiInfImpl(k, x)
% phi = -0.2*sin(k*x);
phi = zeros(size(x));
end

function h = EInf(k)
h = @(x) EInfImpl(k, x);
end

function E = EInfImpl(k, x)
% E = 0.2*k*cos(k*x);
E = zeros(size(x));
end

function h = rhoInf(k, cInf)
h = @(x) rhoInfImpl(k, cInf, x);
end

function rho = rhoInfImpl(k, cInf, x)
rho = cInf * exp(-phiInfImpl(k, x));
end