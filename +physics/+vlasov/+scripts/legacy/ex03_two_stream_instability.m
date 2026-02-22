%==========================================================================
% FileName: ex03h_two_stream_instability
% Author: Yi CAI
% Date: 13/09/2025
% Description:
%   Two-stream instability with Hermite spectral method.
%==========================================================================

clc, clear, close all;

%% EXECUTION
% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = physics.scheme.ConfigParser.parseFile(common);

% Checkpoint information
config.ckptDir = fullfile(fileparts(mfilename('fullpath')), 'ex03');
config.ckptPrefix = [upper(config.schemeName), num2str(config.xBasisOrder)];
if config.useFilter
    config.ckptPrefix = ['F', config.ckptPrefix];
end
config.ckptPrefix = strjoin({'TSH', config.ckptPrefix}, '_');
config.tFinal = 8;
config.ckptTimeStamps = [];
config.ckptFigureName = 'f';

override = false;
strategy1d = physics.visual.Strategy1d();
strategy2d = physics.visual.Strategy2d();
nh = 50;
peLegs = cell(1, length(nh));
for k = 1:length(nh)
    config.epsilon = 1;
    config.tau0 = 1;
    config.nh = nh(k);
    config.ckptPostfix = sprintf('nh%d', config.nh);
    fileName = sprintf('%s_state_%s.mat', config.ckptPrefix, config.ckptPostfix);
%     fileName = sprintf('%s_state_%s_t20p0.mat', config.ckptPrefix, config.ckptPostfix);
%     fileName = sprintf('%s_state_%s_t30p0.mat', config.ckptPrefix, config.ckptPostfix);
%     fileName = sprintf('%s_state_%s_t40p0.mat', config.ckptPrefix, config.ckptPostfix);
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
    caxis([0, 0.3]);
    fileName = sprintf('%s_f_%s.pdf', config.ckptPrefix, config.ckptPostfix);

%     fileName = sprintf('%s_f_%s_t20p0.pdf', config.ckptPrefix, config.ckptPostfix);

%         fileName = sprintf('%s_f_%s_t30p0.pdf', config.ckptPrefix, config.ckptPostfix);

%         fileName = sprintf('%s_f_%s_t40p0.pdf', config.ckptPrefix, config.ckptPostfix);
    fileName = fullfile(config.ckptDir, fileName);
    exportgraphics(ax, fileName, 'ContentType', 'vector');

    % Save potential energy plot
    figure(3);
    ax = gca;
    hold(ax, 'on');
    style = strategy1d.getDefaultLineStyle(k+1, k+1);
    t = state.History.time;
    potentialEnergy = max(1e-16, state.History.potentialEnergy);
    line = plot(t, potentialEnergy, style{:});
    set(line, 'LineWidth', 2);
    set(ax, 'YScale', 'log');
    hold(ax, 'off');
    peLegs{k} = sprintf('$N_H=%d$', config.nh);
    leg = legend(ax, peLegs(1:k));
    set(leg, 'Location', strategy1d.DefaultLegendPosition);
    set(leg, 'FontSize', strategy1d.DefaultLegendFontSize);
    set(leg, 'Interpreter', 'latex');
    xl = xlabel('$t$');
    xlim([0, config.tFinal]);
    set(xl, 'Interpreter', 'latex');
    ylim([1e-8, 1e2]);
    yticks(logspace(-8, 2, 6));
    yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-8, 2, 6), 'Un', 0));
    title('Potential energy');
    fileName = sprintf('%s_potential_energy.pdf', config.ckptPrefix);
    fileName = fullfile(config.ckptDir, fileName);
    exportgraphics(ax, fileName, 'ContentType', 'vector');
end

% Save trend to equilibrium plot
% figure(4);
% ax = gca;
% hold(ax, 'on');
% style1 = strategy1d.getDefaultLineStyle(1, 1);
% style2 = strategy1d.getDefaultLineStyle(2, 2);
% style3 = strategy1d.getDefaultLineStyle(3, 3);
% t = state.History.time;
% L2Entropy = max(1e-25, state.History.L2Entropy);
% L2FDistance = max(1e-25, state.History.L2FDistance);
% L2RhoDistance = max(1e-25, state.History.L2RhoDistance);
% line = plot(t, L2Entropy, style1{:});
% set(line, 'LineWidth', 2);
% line = plot(t, L2FDistance, style2{:});
% set(line, 'LineWidth', 2);
% line = plot(t, L2RhoDistance, style3{:});
% set(line, 'LineWidth', 2);
% distLegs = cell(1, 3);
% distLegs{1} = '$\|f-\rho M\|$';
% distLegs{2} = '$\|f-f_\infty\|$';
% distLegs{3} = '$\|\rho-\rho_\infty\|$';
% leg = legend(ax, distLegs);
% set(leg, 'Location', strategy1d.DefaultLegendPosition);
% set(leg, 'FontSize', strategy1d.DefaultLegendFontSize);
% set(leg, 'Interpreter', 'latex');
% set(ax, 'YScale', 'log');
% hold(ax, 'off');
% xlim([0, config.tFinal]);
% xl = xlabel('$t$');
% set(xl, 'Interpreter', 'latex');
% ylim([1e-8, 1e0]);
% yticks(logspace(-8, 0, 5));
% yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-8, 0, 5), 'Un', 0));
% title('Trend to equilibrium');
% fileName = sprintf('%s_trend2eq.pdf', config.ckptPrefix);
% fileName = fullfile(config.ckptDir, fileName);
% exportgraphics(ax, fileName, 'ContentType', 'vector');

% Save conservation plot
figure(5);
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
% ylim([-1e-5, 1e-5]);
% yticks(-1e-5:5e-6:1e-5);
% yticklabels(arrayfun(@(x) sprintf('%.0e', x), -1e-5:5e-6:1e-5, 'Un', 0));
title('Conservation');
fileName = sprintf('%s_conservation.pdf', config.ckptPrefix);
fileName = fullfile(config.ckptDir, fileName);
exportgraphics(ax, fileName, 'ContentType', 'vector');

%% SIMULATION

function [config, scheme, state] = run(config)

% Experiment-specific configuration struct
config.nDims = 1; % Dimension
config.L = 6; % Length of domain
config.kappa = pi / config.L; % Wave length
config.xBBox = [0, config.L]; % Spatial bounding box
config.vBBox = [-3, 3]; % Velocity bounding box
config.T0 = 1; % Temperature parameter
config.delta = 0.01; % Perturbation
config.alpha1 = 1 / 6;
config.alpha2 = 5 / 6;
config.nx = 64; % Grid resolution
config.phiInf = phiInf(config.kappa); % Electric potential equilibrium
config.cInf = cInf(config.L, config.phiInf, config.T0); % Normalization constant
config.rhoInf = rhoInf(config.kappa, config.cInf, config.T0); % Density equilibrium
config.sqrtRhoInf = sqrtRhoInf(config.kappa, config.cInf, config.T0); % Square root of density equilibrium
config.ic = DInit(config.kappa, config.delta, config.cInf, config.T0); % Initial hermite coefficients
config.bc = []; % Boundary condition
config.exact = []; % Exact solution
config.verbose = 2; % Verbose flag

% Construct and configure scheme
switch config.schemeName
    case 'dghs'
        scheme = physics.vlasov.DghsScheme(config = config);
end

% Configure visualizer
visualizer = scheme.getConfig('visualizer');
visualizer.setDensity(1);
visualizer.setTimeline(scheme.getConfig('tFinal'), 1);
visualizer.setComponents(struct('D', [], 'P', [], 'E', []));
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
xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx);
xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
state = physics.vlasov.HermiteState(xDisc, vDisc);
state = scheme.run(state);
end

%% HELPER FUNCTIONS

function h = DInit(kappa, delta, alpha1, alpha2, cInf, T0)
h = @(x, varargin) DInitImpl(kappa, delta, alpha1, alpha2, cInf, T0, x, varargin{:});
end

function D = DInitImpl(kappa, delta, alpha1, alpha2, cInf, T0, x, varargin)
if nargin >= 1, nh = varargin{1}; end

rhoInf = rhoInfImpl(kappa, cInf, T0, x);

D0 = (alpha1 + alpha2) * (rhoInf + delta * cos(kappa*x)) ./ sqrt(rhoInf);

if nh == 1
    D = D0;
elseif nh == 2
    D = [D0, zeros(1, size(x, 2))];
else
    D2 = sqrt(2) * alpha2 * (rhoInf + delta * cos(kappa*x)) ./ sqrt(rhoInf);
    D = [D0, zeros(1, size(x, 2)), D2, zeros(1, size(x,2)*(nh-3))];
end
end

function h = phiInf(kappa)
h = @(x) phiInfImpl(kappa, x);
end

function phi = phiInfImpl(kappa, x)
phi = 0.1 * (1 - cos(kappa*x));
end

function h = rhoInf(kappa, cInf, T0)
h = @(x) rhoInfImpl(kappa, cInf, T0, x);
end

function rho = rhoInfImpl(kappa, cInf, T0, x)
rho = cInf * exp(-phiInfImpl(kappa, x) / T0);
end

function h = sqrtRhoInf(kappa, cInf, T0)
h = @(x) sqrtRhoInfImpl(kappa, cInf, T0, x);
end

function rho = sqrtRhoInfImpl(kappa, cInf, T0, x)
rho = sqrt(rhoInfImpl(kappa, cInf, T0, x));
end

function c = cInf(L, phiInf, T0)
c = 2 * L / integral(@(x) exp(-phiInf(x)/T0), -L, L);
end