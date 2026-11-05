%==========================================================================
% FileName: ex01_periodic
% Author: Yi Cai
% Date: 2025-07-13
% Description:
%   Solving advection equation with periodic boundary conditions.
%==========================================================================

clc;
clear;
close all;

%% CONFIGURATION
dataDir = fullfile(fileparts(mfilename('fullpath')), 'ex03');
configPath = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = core.io.OptionParser().parse(configPath);
config.nDims = 1; % Dimension
config.xBBox = repmat([0, 1], 1, config.nDims); % Spatial bounding box
config.tFinal = 1.1; % Final time
config.advection = ones(1, config.nDims); % Advection speed
config.omega = repmat(2*pi, 1, config.nDims); % Wave length
config.ic = u0(config.omega); % Initial condition
config.bc = []; % Boundary condition
config.exact = uExact(config.advection, config.omega); % Exact solution
config.nx = repmat(256, 1, config.nDims); % Grid resolution
config.experimentId = sprintf('D%d', config.nDims); % Experiment ID
config.tDiscId = upper(config.tDiscId);
config.xDiscId = sprintf('%s%d', upper(config.xDiscName), config.xBasisOrder);
config.schemeId = strjoin({config.tDiscId, config.xDiscId}, '-');
config.schemeName = strrep(config.schemeId, '-', '_');
config.saveAsRef = 0;
config.verbose = 1; % Verbose flag: 0 - none; 1 - runtime render

%% SIMULATION SETUP

DEFAULT_MARKERS = {'o', 's', '^', 'v', 'd', 'x', '+'};
DEFAULT_COLORS = {'#000000', '#004488', '#BB5566', '#DDAA33'};

refFile = sprintf('%s-REF.mat', config.experimentId);
refPath = fullfile(dataDir, refFile);
if exist(refPath, 'file')
    refDataset = load(refPath).REF;
else
    refDataset = [];
end

if contains(config.xDiscId, 'sl')
    cls = 'SldgScheme';
else
    cls = 'UwdgScheme';
end
scheme = physics.advection.(cls)(config) ...
    .setTimeDiscretization(config.tDiscId, config.tFinal) ...
    .setTimer(config.verbose) ...
    .setVisualizer( ...
    'experimentId', config.experimentId, ...
    'schemeId', config.schemeId, ...
    'nDims', config.nDims, ...
    'nTimeNodes', 10, ...
    'final', config.tFinal, ...
    'density', ones(1, config.nDims), ...
    'components', struct('U', 1)) ...
    .addDataset('REF', 1) ...
    .addStyle('REF', ...
    {'Color', DEFAULT_COLORS{1}, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', 1}) ...
    .addDataset(config.schemeName, 2) ...
    .addStyle(config.schemeName, ...
    {'Color', DEFAULT_COLORS{2}, 'Marker', DEFAULT_MARKERS{1}, 'LineStyle', 'none', 'LineWidth', 1});

%% SIMULATION EXECUTION
if contains(config.xDiscId, 'sl')
    xDisc = approx.discretization.SemiLagrangianFiniteElementDiscretization() ...
        .setDgOrthotopeElement( ...
        'nDims', config.nDims, ...
        'order', config.xBasisOrder, ...
        'basisType', config.xBasisType, ...
        'basisPattern', config.xBasisPattern ...
        ) ...
        .setElementOperator() ...
        .setDerivativeOrder(1) ...
        .setElementClipper(config.advection) ...
        .setSemiLagrangianElementOperator() ...
        .setUniformGrid(config.nx, config.xBBox) ...
        .setMeshSpace();
else
    xDisc = approx.discretization.FiniteElementDiscretization() ...
        .setDgOrthotopeElement( ...
        'nDims', config.nDims, ...
        'order', config.xBasisOrder, ...
        'basisType', config.xBasisType, ...
        'basisPattern', config.xBasisPattern ...
        ) ...
        .setElementOperator() ...
        .setDerivativeOrder(1) ...
        .setUniformGrid(config.nx, config.xBBox) ...
        .setMeshSpace();
end

fine = physics.State(xDisc);
fine = scheme.initialize(fine);
fine = scheme.run(fine);
fine = scheme.finalize(fine);

%% SAVE
if config.saveAsRef
    REF = scheme.visualizer.dataset.(config.schemeName);
    dataFile = sprintf('%s-REF.mat', config.experimentId);
    dataPath = fullfile(dataDir, dataFile);
    save(dataPath, 'REF');
else
    DG = scheme.visualizer.dataset.(config.schemeName);
    dataFile = sprintf('%s-%s', config.experimentId, config.schemeId);
    dataPath = fullfile(dataDir, dataFile);
    save(dataPath, 'DG');
end

%% HELPER FUNCTIONS
function h = u0(omega)
h = @(x) u0Impl(x, omega);
end

function f = u0Impl(x, omega)
z = omega(:).' * x;
s = all((x <= 0.8) & (x >= 0.3), 1);
f = s .* sin(z) + (1 - s) .* (cos(z) - 1 / 2);
end

function h = uExact(c, omega)
h = @(x, t) uExactImpl(x, t, c, omega);
end

function f = uExactImpl(x, t, c, omega)
z = mod(x - c(:) * t, 1);
s = all((z <= 0.8) & (z >= 0.3), 1);
z = omega(:).' * z;
f = s .* sin(z) + (1 - s) .* (cos(z) - 1 / 2);
end