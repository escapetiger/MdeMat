%==========================================================================
% FileName: ex04_gaussian_distribution
% Author: Yi Cai
% Date: 2025-08-01
% Description:
%   Diffusive scaling with Gaussian initial condition.
%==========================================================================

clc;
clear;
close all;

%% CONFIGURATION
dataDir = fullfile(fileparts(mfilename('fullpath')), 'ex03');
if ~exist(dataDir, 'dir'), mkdir(dataDir); end
configPath = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = core.io.OptionParser().parse(configPath);
config.nDims = 2; % Dimension
config.L = 1; % Length
config.xBBox = repmat([-config.L, config.L], 1, config.nDims); % Spatial bounding box
config.tFinal = 0.5; % Final time
config.epsilon = 1; % Scaling parameter
config.alpha = 1; % Time scale
config.beta = 0; % Scattering scale
config.gamma = 0; % Absorption scale
config.T = 3.2e-4; % Variance of Gaussian distribution
config.nx = repmat(64, 1, config.nDims); % Grid resolution
if config.nDims == 1 && config.vReductionType == 1
    config.E = 1 / 2;
    config.M = 1;
    config.nv = 2;
elseif config.nDims == 1 && config.vReductionType == 2
    config.E = 1 / 2;
    config.M = 2;
    config.nv = max(8, config.M);
elseif config.nDims == 2 && config.vReductionType == 1
    config.E = 1 / (2 * pi);
    config.M = 3;
    config.nv = max(16, ceil((3 * config.M + 2)/4)+1);
elseif config.nDims == 2 && config.vReductionType == 2
    config.E = 1 / (4 * pi);
    config.M = 3;
    config.nv = max([4, 4], repmat(ceil(sqrt(config.M))+1, 1, 2));
else
    config.E = 1 / (4 * pi);
    config.M = 4;
    config.nv = max([4, 4], repmat(ceil(sqrt(config.M))+1, 1, 2));
end
config.N = prod(config.nv);
% config.N = 0;
config.ic = fInit(config.T); % Initial condition
config.bc = fBc(config.L); % Boundary condition
config.exact = [];
config.scattering = scattering();
config.absorption = [];
config.source = [];
if config.nDims < 3
    config.experimentId = sprintf('DSG%dR%d', config.nDims, config.vReductionType);
else
    config.experimentId = sprintf('DSG%d', config.nDims);
end
config.vDiscId = sprintf('P%dS%d', config.M, config.N);
config.tDiscId = upper(config.tDiscId);
config.xDiscId = sprintf('%s%d', upper(config.xDiscName), config.xBasisOrder);
config.schemeId = strjoin({config.tDiscId, config.xDiscId, config.vDiscId}, '-');
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
schemeName = strrep(config.schemeId, '-', '_');

scheme = physics.radiation.MacroMicroScheme(config) ...
    .setTimeDiscretization(config.tDiscId, config.tFinal) ...
    .setTimer(config.verbose) ...
    .setVisualizer( ...
    'experimentId', config.experimentId, ...
    'schemeId', config.schemeId, ...
    'nDims', config.nDims, ...
    'nTimeNodes', 10, ...
    'final', config.tFinal, ...
    'density', ones(1, config.nDims), ...
    'components', struct('U', 1, 'G', [])) ...
    .setPattern( ...
    'epsilon', config.epsilon, ...
    'timeScale', config.alpha, ...
    'scatteringScale', config.beta, ...
    'absorptionScale', config.gamma, ...
    'reduction', config.vReductionType) ...
    .addDataset('REF', 0, refDataset) ...
    .addStyle('REF', ...
    {'Color', DEFAULT_COLORS{1}, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', 1}) ...
    .addDataset(config.schemeName, 2) ...
    .addStyle(config.schemeName, ...
    {'Color', DEFAULT_COLORS{2}, 'Marker', DEFAULT_MARKERS{1}, 'LineStyle', 'none', 'LineWidth', 1});

%% SIMULATION EXECUTION

xDisc = approx.discretization.FiniteElementDiscretization() ...
    .setDgOrthotopeElement( ...
    'nDims', config.nDims, ...
    'order', config.xBasisOrder, ...
    'basisType', config.xBasisType, ...
    'basisPattern', config.xBasisPattern ...
    ) ...
    .setElementOperator() ...
    .setDerivativeOrder(1);

vDisc = approx.discretization.AffineDiscretization() ...
    .setC0SphereElement( ...
    'nDims', config.nDims, ...
    'order', config.M, ...
    'reduction', config.vReductionType, ...
    'nPoints', config.nv) ...
    .setDerivativeOrder(0) ...
    .setScaledAffineSpace([]);

xDisc = xDisc.setUniformGrid(config.nx, config.xBBox).setMeshSpace();
state = physics.radiation.MacroMicroState(xDisc, vDisc, config.M, config.N);
state = scheme.initialize(state);
state = scheme.run(state);
state = scheme.finalize(state);

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

function h = fInit(T)
h = @(x, v) fInitImpl(x, v, T);
end

function f = fInitImpl(x, v, T)
f = exp(-sum(x.^2, 1) ./ (2*T)) ./ sqrt(2*pi*T)^(size(x, 1));
end

function h = fBc(L)
h = @(x, v, t) fBcImpl(x, v, t, L);
end

function f = fBcImpl(x, v, t, L)
f = zeros(1, size(x, 2));
end

function h = scattering()
h = @(x) scatteringImpl(x);
end

function C = scatteringImpl(x)
% C = ones(1, size(x, 2));
C = zeros(1, size(x, 2));
end