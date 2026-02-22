%==========================================================================
% FileName: ex03_isotropic_dirichlet
% Author: Yi Cai
% Date: 2025-07-23
% Description:
%   Diffusive scaling with isotropic boundary conditions.
%==========================================================================

clc;
clear;
close all;

%% CONFIGURATION
dataDir = fullfile(fileparts(mfilename('fullpath')), 'ex03');
if ~exist(dataDir, 'dir'), mkdir(dataDir); end
configPath = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = core.io.OptionParser().parse(configPath);
config.nDims = 1; % Dimension
config.L = 1; % Length
config.xBBox = repmat([0, config.L], 1, config.nDims); % Spatial bounding box
config.tFinal = 1.6; % Final time
config.epsilon = 1; % Scaling parameter
config.alpha = 1; % Time scale
config.beta = -1; % Scattering scale
config.gamma = 1; % Absorption scale
config.nx = repmat(4, 1, config.nDims); % Grid resolution
if config.nDims == 1 && config.vReductionType == 1
    config.E = 1 / 2;
    config.M = 1;
    config.nv = 2;
elseif config.nDims == 1 && config.vReductionType == 2
    config.E = 1 / 2;
    config.M = 2;
    config.nv = max(16, config.M);
elseif config.nDims == 2 && config.vReductionType == 1
    config.E = 1 / (2 * pi);
    config.M = 3;
    config.nv = max(8, ceil((3 * config.M + 2)/4)+1);
elseif config.nDims == 2 && config.vReductionType == 2
    config.E = 1 / (4 * pi);
    config.M = 1;
    config.nv = max([3, 3], repmat(ceil(sqrt(config.M))+1, 1, 2));
else
    config.E = 1 / (4 * pi);
    config.M = 4;
    config.nv = max([4, 4], repmat(ceil(sqrt(config.M))+1, 1, 2));
end
config.N = prod(config.nv);
% config.N = 0;
config.ic = fInit(); % Initial condition
config.bc = fBc(config.L); % Boundary condition
config.exact = [];
config.scattering = scattering();
config.absorption = [];
config.source = [];
if config.nDims < 3
    config.experimentId = sprintf('DSID%dR%d', config.nDims, config.vReductionType);
else
    config.experimentId = sprintf('DSID%d', config.nDims);
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

%% COMPARISON

if config.verbose == 2
    dataset = struct();
    style = struct();
    
    dataset.REF = refDataset;
    style.REF = {'Color', DEFAULT_COLORS{1}, 'Marker', 'none', 'LineStyle', '-', 'LineWidth', 1};
    
    candidates = { ...
        'BE', 1, config.M, config.N; ...
        'SDIRK2', 2, config.M, config.N; ...
        'SDIRK3', 3, config.M, config.N;};
    nCandidates = 0;
    for i = 1:size(candidates, 1)
        candidateId = sprintf('%s-%s%d-P%dS%d', candidates{i, 1}, upper(config.xDiscName), ...
            candidates{i, 2}, candidates{i, 3}, candidates{i, 4});
        candidateName = sprintf('%s%d', upper(config.xDiscName), ...
            candidates{i, 2});
        candidateFile = sprintf('%s.mat', strjoin({config.experimentId, candidateId}, '-'));
        candidatePath = fullfile(dataDir, candidateFile);
        if exist(candidatePath, 'file')
            nCandidates = nCandidates + 1;
            candidateDataset = load(candidatePath).DG;
            candidateDataset.type = 0;
            dataset.(candidateName) = candidateDataset;
            style.(candidateName) = {
                'Color', DEFAULT_COLORS{1+nCandidates},...
                'Marker', DEFAULT_MARKERS{nCandidates}, ...
                'LineStyle', 'none', 'LineWidth', 1};
        end
    end

    scheme.visualizer.reset();
    scheme.visualizer.render([], dataset, style, []);
    ax = gca;
    set(ax,'box','off');
    ax.XAxis.TickValues = 0:0.1:1;
    ax.XAxis.MinorTick = 'on';
    ax.XAxis.MinorTickValues = 0:0.02:1;
    ax.YAxis.TickValues = -0.2:0.2:1.4;
    ax.YAxis.MinorTick = 'on';
    ax.YAxis.MinorTickValues = -0.2:0.04:1.4;
    xlim([0, 1]);
    ylim([-0.2, 1.4]);
    xl = xlabel("$x$", "Interpreter","latex");
    yl = ylabel("$\rho$", "Interpreter","latex");
    xl.FontSize = 16;
    yl.FontSize = 16;
end

%% HELPER FUNCTIONS

function h = fInit()
h = @(x, v) fInitImpl(x, v);
end

function f = fInitImpl(x, v)
f = zeros(1, size(x, 2));
end

function h = fBc(L)
h = @(x, v, t) fBcImpl(x, v, t, L);
end

function f = fBcImpl(x, v, t, L)
fL = 1;
fR = 0;
f = any(x == 0, 1) .* fL + any(x == L, 1) .* fR;
end

function h = scattering()
h = @(x) scatteringImpl(x);
end

function C = scatteringImpl(x)
C = ones(1, size(x, 2));
end