%==========================================================================
% FileName: ex01_free_streaming
% Author: Yi CAI
% Date: 13/09/2025
% Description:
%   Free streaming with periodic boundary conditions.
%==========================================================================

clc;
clear;
close all;

%% CONFIGURATION
configPath = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = physics.scheme.ConfigParser.parseFile(configPath);
config.nDims = 2; % Dimension
config.L = 1; % Length
config.xBBox = repmat([-config.L, config.L], 1, config.nDims); % Spatial bounding box
config.tFinal = 0.5; % Final time
config.omega = repmat(pi / config.L, 1, config.nDims); % Wave length
config.ic = fInit(config.omega); % Initial condition
config.bc = []; % Boundary condition
config.exact = fExact(config.omega); % Exact solution for f
% config.exact = [];
config.scattering = [];
config.absorption = [];
config.source = [];
config.verbose = 1; % Verbose flag
config.nx = repmat(4, 1, config.nDims); % Grid resolution
if config.nDims == 1 && config.vReductionType == 1
    config.M = 1;
    config.nv = 2;
elseif config.nDims == 1 && config.vReductionType == 2
    config.M = 2;
    config.nv = max(8, config.M);
elseif config.nDims == 2 && config.vReductionType == 1
    config.M = 3;
    config.nv = max(8, ceil((3*config.M + 2) / 4) + 1);
elseif config.nDims == 2 && config.vReductionType == 2
    config.M = 3;
    config.nv = max([4, 4], repmat(ceil(sqrt(config.M)) + 1, 1, 2));
else
    config.M = 4;
    config.nv = max([4, 4], repmat(ceil(sqrt(config.M)) + 1, 1, 2));
end
config.N = prod(config.nv);
% config.N = 0;
if config.nDims < 3
    config.experimentId = sprintf('FSA%dR%d', config.nDims, config.vReductionType);
else
    config.experimentId = sprintf('FSA%d', config.nDims);
end
config.vDiscId = sprintf('P%dS%d', config.M, config.N);
config.tDiscId = upper(config.tDiscId);
config.xDiscId = sprintf('%s%d', upper(config.xDiscName), config.xBasisOrder);
config.schemeId = strjoin({config.tDiscId, config.xDiscId, config.vDiscId}, '-');

%% SIMULATION SETUP

scheme = physics.radiation.MacroMicroScheme(config) ...
    .setTimeDiscretization(config.tDiscId, config.tFinal) ...
    .setTimer(config.verbose) ...
    .setVisualizer( ...
    'experimentId', config.experimentId, ...
    'schemeId', config.schemeId, ...
    'nDims', config.nDims, ...
    'nTimeNodes', 2, ...
    'final', config.tFinal, ...
    'density', ones(1, config.nDims), ...
    'components', struct('U', 1, 'G', [])) ...
    .setAnalyzer( ...
    'nLevels', 5, ...
    'reductions', {'L1', 'L2', 'Lx'}, ...
    'density', repmat(4, 1, config.nDims), ...
    'components', struct('U', 1, 'G', [])) ...
    .setPattern( ...
    'reduction', config.vReductionType) ...
    .addDataset('DG', 2) ...
    .addStyle('DG', {'Color', 'b', 'Marker', 'o', 'LineStyle', 'none'}) ...
    .addDataset('REF', ~isempty(config.exact)) ...
    .addStyle('REF', {'Color', 'r', 'Marker', 'none', 'LineStyle', '-'});

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
    .setScaledAffineSpace(scheme.pattern.scales);

if isempty(config.exact)
    xDisc = xDisc.setUniformGrid(config.nx/2, config.xBBox).setMeshSpace();
    coarse = physics.radiation.MacroMicroState(xDisc, vDisc, config.M, config.N);
    coarse = scheme.initialize(coarse);
    coarse = scheme.run(coarse);
    xDisc = xDisc.refineMesh(2);
else
    xDisc = xDisc.setUniformGrid(config.nx, config.xBBox).setMeshSpace();
    coarse = [];
end

for i = 1:scheme.analyzer.profiler.nLevels
    fine = physics.radiation.MacroMicroState(xDisc, vDisc, config.M, config.N);
    fine = scheme.initialize(fine);
    fine = scheme.run(fine);
    if isempty(coarse)
        fine = scheme.finalize(fine);
    else
        fine = scheme.finalize(fine, coarse);
    end
    if ~isempty(coarse), coarse = fine; end
    xDisc = xDisc.refineMesh(2);
end

%% ANALYSIS
results = scheme.analyzer.analyze();
varNames = fieldnames(results);
for iVar = 1:length(varNames)
    varName = varNames{iVar};
    disp(results.(varName));
end

%% HELPER FUNCTIONS

function h = fInit(omega)
h = @(x, v) fInitImpl(x, v, omega);
end

function f = fInitImpl(x, v, omega)
z = omega(:).' * x;
f = sin(z);
end

function h = fExact(omega)
h = @(x, v, t) fExactImpl(x, v, t, omega);
end

function f = fExactImpl(x, v, t, omega)
z = omega(:).' * (x - v .* t);
f = sin(z);
end