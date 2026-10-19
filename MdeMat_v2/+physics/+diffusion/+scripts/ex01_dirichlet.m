%==========================================================================
% FileName: ex01_dirichlet
% Author: Yi Cai
% Date: 2025-07-13
% Description:
%   Solving diffusion equation with Dirichlet boundary conditions.
%==========================================================================

clc;
clear;
close all;

%% CONFIGURATION
configPath = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = core.io.OptionParser().parse(configPath);
config.nDims = 2; % Dimension
config.L = 1; % Length
config.xBBox = repmat([0, config.L], 1, config.nDims); % Spatial bounding box
config.tFinal = 0.5; % Final time
if config.nDims == 1
    config.diffusion = 0.1;
elseif config.nDims == 2
    config.diffusion = [0.06, 0.01; 0.02, 0.07]; 
else
    config.diffusion = [0.05, -0.01, 0.01; 0.01, 0.06, 0.01;0.01, 0.01, 0.07];
end 
config.omega = repmat(pi / config.L, 1, config.nDims); % Wave length
config.lambda = config.omega * config.diffusion * config.omega.';
config.ic = u0(config.omega); % Initial condition
config.bc = uBc(config.lambda, config.omega); % Boundary condition
config.exact = uExact(config.lambda, config.omega); % Exact solution
% config.exact = [];
config.nx = repmat(4, 1, config.nDims); % Grid resolution
config.experimentId = sprintf('PA%d', config.nDims); % Experiment ID
config.tDiscId = upper(config.tDiscId);
config.xDiscId = sprintf('%s%d', upper(config.xDiscName), config.xBasisOrder);
config.schemeId = strjoin({config.tDiscId, config.xDiscId}, '-');

%% SIMULATION SETUP

scheme = physics.diffusion.LdgScheme(config) ...
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
    .setAnalyzer( ...
    'nLevels', 5, ...
    'reductions', {'L1', 'L2', 'Lx'}, ...
    'density', repmat(4, 1, config.nDims), ...
    'components', struct('U', 1)) ...
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

if isempty(config.exact)
    xDisc = xDisc.setUniformGrid(config.nx/2, config.xBBox).setMeshSpace();
    coarse = physics.State(xDisc);
    coarse = scheme.initialize(coarse);
    coarse = scheme.run(coarse);
    xDisc = xDisc.refineMesh(2);
else
    xDisc = xDisc.setUniformGrid(config.nx, config.xBBox).setMeshSpace();
    coarse = [];
end

for i = 1:scheme.analyzer.profiler.nLevels
    fine = physics.State(xDisc);
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
function h = u0(omega)
h = @(x) u0Impl(x, omega);
end

function f = u0Impl(x, omega)
z = omega(:).' * x;
f = sin(z);
end

function h = uExact(lambda, omega)
h = @(x, t) uExactImpl(x, t, lambda, omega);
end

function f = uExactImpl(x, t, lambda, omega)
    z = omega(:).' * x;
    f = exp(-lambda * t) * sin(z);
end

function h = uBc(lambda, omega)
h = @(i, x, t) uExactImpl(x, t, lambda, omega);
end