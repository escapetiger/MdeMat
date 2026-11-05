%==========================================================================
% FileName: ex01_free_streaming
% Author: Yi Cai
% Date: 2025-07-23
% Description:
%   Free streaming with periodic boundary conditions.
%==========================================================================

clc;
clear;
close all;

%% CONFIGURATION
configPath = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = core.io.OptionParser().parse(configPath);
config.nDims = 1; % Dimension
config.L = 1; % Length
config.xBBox = repmat([-config.L, config.L], 1, config.nDims); % Spatial bounding box
config.tFinal = 0.5; % Final time
config.omega = repmat(pi / config.L, 1, config.nDims); % Wave length
config.ic.f = fInit(config.omega); % Initial condition
config.bc.f = []; % Boundary condition
config.exacts.f = fExact(config.omega); % Exact solution for f
% config.exacts = [];
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
    config.nv = max([3, 3], repmat(ceil(sqrt(config.M)) + 1, 1, 2));
else
    config.M = 4;
    config.nv = max([3, 3], repmat(ceil(sqrt(config.M)) + 1, 1, 2));
end
config.N = prod(config.nv);
config.experimentId = sprintf('FSA%dR%d', config.nDims, config.vReductionType);
config.vDiscId = sprintf('P%dS%d', config.M, config.N);
config.tDiscId = upper(config.tDiscId);
config.xDiscId = sprintf('%s%d', upper(config.xDiscId), config.xBasisOrder);
config.schemeId = strjoin({config.tDiscId, config.xDiscId, config.vDiscId}, '-');

%% SIMULATION SETUP

switch lower(config.tDiscId)
    case 'be'
        cls = 'BeIntegrator';
    case 'bdf2'
        cls = 'Bdf2Integrator';
    case 'bdf3'
        cls = 'Bdf3Integrator';
    case 'sdirk2'
        cls = 'Sdirk2Integrator';
    case 'sdirk3'
        cls = 'Sdirk3Integrator';
end
tDisc = approx.odeint.(cls)(config.tFinal);

visualizer = profilers.visual.Visualizer( ...
    'experimentId', config.experimentId, ...
    'schemeId', config.schemeId, ...
    'nDims', config.nDims, ...
    'nTimeNodes', 2, ...
    'final', config.tFinal, ...
    'density', ones(1, config.nDims), ...
    'components', struct('U', 1, 'G', []), ...
    'exacts', config.exacts ...
    );

analyzer = profilers.analysis.Analyzer( ...
    'nLevels', 5, ...
    'reductions', {'L1', 'L2', 'Lx'}, ...
    'density', repmat(4, 1, config.nDims), ...
    'components', struct('U', 1, 'G', []), ...
    'exacts', config.exacts ...
    );

pattern = physics.radiation.MacroMicroPattern( ...
    'reduction', config.vReductionType);
pattern.setDecompositionPattern(config.nDims, config.M);

scheme = physics.radiation.MacroMicroScheme(config, tDisc, visualizer, analyzer, pattern);

%% SIMULATION EXECUTION

disc = struct();
disc.v.fe = fem.element.C0SphereFiniteElement( ...
    'nDims', config.nDims, ...
    'order', config.M, ...
    'reduction', config.vReductionType, ...
    'nPoints', config.nv);
disc.v.space = fem.space.ScaledAffineSpace(disc.v.fe, pattern.scales);
disc.x.fe = fem.element.L2OrthotopeFiniteElement( ...
    'nDims', config.nDims, ...
    'order', config.xBasisOrder, ...
    'id', config.xBasisId, ...
    'pattern', config.xBasisPattern ...
    );
disc.x.op = fem.element.L2FiniteElementOperator(disc.x.fe);

if isempty(config.exacts)
    disp(['[C] nx = ', mat2str(config.nx/2)]);
    disc.x.mesh = approx.mesh.UniformGrid(config.nx/2, config.xBBox);
    disc.x.space = fem.space.MeshSpace(disc.x.fe, disc.x.mesh);
    coarse = physics.radiation.MacroMicroState(disc);
    coarse = scheme.initialize(coarse);
    coarse = scheme.run(coarse);
else
    coarse = [];
end

for i = 1:analyzer.profiler.nLevels
    disp(['[C] nx = ', mat2str(config.nx)]);
    disc.x.mesh = approx.mesh.UniformGrid(config.nx, config.xBBox);
    disc.x.space = fem.space.MeshSpace(disc.x.fe, disc.x.mesh);
    fine = physics.radiation.MacroMicroState(disc);
    fine = scheme.initialize(fine);
    fine = scheme.run(fine);
    if isempty(coarse)
        fine = scheme.finalize(fine);
    else
        fine = scheme.finalize(fine, coarse);
    end
    if ~isempty(coarse), coarse = fine; end
    config.nx = config.nx * 2;
end

%% ANALYSIS
results = analyzer.analyze();
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