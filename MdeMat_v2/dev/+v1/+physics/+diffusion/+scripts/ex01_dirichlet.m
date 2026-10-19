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
config.nDims = 1; % Dimension
config.experimentId = sprintf('A%d', config.nDims); % Experiment ID
config.schemeId = sprintf('%s%s%d', ...
    upper(config.tDiscId), upper(config.xDiscId), config.xBasisOrder); % Scheme ID
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
config.bc = uExact(config.lambda, config.omega); % Boundary condition
config.exacts.U = uExact(config.lambda, config.omega); % Exact solution
% config.exacts = [];
config.nx = repmat(4, 1, config.nDims); % Grid resolution

%% SIMULATION SETUP
visualizer = profilers.visual.Visualizer( ...
    'nDims', config.nDims, ...
    'nTimeNodes', 10, ...
    'final', config.tFinal, ...
    'density', ones(1, config.nDims), ...
    'components', struct('U', 1), ...
    'exacts', config.exacts ...
    );

analyzer = profilers.analysis.Analyzer( ...
    'nLevels', 5, ...
    'reductions', {'L1', 'L2', 'Lx'}, ...
    'density', repmat(4, 1, config.nDims), ...
    'components', struct('U', 1), ...
    'exacts', config.exacts ...
    );

switch config.tDiscId
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

scheme = physics.diffusion.LdgScheme(config, tDisc, visualizer, analyzer);

%% SIMULATION EXECUTION

disc.x.fe = fem.element.L2OrthotopeFiniteElement( ...
    'nDims', config.nDims, ...
    'order', config.xBasisOrder, ...
    'id', config.xBasisId, ...
    'pattern', config.xBasisPattern ...
    );
disc.x.op = fem.element.L2FiniteElementOperator(disc.x.fe);

if isempty(config.exacts)
    disc.x.mesh = approx.mesh.UniformGrid(config.nx/2, config.xBBox);
    disc.x.space = fem.space.MeshSpace(disc.x.fe, disc.x.mesh);
    coarse = physics.State(disc);
    coarse = scheme.initialize(coarse);
    coarse = scheme.run(coarse);
else
    coarse = [];
end

for i = 1:analyzer.profiler.nLevels
    disc.x.mesh = approx.mesh.UniformGrid(config.nx, config.xBBox);
    disc.x.space = fem.space.MeshSpace(disc.x.fe, disc.x.mesh);
    fine = physics.State(disc);
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
function h = u0(omega)
h = @(x) u0Impl(x, omega);

    function f = u0Impl(x, omega)
    z = omega(:).' * x;
    f = sin(z);
    end
end

function h = uExact(lambda, omega)
h = @(x, t) uExactImpl(x, t, lambda, omega);
    function f = uExactImpl(x, t, lambda, omega)
        z = omega(:).' * x;
        f = exp(-lambda * t) * sin(z);
    end
end