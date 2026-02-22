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
configPath = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = core.io.OptionParser().parse(configPath);
config.nDims = 2; % Dimension
config.experimentId = sprintf('A%d', config.nDims); % Experiment ID
config.schemeId = sprintf('%s%s%d', ...
    upper(config.tDiscId), upper(config.xDiscId), config.xBasisOrder); % Scheme ID
config.L = 1; % Length
config.xBBox = repmat([-config.L, config.L], 1, config.nDims); % Spatial bounding box
config.tFinal = 0.5; % Final time
config.advection = ones(1, config.nDims); % Advection speed
config.omega = repmat(pi / config.L, 1, config.nDims); % Wave length
config.ampl = [0, 1]; % Wave amplitude [sine, cosine]
config.ic = u0(config.omega, config.ampl); % Initial condition
config.bc = []; % Boundary condition
config.exacts.U = uExact(config.advection, config.omega, config.ampl); % Exact solution
% config.exacts = [];
config.nx = repmat(4, 1, config.nDims); % Grid resolution
config.verbose = 1; % Verbose flag

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

switch config.xAdvectionType
    case 1
        cls = 'UwdgScheme';
    case 2
        cls = 'SldgScheme';
end
scheme = physics.advection.(cls)(config, tDisc, visualizer, analyzer);

%% SIMULATION EXECUTION
disc = struct();
disc.x.fe = fem.element.L2OrthotopeFiniteElement( ...
    'nDims', config.nDims, ...
    'order', config.xBasisOrder, ...
    'id', config.xBasisId, ...
    'pattern', config.xBasisPattern ...
    );
disc.x.op = fem.element.L2FiniteElementOperator(disc.x.fe);
if config.xAdvectionType == 2
    disc.x.clp = fem.element.L2FiniteElementClipper(disc.x.fe, config.advection);
    disc.x.slop = fem.element.L2SemiLagrangianFiniteElementOperator(disc.x.fe, disc.x.clp);
end

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
function h = u0(omega, ampl)
h = @(x) u0Impl(x, omega, ampl);

    function f = u0Impl(x, omega, ampl)
    z = omega(:).' * x;
    f = ampl(1) * cos(z) + ampl(2) * sin(z);
    end
end

function h = uExact(c, omega, ampl)
h = @(x, t) uExactImpl(x, t, c, omega, ampl);
    function f = uExactImpl(x, t, c, omega, ampl)
        z = omega(:).' * (x - c(:) * t);
        f = ampl(1) * cos(z) + ampl(2) * sin(z);
    end
end