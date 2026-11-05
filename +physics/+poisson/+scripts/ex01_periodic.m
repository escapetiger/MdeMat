%==========================================================================
% FileName: ex01_periodic
% Author: Yi Cai
% Date: 2025-08-21
% Description:
%   Solving poisson equation with periodic boundary conditions.
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
config.omega = repmat(pi / config.L, 1, config.nDims); % Wave length
config.force = force(config.omega);
config.bc = []; % Boundary condition
config.uExact = uExact(config.omega); % Exact solution
config.nx = repmat(4, 1, config.nDims); % Grid resolution
config.verbose = 1; % Verbose flag
config.XID = sprintf('%s%d', upper(config.xDiscName), config.xBasisOrder);
config.EID = sprintf('PA%d', config.nDims); % Experiment ID

config.records = {'Integration', 'Visualization', 'Analysis'};
config.timer = core.chrono.Timer(config.records, config.verbose);

config.visualizer = physics.Visualizer(config.nDims);
config.visualizer.setDensity(ones(1, config.nDims));
config.visualizer.setComponents(struct('U', 1));
config.visualizer.setTitleHead(sprintf('%s with %s', config.EID, config.XID));
config.visualizer.addExact('U', config.uExact);
config.visualizer.addDataset('REF', 1);
config.visualizer.addDataset(config.XID, 2);

switch config.nDims
    case 1
        config.visualizer.addPlotter('profile', '1d');
    case 2
%                 config.visualizer.addPlotter('profile', '2d');
        config.visualizer.addPlotter('slice', 'slice1d');
    case 3
        %         config.visualizer.addPlotter('profile', '3d');
        config.visualizer.addPlotter('slice', 'slice1d');
end

config.analyzer = profilers.Analyzer(8);
config.analyzer.setReductions({'L1', 'L2', 'Lx'});
config.analyzer.setDensity(repmat(4, 1, config.nDims));
config.analyzer.setComponents(struct('U', 1));
config.analyzer.addExact('U', config.uExact);

xElement = approx.element.L2OrthotopeElement(config.nDims, ...
    config.xBasisOrder, config.xBasisType, config.xBasisPattern);
scheme = physics.poisson.PoissonScheme(config);

if ~config.analyzer.hasExact
    %< Mesh initialization
    xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx/2);
    xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
    state = physics.SpatialState(xDisc);
    state = scheme.converge(state);
else
    %< Mesh initialization
    xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx);
    xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
    state = physics.SpatialState(xDisc);
    state = scheme.converge(state);
end

%% HELPER FUNCTIONS
function h = uExact(omega)
h = @(x) uExactImpl(x, omega);
end

function f = uExactImpl(x, omega)
z = omega(:).' * x;
f = sin(z);
end

function h = force(omega)
h = @(x) forceImpl(x, omega);
end

function f = forceImpl(x, omega)
z = omega(:).' * x;
f = omega(:).' * omega(:) * sin(z);
end