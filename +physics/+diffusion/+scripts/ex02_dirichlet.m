%==========================================================================
% FileName: ex02_dirichlet
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
config.uExact = uExact(config.lambda, config.omega); % Exact solution
config.nx = repmat(4, 1, config.nDims); % Grid resolution
config.verbose = 1; % Verbose flag
config.TID = upper(config.tOdeIntName);
config.XID = sprintf('%s%d', upper(config.xDiscName), config.xBasisOrder);
config.EID = sprintf('PA%d', config.nDims); % Experiment ID
config.SID = strjoin({config.TID, config.XID}, '-'); % Scheme ID

config.records = {'Integration', 'Visualization', 'Analysis'};
config.timer = core.chrono.Timer(config.records, config.verbose);

config.visualizer = physics.Visualizer(config.nDims);
config.visualizer.setDensity(ones(1, config.nDims));
config.visualizer.setTimeline(config.tFinal, 10);
config.visualizer.setComponents(struct('U', 1));
config.visualizer.setTitleHead(sprintf('%s with %s', config.EID, config.SID));
config.visualizer.addExact('U', config.uExact);
config.visualizer.addDataset('REF', 1);
config.visualizer.addDataset(config.XID, 2);

switch config.nDims
    case 1
        config.visualizer.addPlotter('profile', '1d');
    case 2
        %         config.visualizer.addPlotter('profile', '2d');
        config.visualizer.addPlotter('slice', 'slice1d');
    case 3
        %         config.visualizer.addPlotter('profile', '3d');
        config.visualizer.addPlotter('slice', 'slice1d');
end

config.analyzer = profilers.Analyzer(5);
config.analyzer.setReductions({'L1', 'L2', 'Lx'});
config.analyzer.setDensity(repmat(4, 1, config.nDims));
config.analyzer.setComponents(struct('U', 1));
config.analyzer.addExact('U', config.uExact);


xElement = approx.element.L2OrthotopeElement(config.nDims, ...
    config.xBasisOrder, config.xBasisType, config.xBasisPattern);
scheme = physics.diffusion.LocalScheme(config);

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
h = @(x, t) uExactImpl(x, t, lambda, omega);
end