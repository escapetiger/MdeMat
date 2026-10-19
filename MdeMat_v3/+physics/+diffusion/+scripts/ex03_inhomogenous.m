%==========================================================================
% FileName: ex03_inhomogenous
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
config.diffusion = 0.1;
config.ic = u0(); % Initial condition
config.bc = uBc(); % Boundary condition
config.nx = repmat(4, 1, config.nDims); % Grid resolution
config.verbose = 1; % Verbose flag
config.TID = upper(config.tOdeIntName);
config.XID = sprintf('%s%d', upper(config.xDiscName), config.xBasisOrder);
config.EID = sprintf('IHA%d', config.nDims); % Experiment ID
config.SID = strjoin({config.TID, config.XID}, '-'); % Scheme ID


config.records = {'Integration', 'Visualization', 'Analysis'};
config.timer = core.chrono.Timer(config.records, config.verbose);

config.visualizer = physics.Visualizer(config.nDims);
config.visualizer.setDensity(ones(1, config.nDims));
config.visualizer.setTimeline(config.tFinal, 10);
config.visualizer.setComponents(struct('U', 1));
config.visualizer.setTitleHead(sprintf('%s with %s', config.EID, config.SID));
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

scheme = physics.diffusion.LocalScheme(config);

xElement = approx.element.L2OrthotopeElement(config.nDims, ...
    config.xBasisOrder, config.xBasisType, config.xBasisPattern);
xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx/2);
xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
state = physics.SpatialState(xDisc);
state = scheme.converge(state);

%% HELPER FUNCTIONS
function h = u0()
h = @(x) u0Impl(x);
end

function f = u0Impl(x)
f = zeros(1, size(x, 2));
end

function h = uBc()
h = @(x, t) uBcImpl(x, t);
end

function f = uBcImpl(x, t)
f = any(x == 0, 1) .* 1 + any(x == 1, 1) .* 0;
end