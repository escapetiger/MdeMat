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

%% CONFIGURATION SETUP
configPath = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = core.io.OptionParser().parse(configPath);
config.nDims = 1; % Dimension
config.L = 1; % Length
config.xBBox = repmat([-config.L, config.L], 1, config.nDims); % Spatial bounding box
config.tFinal = 0.5; % Final time
config.advection = ones(1, config.nDims); % Advection speed
config.omega = repmat(pi / config.L, 1, config.nDims); % Wave length
config.ampl = [0, 1]; % Wave amplitude [sine, cosine]
config.ic = u0(config.omega, config.ampl); % Initial condition
config.bc = []; % Boundary condition
config.uExact = uExact(config.advection, config.omega, config.ampl); % Exact solution
config.nx = repmat(4, 1, config.nDims); % Grid resolution
config.records = {'Integration', 'Visualization', 'Analysis'};
config.verbose = 1; % Verbose flag
config.TID = upper(config.tOdeIntName);
config.XID = sprintf('%s%d', upper(config.xDiscName), config.xBasisOrder);
config.EID = sprintf('PA%d', config.nDims); % Experiment ID
config.SID = strjoin({config.TID, config.XID}, '-'); % Scheme ID

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

%% SIMULATION EXECUTION

if contains(config.XID, 'SL')
    xElement = approx.element.L2SLOrthotopeElement(config.nDims, ...
        config.xBasisOrder, config.xBasisType, config.xBasisPattern);
    scheme = physics.advection.SemiLagrangianScheme(config);
else
    xElement = approx.element.L2OrthotopeElement(config.nDims, ...
        config.xBasisOrder, config.xBasisType, config.xBasisPattern);
    scheme = physics.advection.UpwindScheme(config);
end

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
function h = u0(omega, ampl)
h = @(x) u0Impl(x, omega, ampl);
end

function f = u0Impl(x, omega, ampl)
z = omega(:).' * x;
f = ampl(1) * cos(z) + ampl(2) * sin(z);
end

function h = uExact(c, omega, ampl)
h = @(x, t) uExactImpl(x, t, c, omega, ampl);
end

function f = uExactImpl(x, t, c, omega, ampl)
z = omega(:).' * (x - c(:) * t);
f = ampl(1) * cos(z) + ampl(2) * sin(z);
end