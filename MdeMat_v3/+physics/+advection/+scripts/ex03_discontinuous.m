%==========================================================================
% FileName: ex03_discontinuous
% Author: Yi Cai
% Date: 2025-07-13
% Description:
%   Solving advection equation with discontinuous initial condition.
%==========================================================================

clc;
clear;
close all;

%% CONFIGURATION
dataDir = fullfile(fileparts(mfilename('fullpath')), 'ex03');
configPath = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = core.io.OptionParser().parse(configPath);
config.nDims = 1; % Dimension
config.xBBox = repmat([0, 1], 1, config.nDims); % Spatial bounding box
config.tFinal = 1.1; % Final time
config.advection = ones(1, config.nDims); % Advection speed
config.omega = repmat(2*pi, 1, config.nDims); % Wave length
config.ic = u0(config.omega); % Initial condition
config.bc = []; % Boundary condition
config.exact = uExact(config.advection, config.omega); % Exact solution
config.nx = repmat(128, 1, config.nDims); % Grid resolution
config.records = {'Integration', 'Visualization'};
config.verbose = 1; % Verbose flag
config.TID = upper(config.tOdeIntName);
config.XID = sprintf('%s%d', upper(config.xDiscName), config.xBasisOrder);
config.EID = sprintf('DA%d', config.nDims); % Experiment ID
config.SID = strjoin({config.TID, config.XID}, '-'); % Scheme ID
config.saveAsRef = 0;

%% SIMULATION SETUP

refFile = sprintf('%s-REF.mat', config.EID);
refPath = fullfile(dataDir, refFile);
if exist(refPath, 'file')
    refData = load(refPath).dataset.data;
else
    refData = [];
end

config.timer = core.chrono.Timer(config.records, config.verbose);

config.visualizer = physics.Visualizer(config.nDims);
config.visualizer.setDensity(ones(1, config.nDims));
config.visualizer.setTimeline(config.tFinal, 10);
config.visualizer.setComponents(struct('U', 1));
config.visualizer.setTitleHead(sprintf('%s with %s', config.EID, config.SID));
config.visualizer.addDataset('REF', 0, refData);
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

xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx);
xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
state = physics.SpatialState(xDisc);
state = scheme.run(state);
if config.verbose > 0, config.timer.report(); end

%% SAVE
if config.saveAsRef
    dataset = config.visualizer.database.datasets.(config.XID);
    dataFile = sprintf('%s-REF.mat', config.EID);
    dataPath = fullfile(dataDir, dataFile);
    save(dataPath, 'dataset');
else
    dataset = config.visualizer.database.datasets.(config.XID);
    dataFile = sprintf('%s-%s', config.EID, config.XID);
    dataPath = fullfile(dataDir, dataFile);
    save(dataPath, 'dataset');
end

%% HELPER FUNCTIONS
function h = u0(omega)
h = @(x) u0Impl(x, omega);
end

function f = u0Impl(x, omega)
z = omega(:).' * x;
s = all((x <= 0.8) & (x >= 0.3), 1);
f = s .* sin(z) + (1 - s) .* (cos(z) - 1 / 2);
end

function h = uExact(c, omega)
h = @(x, t) uExactImpl(x, t, c, omega);
end

function f = uExactImpl(x, t, c, omega)
z = mod(x-c(:)*t, 1);
s = all((z <= 0.8) & (z >= 0.3), 1);
z = omega(:).' * z;
f = s .* sin(z) + (1 - s) .* (cos(z) - 1 / 2);
end