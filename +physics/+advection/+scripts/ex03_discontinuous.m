%==========================================================================
% FileName: ex03_discontinuous
% Author: Yi CAI
% Date: 13/09/2025
% Description:
%   Solving advection equation with discontinuous initial condition.
%==========================================================================

clc, clear, close all;

%% CONFIGURATION

% Data directory
dataDir = fullfile(fileparts(mfilename('fullpath')), 'ex03');

% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = physics.scheme.ConfigParser.parseFile(common);

% Experiment-specific configuration struct
config.nDims = 1; % Dimension
config.xBBox = repmat([0, 1], 1, config.nDims); % Spatial bounding box
config.tFinal = 1.1; % Final time
config.advection = ones(1, config.nDims); % Advection speed
config.omega = repmat(2*pi, 1, config.nDims); % Wave length
config.ic = u0(config.omega); % Initial condition
config.bc = []; % Boundary condition
config.exact = uExact(config.advection, config.omega); % Exact solution
config.nx = repmat(128, 1, config.nDims); % Grid resolution
config.verbose = 1; % Verbose flag
config.eId = sprintf('DA%d', config.nDims); % Experiment ID
config.saveAsRef = 0;

% Construct and configure scheme
switch config.schemeName
    case 'sldg'
        scheme = physics.advection.SldgScheme(config = config);
    case 'uwdg'
        scheme = physics.advection.UwdgScheme(config = config);
end

%% SIMULATION SETUP

refFile = sprintf('%s-REF.mat', scheme.getConfig('eId'));
refPath = fullfile(dataDir, refFile);
if exist(refPath, 'file')
    refData = load(refPath).dataset.Data;
else
    refData = [];
end

% Configure visualizer
visualizer = scheme.getConfig('visualizer');
visualizer.setDensity(1);
visualizer.setTimeline(scheme.getConfig('tFinal'), 10);
visualizer.addDataset('REF', 'fixed', refData);
visualizer.addDataset(scheme.getConfig('xId'), 'numeric');

switch config.nDims
    case 1
        visualizer.addPlotter('profile', '1d');
    case 2
        %         visualizer.addPlotter('profile', '2d');
        visualizer.addPlotter('slice', 'slice1d');
    case 3
        %         visualizer.addPlotter('profile', '3d');
        visualizer.addPlotter('slice', 'slice1d');
end

%% SIMULATION EXECUTION

switch config.schemeName
    case 'sldg'
        switch scheme.getConfig('xBasisType')
            case 'modal'
                xElement = approx.element.BH1SLOrthotopeElement.modal(config.nDims, ...
                    scheme.getConfig('xBasisOrder'), ...
                    pattern=scheme.getConfig('xBasisPattern'));
            case 'nodal'
                xElement = approx.element.BH1SLOrthotopeElement.nodal(config.nDims, ...
                    scheme.getConfig('xBasisOrder'));
        end
    case 'uwdg'
        switch scheme.getConfig('xBasisType')
            case 'modal'
                xElement = approx.element.BH1OrthotopeElement.modal(config.nDims, ...
                    scheme.getConfig('xBasisOrder'), ...
                    pattern=scheme.getConfig('xBasisPattern'));
            case 'nodal'
                xElement = approx.element.BH1OrthotopeElement.nodal(config.nDims, ...
                    scheme.getConfig('xBasisOrder'));
        end
end

xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx);
xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
state = physics.state.SpatialState(xDisc);
state = scheme.run(state);

%% SAVE
if config.saveAsRef
    dataset = visualizer.Database.getDataset(scheme.getConfig('xId'));
    dataFile = sprintf('%s-REF.mat', scheme.getConfig('eId'));
    dataPath = fullfile(dataDir, dataFile);
    save(dataPath, 'dataset');
else
    dataset = visualizer.Database.getDataset(scheme.getConfig('xId'));
    dataFile = sprintf('%s-%s', scheme.getConfig('eId'), scheme.getConfig('xId'));
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
f = s .* sin(z) + (1 - s) .* (cos(z) - 1/2);
end

function h = uExact(c, omega)
h = @(x, varargin) uExactImpl(c, omega, x, varargin{:});
end

function f = uExactImpl(c, omega, x, options)
arguments
    c
    omega
    x
    options.t
end
z = mod(x - c(:) * options.t, 1);
s = all((z <= 0.8) & (z >= 0.3), 1);
z = omega(:).' * z;
f = s .* sin(z) + (1 - s) .* (cos(z) - 1/2);
end