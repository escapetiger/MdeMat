%==========================================================================
% FileName: ex03_inhomogenous
% Author: Yi CAI
% Date: 13/09/2025
% Description:
%   Solving inhomogeneous diffusion equation with manufactured solution.
%==========================================================================

clc, clear, close all;

%% CONFIGURATION

% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = physics.scheme.ConfigParser.parseFile(common);

% Experiment-specific configuration struct
config.nDims = 1; % Dimension
config.L = 1; % Length
config.xBBox = repmat([0, config.L], 1, config.nDims); % Spatial bounding box
config.tFinal = 0.5; % Final time
config.diffusion = 0.1;
config.ic = u0(); % Initial condition
config.bc = uBc(); % Boundary condition
config.nx = repmat(4, 1, config.nDims); % Grid resolution
config.verbose = 1; % Verbose flag
config.eId = sprintf('IHA%d', config.nDims); % Experiment ID

% Construct and configure scheme
switch config.schemeName
    case 'ldg'
        scheme = physics.diffusion.LdgScheme(config = config);
end

% Configure visualizer
visualizer = scheme.getConfig('visualizer');
visualizer.setDensity(1);
visualizer.setTimeline(scheme.getConfig('tFinal'), 10);
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

% Configure analyzer
analyzer = scheme.getConfig('analyzer');
analyzer.setNLevels(5);
analyzer.setDensity(4^config.nDims);

%% SIMULATION EXECUTION

switch scheme.getConfig('xBasisType')
    case 'modal'
        xElement = approx.element.BH1OrthotopeElement.modal(config.nDims, ...
            scheme.getConfig('xBasisOrder'), ...
            pattern=scheme.getConfig('xBasisPattern'));
    case 'nodal'
        xElement = approx.element.BH1OrthotopeElement.nodal(config.nDims, ...
            scheme.getConfig('xBasisOrder'));
end

xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx/2);
xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
state = physics.state.SpatialState(xDisc);
state = scheme.converge(state);

%% HELPER FUNCTIONS
function h = u0()
h = @(x) u0Impl(x);
end

function f = u0Impl(x)
f = zeros(1, size(x, 2));
end

function h = uBc()
h = @(x, varargin) uBcImpl(x, varargin{:});
end

function f = uBcImpl(x, options)
arguments
    x
    options.t
end
f = any(x == 0, 1) .* 1 + any(x == 1, 1) .* 0;
end