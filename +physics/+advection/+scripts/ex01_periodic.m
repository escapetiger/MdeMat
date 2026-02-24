%==========================================================================
% FileName: ex01_periodic
% Author: Yi CAI
% Date: 13/09/2025
% Description:
%   Solving advection equation with periodic boundary conditions.
%==========================================================================

clc, clear, close all;

%% CONFIGURATION

% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = physics.scheme.ConfigParser.parseFile(common);

% Experiment-specific configuration struct
config.nDims = 2; % Dimension
config.L = 1; % Length
config.xBBox = repmat([-config.L, config.L], 1, config.nDims); % Spatial bounding box
config.tFinal = 0.5; % Final time
config.advection = ones(1, config.nDims); % Advection speed
config.omega = repmat(pi / config.L, 1, config.nDims); % Wave length
config.ampl = [0, 1]; % Wave amplitude [sine, cosine]
config.ic = u0(config.omega, config.ampl); % Initial condition
config.bc = []; % Boundary condition
config.exact = uExact(config.advection, config.omega, config.ampl); % Exact solution
config.nx = repmat(4, 1, config.nDims); % Grid resolution
config.verbose = 1; % Verbose flag
config.eId = sprintf('PA%d', config.nDims); % Experiment ID

% Construct and configure scheme
switch config.schemeName
    case 'sldg'
        scheme = physics.advection.SldgScheme(config = config);
    case 'uwdg'
        scheme = physics.advection.UwdgScheme(config = config);
end

% Configure visualizer
visualizer = scheme.getConfig('visualizer');
visualizer.setDensity(1);
visualizer.setTimeline(scheme.getConfig('tFinal'), 10);
visualizer.addDataset('REF', 'exact');
visualizer.addDataset(scheme.getConfig('xId'), 'numeric');

switch config.nDims
    case 1
        visualizer.addPlotter('profile', '1d');
    case 2
        visualizer.addPlotter('profile', '2d');
        visualizer.addPlotter('slice', 'slice1d');
    case 3
        % visualizer.addPlotter('profile', '3d');
        visualizer.addPlotter('slice', 'slice1d');
end

% Configure analyzer
analyzer = scheme.getConfig('analyzer');
analyzer.setNLevels(3);
analyzer.setDensity(4^config.nDims);

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

if ~analyzer.HasExact
    %< Mesh initialization
    xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx/2);
    xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
    state = physics.state.SpatialState(xDisc);
    state = scheme.converge(state);
else
    %< Mesh initialization
    xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx);
    xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
    state = physics.state.SpatialState(xDisc);
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
h = @(x, varargin) uExactImpl(c, omega, ampl, x, varargin{:});
end

function f = uExactImpl(c, omega, ampl, x, options)
arguments
    c
    omega
    ampl
    x
    options.t
end
z = omega(:).' * (x - c(:) * options.t);
f = ampl(1) * cos(z) + ampl(2) * sin(z);
end