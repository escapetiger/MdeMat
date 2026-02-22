%==========================================================================
% FileName: ex01_periodic
% Author: Yi CAI
% Date: 13/09/2025
% Description:
%   Solving Poisson equation with periodic boundary conditions.
%==========================================================================

clc, clear, close all;

%% CONFIGURATION

% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = physics.scheme.ConfigParser.parseFile(common);

% Experiment-specific configuration struct
config.nDims = 1; % Dimension
config.L = 1; % Length
config.xBBox = repmat([-config.L, config.L], 1, config.nDims); % Spatial bounding box
config.omega = repmat(pi / config.L, 1, config.nDims); % Wave length
config.force = force(config.omega); % Force function
config.bc = []; % Boundary condition (periodic)
config.exact.U = uExact(config.omega); % Exact solution for u
config.exact.Q = qExact(config.omega); % Exact solution for q
config.nx = repmat(8, 1, config.nDims); % Grid resolution
config.verbose = 1; % Verbose flag
config.eId = sprintf('PA%d', config.nDims); % Experiment ID

% Construct and configure scheme
switch config.schemeName
    case 'pdg'
        scheme = physics.poisson.PdgScheme(config = config);
    case 'ldg'
        scheme = physics.poisson.LdgScheme(config = config);
end

% Configure visualizer
visualizer = scheme.getConfig('visualizer');
visualizer.setDensity(1);
visualizer.addDataset('REF', 'exact');
visualizer.addDataset(scheme.getConfig('xId'), 'numeric');
switch config.nDims
    case 1
        visualizer.addPlotter('profile', '1d');
    case 2
        % visualizer.addPlotter('profile', '2d');
        visualizer.addPlotter('slice', 'slice1d');
    case 3
        % visualizer.addPlotter('profile', '3d');
        visualizer.addPlotter('slice', 'slice1d');
end

% Configure analyzer
analyzer = scheme.getConfig('analyzer');
analyzer.setNLevels(5);
analyzer.setDensity(4^config.nDims);

%% SIMULATION

switch scheme.getConfig('xBasisType')
    case 'modal'
        xElement = approx.element.BH1OrthotopeElement.modal(config.nDims, ...
            scheme.getConfig('xBasisOrder'), ...
            pattern=scheme.getConfig('xBasisPattern'));
    case 'nodal'
        xElement = approx.element.BH1OrthotopeElement.nodal(config.nDims, ...
            scheme.getConfig('xBasisOrder'));
end

if ~analyzer.HasExact
    xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx/2);
    xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
    state = physics.state.SpatialState(xDisc);
    state = scheme.converge(state);
else
    xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx);
    xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
    state = physics.state.SpatialState(xDisc);
    state = scheme.converge(state);
end

%% HELPER FUNCTIONS
function h = uExact(omega)
h = @(x) uExactImpl(omega, x);
end

function f = uExactImpl(omega, x)
z = omega(:).' * x;
f = sin(z);
end

function h = qExact(omega)
h = @(x) qExactImpl(omega, x);
end

function f = qExactImpl(omega, x)
z = omega(:).' * x;
f = omega(:) .* cos(z);
f = f.';
end

function h = force(omega)
h = @(x) forceImpl(omega, x);
end

function f = forceImpl(omega, x)
z = omega(:).' * x;
f = omega(:).' * omega(:) * sin(z);
end