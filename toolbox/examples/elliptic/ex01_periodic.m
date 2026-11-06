function [state, result] = ex01_periodic(options)
%   Smoke=true executes this original model at bounded resolution and time,
%   returns actual results, and skips figures/files. Defaults retain the full
%   research parameters and sweeps below.
% Run the retained research study; return its final workspace state and assessment.
arguments
    options.Smoke (1, 1) logical = false
    options.Plot (1, 1) logical = true
    options.Write (1, 1) logical = true
end
if options.Smoke
    options.Plot = false;
    options.Write = false;
end


%==========================================================================
% FileName: ex01_periodic
% Author: Yi CAI
% Date: 13/09/2025
% Description:
%   Solving Poisson equation with periodic boundary conditions.
%==========================================================================


%% CONFIGURATION

% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = mdemat.tools.environment.readConfig(common);

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

config = mdemat.internal.execution.studyConfig(config, options);

% Construct and configure scheme
switch config.schemeName
    case 'pdg'
        scheme = mdemat.solver.elliptic.PdgScheme(config = config);
    case 'ldg'
        scheme = mdemat.solver.elliptic.LdgScheme(config = config);
end

% Configure visualizer
plotOptions = scheme.Observation.Plot;
plotOptions.PointsPerDimension = 1;
plotOptions.DatasetTypes.('REF') = 'exact';
plotOptions.DatasetTypes.(scheme.getConfig('xId')) = 'numeric';
switch config.nDims
    case 1
        if options.Plot, plotOptions.Modes{end+1} = '1d'; end
    case 2
        % plotOptions.Modes{end+1} = '2d';
        if options.Plot, plotOptions.Modes{end+1} = 'slice1d'; end
    case 3
        % plotOptions.Modes{end+1} = '3d';
        if options.Plot, plotOptions.Modes{end+1} = 'slice1d'; end
end

% Configure analyzer
analysisOptions = scheme.Observation.Analysis;
if options.Smoke, analysisOptions.NLevels = 1; else, analysisOptions.NLevels = 5; end
analysisOptions.PointsPerDimension = 4;
scheme.setObservation('plot', plotOptions);
scheme.setObservation('analysis', analysisOptions);


%% SIMULATION

switch scheme.getConfig('xBasisType')
    case 'modal'
        xElement = mdemat.common.fem.BH1OrthotopeElement.modal(config.nDims, ...
            scheme.getConfig('xBasisOrder'), ...
            pattern=scheme.getConfig('xBasisPattern'));
    case 'nodal'
        xElement = mdemat.common.fem.BH1OrthotopeElement.nodal(config.nDims, ...
            scheme.getConfig('xBasisOrder'));
end

if isempty(fieldnames(analysisOptions.Exacts))
    xMesh = mdemat.data.mesh.UniformGrid(config.xBBox, config.nx/2);
    xDisc = mdemat.common.fem.FiniteElementSpace(xMesh, xElement);
    state = mdemat.common.fem.SpatialState(xDisc);
    [state, result] = scheme.converge(state);
else
    xMesh = mdemat.data.mesh.UniformGrid(config.xBBox, config.nx);
    xDisc = mdemat.common.fem.FiniteElementSpace(xMesh, xElement);
    state = mdemat.common.fem.SpatialState(xDisc);
    [state, result] = scheme.converge(state);
end

%% HELPER FUNCTIONS

end

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
end

function h = force(omega)
h = @(x) forceImpl(omega, x);
end

function f = forceImpl(omega, x)
z = omega(:).' * x;
f = omega(:).' * omega(:) * sin(z);
end
