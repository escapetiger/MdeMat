function [state, result] = ex03_inhomogenous(options)
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
% FileName: ex03_inhomogenous
% Author: Yi CAI
% Date: 13/09/2025
% Description:
%   Solving inhomogeneous diffusion equation with manufactured solution.
%==========================================================================


%% CONFIGURATION

% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = mdemat.tools.environment.readConfig(common);

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

config = mdemat.internal.execution.studyConfig(config, options);

% Construct and configure scheme
switch config.schemeName
    case 'ldg'
        scheme = mdemat.solver.parabolic.LdgScheme(config = config);
end

% Configure visualizer
plotOptions = scheme.Observation.Plot;
plotOptions.PointsPerDimension = 1;
plotOptions.Times = linspace(0, scheme.getConfig('tFinal'), 10+1);
plotOptions.DatasetTypes.(scheme.getConfig('xId')) = 'numeric';

switch config.nDims
    case 1
        if options.Plot, plotOptions.Modes{end+1} = '1d'; end
    case 2
        %         plotOptions.Modes{end+1} = '2d';
        if options.Plot, plotOptions.Modes{end+1} = 'slice1d'; end
    case 3
        %         plotOptions.Modes{end+1} = '3d';
        if options.Plot, plotOptions.Modes{end+1} = 'slice1d'; end
end

% Configure analyzer
analysisOptions = scheme.Observation.Analysis;
if options.Smoke, analysisOptions.NLevels = 1; else, analysisOptions.NLevels = 5; end
analysisOptions.PointsPerDimension = 4;
scheme.setObservation('plot', plotOptions);
scheme.setObservation('analysis', analysisOptions);


%% SIMULATION EXECUTION

switch scheme.getConfig('xBasisType')
    case 'modal'
        xElement = mdemat.common.fem.BH1OrthotopeElement.modal(config.nDims, ...
            scheme.getConfig('xBasisOrder'), ...
            pattern=scheme.getConfig('xBasisPattern'));
    case 'nodal'
        xElement = mdemat.common.fem.BH1OrthotopeElement.nodal(config.nDims, ...
            scheme.getConfig('xBasisOrder'));
end

xMesh = mdemat.data.mesh.UniformGrid(config.xBBox, config.nx/2);
xDisc = mdemat.common.fem.FiniteElementSpace(xMesh, xElement);
state = mdemat.common.fem.SpatialState(xDisc);
[state, result] = scheme.converge(state);

%% HELPER FUNCTIONS

end

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
