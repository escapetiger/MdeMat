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
%   Solving diffusion equation with periodic boundary conditions.
%==========================================================================


%% CONFIGURATION

% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = mdemat.tools.environment.readConfig(common);

% Experiment-specific configuration struct
config.nDims = 1; % Dimension
config.L = 1; % Length
config.xBBox = repmat([-config.L, config.L], 1, config.nDims); % Spatial bounding box
config.tFinal = 0.5; % Final time
if config.nDims == 1
    config.diffusion = 0.1;
elseif config.nDims == 2
    config.diffusion = [0.05, -0.02; -0.02, 0.05]; 
else
    config.diffusion = [0.05, -0.01, 0.01; 0.01, 0.06, 0.01;0.01, 0.01, 0.07];
end 
config.omega = repmat(pi / config.L, 1, config.nDims); % Wave length
config.lambda = config.omega * config.diffusion * config.omega.';
config.ampl = [0, 1]; % Wave amplitude [sine, cosine]
config.ic = u0(config.omega, config.ampl); % Initial condition
config.bc = []; % Boundary condition
config.exact = uExact(config.lambda, config.omega, config.ampl); % Exact solution
config.nx = repmat(4, 1, config.nDims); % Grid resolution
config.verbose = 1; % Verbose flag
config.eId = sprintf('PA%d', config.nDims); % Experiment ID

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
plotOptions.DatasetTypes.('REF') = 'exact';
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

if isempty(fieldnames(analysisOptions.Exacts))
    %< Mesh initialization
    xMesh = mdemat.data.mesh.UniformGrid(config.xBBox, config.nx/2);
    xDisc = mdemat.common.fem.FiniteElementSpace(xMesh, xElement);
    state = mdemat.common.fem.SpatialState(xDisc);
    [state, result] = scheme.converge(state);
else
    %< Mesh initialization
    xMesh = mdemat.data.mesh.UniformGrid(config.xBBox, config.nx);
    xDisc = mdemat.common.fem.FiniteElementSpace(xMesh, xElement);
    state = mdemat.common.fem.SpatialState(xDisc);
    [state, result] = scheme.converge(state);
end

%% HELPER FUNCTIONS

end

function h = u0(omega, ampl)
h = @(x) u0Impl(x, omega, ampl);
end

function f = u0Impl(x, omega, ampl)
z = omega(:).' * x;
f = ampl(1) * cos(z) + ampl(2) * sin(z);
end

function h = uExact(lambda, omega, ampl)
h = @(x, varargin) uExactImpl(lambda, omega, ampl, x, varargin{:});
end

function f = uExactImpl(lambda, omega, ampl, x, options)
arguments
    lambda
    omega
    ampl
    x
    options.t
end
z = omega(:).' * x;
f = exp(-lambda * options.t) * (ampl(1) * cos(z) + ampl(2) * sin(z));
end
