function [state, result] = ex03_discontinuous(options)
%   Smoke=true executes this original model at bounded resolution and time,
%   returns actual results, and skips figures/files. Defaults retain the full
%   research parameters and sweeps below.
% Run the retained research study; return its final workspace state and assessment.
arguments
    options.Smoke (1, 1) logical = false
    options.Plot (1, 1) logical = true
    options.Write (1, 1) logical = true
    options.OutputFile (1, 1) string = ""
end
if options.Smoke
    options.Plot = false;
    options.Write = false;
end


%==========================================================================
% FileName: ex03_discontinuous
% Author: Yi CAI
% Date: 13/09/2025
% Description:
%   Solving advection equation with discontinuous initial condition.
%==========================================================================


%% CONFIGURATION

% Data directory
dataDir = fullfile(fileparts(mfilename('fullpath')), 'ex03');

% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = mdemat.tools.environment.readConfig(common);

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

config = mdemat.internal.execution.studyConfig(config, options);

% Construct and configure scheme
switch config.schemeName
    case 'sldg'
        scheme = mdemat.solver.hyperbolic.SldgScheme(config = config);
    case 'uwdg'
        scheme = mdemat.solver.hyperbolic.UwdgScheme(config = config);
end

%% SIMULATION SETUP

refFile = sprintf('%s-REF.mat', scheme.getConfig('eId'));
refPath = fullfile(dataDir, refFile);
if exist(refPath, 'file')
    sample = mdemat.internal.execution.readStudySamples(refPath);
    refData = sample.Values;
    refData.x = sample.Points;
else
    refData = struct();
end

% Configure visualizer
plotOptions = scheme.Observation.Plot;
plotOptions.PointsPerDimension = 1;
plotOptions.Times = linspace(0, scheme.getConfig('tFinal'), 10+1);
plotOptions.DatasetTypes.('REF') = 'fixed';
plotOptions.FixedData.('REF') = mdemat.internal.execution.sampledDataset(refData);
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


scheme.setObservation('plot', plotOptions);
%% SIMULATION EXECUTION

switch config.schemeName
    case 'sldg'
        switch scheme.getConfig('xBasisType')
            case 'modal'
                xElement = mdemat.common.fem.BH1SLOrthotopeElement.modal(config.nDims, ...
                    scheme.getConfig('xBasisOrder'), ...
                    pattern=scheme.getConfig('xBasisPattern'));
            case 'nodal'
                xElement = mdemat.common.fem.BH1SLOrthotopeElement.nodal(config.nDims, ...
                    scheme.getConfig('xBasisOrder'));
        end
    case 'uwdg'
        switch scheme.getConfig('xBasisType')
            case 'modal'
                xElement = mdemat.common.fem.BH1OrthotopeElement.modal(config.nDims, ...
                    scheme.getConfig('xBasisOrder'), ...
                    pattern=scheme.getConfig('xBasisPattern'));
            case 'nodal'
                xElement = mdemat.common.fem.BH1OrthotopeElement.nodal(config.nDims, ...
                    scheme.getConfig('xBasisOrder'));
        end
end

xMesh = mdemat.data.mesh.UniformGrid(config.xBBox, config.nx);
xDisc = mdemat.common.fem.FiniteElementSpace(xMesh, xElement);
state = mdemat.common.fem.SpatialState(xDisc);
state = scheme.run(state);
result = mdemat.internal.execution.studyResult(scheme, state);
plan = mdemat.common.fem.plotSamplingPlan(state.XDisc, plotOptions.PointsPerDimension);
points = plan.Coordinates;
values = struct();
names = fieldnames(plotOptions.Components);
for field = 1:numel(names)
    name = names{field}; ids = plotOptions.Components.(name);
    if isempty(ids) || ~isfield(state.Dofs,name),continue;end
    values.(name) = state.XDisc.eval(plan.Nodes,state.Dofs.(name)(:,ids));
end
result.Samples = struct('type', "mdemat.sampledStudy", 'schemaVersion', 1, ...
    'Points', points, 'Values', values);

%% SAVE
if options.Write || strlength(options.OutputFile) > 0
if strlength(options.OutputFile) == 0 && ~isfolder(dataDir), mkdir(dataDir); end
if config.saveAsRef
    dataFile = sprintf('%s-REF.mat', scheme.getConfig('eId'));
    dataPath = fullfile(dataDir, dataFile);
else
    dataFile = sprintf('%s-%s.mat', scheme.getConfig('eId'), scheme.getConfig('xId'));
    dataPath = fullfile(dataDir, dataFile);
end
if strlength(options.OutputFile) > 0, dataPath = options.OutputFile; end
mdemat.tools.io.writeRecord(dataPath, result.Samples, createDirectory=true);

end

%% HELPER FUNCTIONS

end

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
