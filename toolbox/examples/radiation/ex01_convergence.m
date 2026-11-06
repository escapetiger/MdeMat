function [state, result] = ex01_convergence(options)
%   Smoke=true executes this original model at bounded resolution and time,
%   returns actual results, and skips figures/files. Defaults retain the full
%   research parameters and sweeps below.
% Run the retained research study; return its final workspace state and assessment.
arguments
    options.Smoke (1, 1) logical = false
end
options.Plot = ~options.Smoke;
options.Write = ~options.Smoke;
state = []; result = [];

%==========================================================================
% FileName: ex01_convergence
% Author: Yi CAI
% Description:
%   Convergence test for radiation transport.
%==========================================================================


%% EXECUTION

% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = mdemat.tools.environment.readConfig(common);

% Problem setup
config.id = 'ex01_2dt_ds';


% Parse settings from problem ID
tok = regexp(config.id, '_(\d+)d_|_(\d+)d([ts])_', 'tokens');
assert(~isempty(tok), 'String format not recognized.');
config.nDims = str2double(tok{1}{1});
if length(tok{1}) >= 2
    switch tok{1}{2}
        case 't'
            config.vDimReduction = 'topology';
        case 's'
            config.vDimReduction = 'symmetry';
    end
else
    config.vDimReduction = '';
end

% Equilibrium
if config.nDims == 1
    config.E = 1 / 2;
elseif config.nDims == 2 && strcmpi(config.vDimReduction, 'topology')
    config.E = 1 / (2*pi);
else
    config.E = 1 / (4*pi);
end

% Spatial domain [-L, L]^d
config.L = 1; % Length
config.xBBox = repmat([-config.L, config.L], 1, config.nDims); % Spatial bounding box
config.nx = repmat(4, 1, config.nDims); % Grid resolution

% Velocity discretization parameters based on dimensions
switch config.nDims
    case 1
        if strcmp(config.vDimReduction, 'topology')
            config.nu = 1; % order: 0 or 1
            config.nv = 2; % number of velocities
        else
            config.nu = 4; % order: 0 for none; nu - 1 for degree
            config.nv = 8; % number of velocities
        end
    case 2
        if strcmp(config.vDimReduction, 'topology')
            config.nu = 1; % order: 0 for none; nu - 1 for degree
            config.nv = 16; % number of velocities
        else
            config.nu = 2; % order: 0 for none; nu - 1 for degree
            config.nv = max([4, 4], [config.nu + 1, 2*(config.nu-1)+2]); % number of velocities
        end
    case 3
        config.nu = 4; % order: 0 for none; nu - 1 for degree
        config.nv = max(1, [config.nu + 1, 2*(config.nu-1)+2]); % number of velocities
end

% Checkpoint information
config.ckptDir = fullfile(fileparts(mfilename('fullpath')), config.id);
if options.Write && ~isfolder(config.ckptDir), mkdir(config.ckptDir); end

config.verbose = 1;
config.tFinal = 0.5;
config.cfl = 0.5;
config.dt = [];
config.useFilter = false;
config.usePositivityLimiter = false;

if options.Smoke
    config.epsilon = 1;
    config.tau0 = 1;
    [config, scheme, state, result] = run(config, options); %#ok<ASGLU>
    return;
end

% Parameter sweep
epsilon = 10.^0;
order = [2,3];
for k = 1:length(epsilon)
    config.epsilon = epsilon(k);
    config.ckptPostfix = sprintf('epsilon%.0e', epsilon(k));

    for j = 1:length(order)
        switch order(j)
            case 1
                config.tOdeIntName = 'be';
                config.xBasisOrder = 1;
            case 2
                config.tOdeIntName = 'sdirk2';
                config.xBasisOrder = 2;
            case 3
                config.tOdeIntName = 'sdirk3';
                config.xBasisOrder = 3;
        end
        config.ckptPrefix = [upper(config.schemeName), num2str(config.xBasisOrder)];
        switch config.id
            case 'ex01_1dt_ft'
                config.ckptPrefix = strjoin({'FT1DT', config.ckptPrefix}, '_');
            case 'ex01_1ds_ft'
                config.ckptPrefix = strjoin({'FT1DS', config.ckptPrefix}, '_');
            case 'ex01_2dt_ft'
                config.ckptPrefix = strjoin({'FT2DT', config.ckptPrefix}, '_');
            case 'ex01_2ds_ft'
                config.ckptPrefix = strjoin({'FT2DS', config.ckptPrefix}, '_');
            case 'ex01_3d_ft'
                config.ckptPrefix = strjoin({'FT3D', config.ckptPrefix}, '_');
            case 'ex01_1dt_ds'
                config.ckptPrefix = strjoin({'DS1DT', config.ckptPrefix}, '_');
            case 'ex01_1ds_ds'
                config.ckptPrefix = strjoin({'DS1DS', config.ckptPrefix}, '_');
            case 'ex01_2dt_ds'
                config.ckptPrefix = strjoin({'DS2DT', config.ckptPrefix}, '_');
        end
        config.ckptTimeStamps = [];
        fileName = sprintf('%s_state_%s.mat', config.ckptPrefix, config.ckptPostfix);
        fileName = fullfile(config.ckptDir, fileName);
    
        [config, ~, state, result] = run(config, options);
        state.save(fileName, omitCallbacks=true, createDirectory=true);
    
        fieldNames = fieldnames(result);
        for iName = 1:length(fieldNames)
            fieldName = fieldNames{iName};
            fileName = sprintf('%s_%s_%s.csv', config.ckptPrefix, fieldName, config.ckptPostfix);
            fileName = fullfile(config.ckptDir, fileName);
            writetable(result.(fieldName), fileName);
        end
    end
end

%% SIMULATION


if isempty(result) && ~isempty(state)
    result = struct('schemaVersion', 1, 'completed', true, 'diagnostics', state.History);
end

end

function [config, scheme, state, result] = run(config, options)

% Physical parameters
config.omega = repmat(pi / config.L, 1, config.nDims); % Wave number
if contains(config.id, 'ft')
    config.decomposition = '';
    config.ic = fInitFs(config);
    config.bc = [];
    config.exact = fExactFs(config);
    config.scattering = [];
    config.absorption = [];
    config.source = [];
else
    config.lambda = 0.1;
    config.decomposition = '';
    config.timeScale = 1;
    config.scatteringScale = -1;
    config.absorptionScale = 1;
    config.ic = fInitDs(config);
    config.bc = [];
    % config.exact = [];
    % config.scattering = scatteringDs(config);
    % config.absorption = [];
    % config.source = [];
    config.exact = fExactDs(config);
    config.scattering = scatteringDs(config);
    config.absorption = absorptionDs(config);
    config.source = sourceDs(config);
end

config = mdemat.internal.execution.studyConfig(config, options);

% Construct and configure scheme
switch config.schemeName
    case 'mmdg'
        scheme = mdemat.application.radiation.MmDgScheme(config=config);
end

% Configure visualizer
plotOptions = scheme.Observation.Plot;
plotOptions.PointsPerDimension = 1;
plotOptions.Times = linspace(0, scheme.getConfig('tFinal'), 1+1);
plotOptions.Components = struct('U', [], 'G', [], 'F', 1);
plotOptions.DatasetTypes.('REF') = 'exact';
plotOptions.DatasetTypes.(scheme.getConfig('sId')) = 'numeric';
switch config.nDims
    case 1
        if options.Plot, plotOptions.Modes{end+1} = '1d'; end
    case 2
        % plotOptions.Modes{end+1} = '2d';
        if options.Plot, plotOptions.Modes{end+1} = 'slice1d'; end
    case 3
        if options.Plot, plotOptions.Modes{end+1} = 'slice1d'; end
end

% Configure analyzer
analysisOptions = scheme.Observation.Analysis;
if options.Smoke, analysisOptions.NLevels = 1; else, analysisOptions.NLevels = 5; end
analysisOptions.PointsPerDimension = 4;
M = scheme.getConfig('M');
N = scheme.getConfig('N');
analysisOptions.Components = struct('U', 1:M, 'G', 1:N, 'F', []);
analysisOptions.Reductions.U.v = 'sequence';
analysisOptions.Reductions.G.v = 'sequence';
analysisOptions.Reductions.F.v = 'sequence';
scheme.setObservation('plot', plotOptions);
scheme.setObservation('analysis', analysisOptions);


% Build velocity discretization (SumSpace)
if config.nu > 0
    vMacroElement = mdemat.common.fem.L2SphereElement.modal(config.nDims, ...
        config.nu, reduction=config.vDimReduction, np=config.nv);
    vMacroDisc = mdemat.common.fem.SpectralSpace(vMacroElement);
else
    vMacroDisc = [];
end

if all(config.nv > 0)
    if config.nDims == 2 && strcmp(config.vDimReduction, 'symmetry')
        nVDims = 3;
    else
        nVDims = config.nDims;
    end
    vMicroElement = mdemat.common.fem.L2SphereElement.nodal(nVDims, ...
        config.nv, reduction=config.vDimReduction);
    vMicroDisc = mdemat.common.fem.SpectralSpace(vMicroElement);
else
    vMicroDisc = [];
end

vDisc = mdemat.common.fem.SumSpace(vMacroDisc, vMicroDisc);

% Build spatial discretization
switch scheme.getConfig('xBasisType')
    case 'modal'
        xElement = mdemat.common.fem.BH1OrthotopeElement.modal(config.nDims, ...
            scheme.getConfig('xBasisOrder'), ...
            pattern=scheme.getConfig('xBasisPattern'));
    case 'nodal'
        xElement = mdemat.common.fem.BH1OrthotopeElement.nodal(config.nDims, ...
            scheme.getConfig('xBasisOrder'));
end
xMesh = mdemat.data.mesh.UniformGrid(config.xBBox, config.nx);
xDisc = mdemat.common.fem.FiniteElementSpace(xMesh, xElement);

% Create state and run simulation
state = mdemat.application.radiation.MacroMicroState(xDisc, vDisc);
[state, result] = scheme.converge(state);
end

%% HELPER FUNCTIONS - FREE STREAMING

function h = fInitFs(config)
h = @(x, v) fInitFsImpl(x, v, config);
end

function f = fInitFsImpl(x, v, config)
E = config.E;
omega = config.omega;
z = omega(:).' * x;
f = (1 + 0.05 * cos(z)) .* E;
end

function h = fExactFs(config)
h = @(x, v, t) fExactFsImpl(x, v, t, config);
end

function f = fExactFsImpl(x, v, t, config)
E = config.E;
omega = config.omega;
z = omega(:).' * (x - v .* t);
f = (1 + 0.05 * cos(z)) .* E;
end

%% HELPER FUNCTIONS - DIFFUSIVE SCALING

function h = fInitDs(config)
h = @(x, v) fInitDsImpl(x, v, config);
end

function f = fInitDsImpl(x, v, config)
E = config.E;
omega = config.omega;
epsilon = config.epsilon;
sigma_s = scatteringDsImpl(x, config);
z = omega(:).' * x;
y = omega(:).' * v;
f = (sin(z) - epsilon * y .* cos(z) ./ sigma_s) .* E;
end

function h = fExactDs(config)
h = @(x, v, t) fExactDsImpl(x, v, t, config);
end

function f = fExactDsImpl(x, v, t, config)
E = config.E;
omega = config.omega;
lambda = config.lambda;
epsilon = config.epsilon;
sigma_s = scatteringDsImpl(x, config);
z = omega(:).' * x;
y = omega(:).' * v;
f = exp(-lambda*t) * (sin(z) - epsilon * y .* cos(z) ./ sigma_s) .* E;
end

function h = scatteringDs(config)
h = @(x) scatteringDsImpl(x, config);
end

function C = scatteringDsImpl(x, config)
C = ones(1, size(x, 2));
end

function h = absorptionDs(config)
h = @(x) absorptionDsImpl(x, config);
end

function C = absorptionDsImpl(x, config)
C = zeros(1, size(x, 2));
end

function h = sourceDs(config)
h = @(x, v, t) sourceDsImpl(x, v, t, config);
end

function S = sourceDsImpl(x, v, t, config)
E = config.E;
omega = config.omega;
lambda = config.lambda;
epsilon = config.epsilon;
alpha = config.timeScale;
beta = config.scatteringScale;
gamma = config.absorptionScale;
z = omega(:).' * x;
y = omega(:).' * v;
ea = epsilon^alpha;
eb = epsilon^beta;
eg = epsilon^gamma;
sigma_s = scatteringDsImpl(x, config);
sigma_a = absorptionDsImpl(x, config);

% \partial_t f
T1 = -lambda * exp(-lambda*t) * (sin(z) - epsilon * y .* cos(z) ./ sigma_s) .* E;
% v \cdot \nabla_x f
T2 = exp(-lambda*t) * y .* (cos(z) + epsilon * y .* sin(z) ./ sigma_s) .* E;
% -\sigma_s * g
T3 = -epsilon * exp(-lambda*t) * y .* cos(z) .* E;
% \sigma_a * f
T4 = exp(-lambda*t) * sigma_a .* (sin(z) - epsilon * y .* cos(z) ./ sigma_s) .* E;

S = (ea * T1 + T2 + eb * T3 + eg * T4) ./ eg;
end
