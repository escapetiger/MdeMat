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
% Date: 13/09/2025
% Description:
%   Convergence test.
%==========================================================================


%% EXECUTION
% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = mdemat.tools.environment.readConfig(common);

% Checkpoint information
config.id = 'ex01ldh';
config.ckptDir = fullfile(fileparts(mfilename('fullpath')), config.id);
config.ckptPrefix = [upper(config.schemeName), num2str(config.xBasisOrder)];
if config.useFilter
    config.ckptPrefix = ['F', config.ckptPrefix];
end
switch config.id
    case 'ex01ldh'
        config.ckptPrefix = strjoin({'LDH', config.ckptPrefix}, '_');
    case 'ex01ldi'
        config.ckptPrefix = strjoin({'LDI', config.ckptPrefix}, '_');
end
config.tFinal = 1;
config.vBBox = [-3.2, 3.2];
config.ckptTimeStamps = [];

if options.Smoke
    config.epsilon = 1;
    config.tau0 = 1;
    [config, scheme, state, result] = run(config, options); %#ok<ASGLU>
    return;
end

% Test the convergence rate from short to long time scale
epsilon = 10.^0;
tau0 = 10.^[1,3,5];
for k1 = 1:length(epsilon)
    for k2 = 1:length(tau0)
        config.epsilon = epsilon(k1);
        config.tau0 = tau0(k2);
        % config.tFinal = epsilon(k1);
        config.ckptPostfix = sprintf('epsilon%.0e_tau0%.0e', epsilon(k1), tau0(k2));
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

config.nDims = 1; % Dimension
switch config.id
    case 'ex01ldh'
        config.L = 4*pi; % Length of domain
        config.kappa = 2*pi / config.L; % Wave length
        config.xBBox = [0, config.L]; % Spatial bounding box
    case 'ex01ldi'
        config.L = 6; % Length of domain
        config.kappa = pi / config.L; % Wave length
        config.xBBox = [-config.L, config.L]; % Spatial bounding box
end
config.T0 = 1; % Temperature parameter
config.delta = 0.05; % Perturbation
config.nx = 16; % Grid resolution
config.nh = 16; % Number of Hermite modes - 1
config.tau = config.tau0 * config.epsilon;
config.phiInf = phiInf(config); % Electric potential equilibrium
config.EInf = EInf(config); % Electric field equilibrium
config.cInf = cInf(config); % Normalization constant
config.rhoInf = rhoInf(config); % Density equilibrium
config.sqrtRhoInf = sqrtRhoInf(config); % Square root of density equilibrium
config.dSqrtRhoInf = dSqrtRhoInf(config); % Derivative of square root of density equilibrium
config.ic = DInit(config); % Initial hermite coefficients
config.bc = []; % Boundary condition
config.exact = []; % Exact solution
config.verbose = 1; % Verbose flag

config = mdemat.internal.execution.studyConfig(config, options);

% Construct and configure scheme
scheme = mdemat.application.plasma.ApDghsScheme(config = config);

% Configure visualizer
plotOptions = scheme.Observation.Plot;
plotOptions.PointsPerDimension = 1;
plotOptions.Times = linspace(0, scheme.getConfig('tFinal'), 1+1);
plotOptions.Components = struct('D', [], 'Psi', [], 'Q', [], 'F', 1);
plotOptions.DatasetTypes.(scheme.getConfig('sId')) = 'numeric';
if options.Plot, plotOptions.Modes{end+1} = '2d'; end

% Configure analyzer
analysisOptions = scheme.Observation.Analysis;
if options.Smoke, analysisOptions.NLevels = 1; else, analysisOptions.NLevels = 6; end
analysisOptions.PointsPerDimension = 4;
analysisOptions.Components = struct('D', 'all', 'Psi', 1, 'Q', 1);
analysisOptions.Reductions.D.v = 'sequence';
scheme.setObservation('plot', plotOptions);
scheme.setObservation('analysis', analysisOptions);


vElement = mdemat.common.fem.L2OrthotopeElement.hermite(1, config.nh, T = config.T0);
vDisc = mdemat.common.fem.SpectralSpace(vElement);

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
state = mdemat.application.plasma.ApHermiteState(xDisc, vDisc);
[state, result] = scheme.converge(state);
end

%% HELPER FUNCTIONS

function h = DInit(config)
h = @(x, varargin) DInitImpl(config, x, varargin{:});
end

function D = DInitImpl(config, x, varargin)
if isempty(varargin), nh = config.nh; else, nh = varargin{1}; end
kappa = config.kappa;
delta = config.delta;
rhoInf = rhoInfImpl(config, x);
D0 = (rhoInf + delta * cos(kappa*x)) ./ sqrt(rhoInf);
D = [D0; zeros(nh-1, size(x, 2))];
end

function h = phiInf(config)
h = @(x) phiInfImpl(config, x);
end

function phi = phiInfImpl(config, x)
switch config.id
    case 'ex01ldh'
        phi = zeros(1, size(x, 2));
    case 'ex01ldi'
        kappa = config.kappa;
        phi = 0.2*sin(kappa*x);
end
end

function h = rhoInf(config)
h = @(x) rhoInfImpl(config, x);
end

function rho = rhoInfImpl(config, x)
cInf = config.cInf;
T0 = config.T0;
rho = cInf * exp(-phiInfImpl(config, x) / T0);
end

function h = sqrtRhoInf(config)
h = @(x) sqrtRhoInfImpl(config, x);
end

function rho = sqrtRhoInfImpl(config, x)
rho = sqrt(rhoInfImpl(config, x));
end

function h = EInf(config)
% Automatic symbolic differentiation of sqrt(rhoInf)
syms xs;
phiSym = -phiInfImpl(config, xs);
ESym = diff(phiSym, xs);

% Check if derivative is a constant (independent of xs)
if ~has(ESym, xs)
    % Derivative is constant - convert to numeric value
    constVal = double(ESym);
    % Return function that returns constant with same size as input
    h = @(x) constVal * ones(1, size(x, 2));
else
    % Derivative depends on x - convert symbolic expression to function handle
    h = matlabFunction(ESym, 'Vars', xs);
end
end

function h = dSqrtRhoInf(config)
% Automatic symbolic differentiation of sqrt(rhoInf)
syms xs;
cInf = config.cInf;
T0 = config.T0;
phiSym = phiInfImpl(config, xs);
sqrtRhoSym = sqrt(cInf) * exp(-phiSym / (2*T0));
dSqrtRhoSym = diff(sqrtRhoSym, xs);

% Check if derivative is a constant (independent of xs)
if ~has(dSqrtRhoSym, xs)
    % Derivative is constant - convert to numeric value
    constVal = double(dSqrtRhoSym);
    % Return function that returns constant with same size as input
    h = @(x) constVal * ones(1, size(x, 2));
else
    % Derivative depends on x - convert symbolic expression to function handle
    h = matlabFunction(dSqrtRhoSym, 'Vars', xs);
end
end

function c = cInf(config)
L = config.L;
T0 = config.T0;
phiInf = config.phiInf;
c = 2 * L / integral(@(x) exp(-phiInf(x)/T0), -L, L);
end
