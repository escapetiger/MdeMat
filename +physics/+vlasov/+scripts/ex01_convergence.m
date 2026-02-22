%==========================================================================
% FileName: ex01_convergence
% Author: Yi CAI
% Date: 13/09/2025
% Description:
%   Convergence test.
%==========================================================================

clc, clear, close all;

%% EXECUTION
% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = physics.scheme.ConfigParser.parseFile(common);

% Checkpoint information
config.id = 'ex01ldi';
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
config.tFinal = 0.1;
config.vBBox = [-3.2, 3.2];
config.ckptTimeStamps = [];

% Test the convergence rate from short to long time scale
epsilon = 10.^[0, -1, -2, -3, -4, -5, -6];
tau0 = 10.^[1, 3, 5];
for k1 = 1:length(epsilon)
    for k2 = 1:length(tau0)
        config.epsilon = epsilon(k1);
        config.tau0 = tau0(k2);
        config.tFinal = epsilon(k1);
        config.ckptPostfix = sprintf('epsilon%.0e_tau0%.0e', epsilon(k1), tau0(k2));
        fileName = sprintf('%s_state_%s.mat', config.ckptPrefix, config.ckptPostfix);
        fileName = fullfile(config.ckptDir, fileName);
    
        [config, scheme, state, result] = run(config);
        state.save(fileName);
    
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

function [config, scheme, state, result] = run(config)

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

% Construct and configure scheme
scheme = physics.vlasov.ApDghsScheme(config = config);

% Configure visualizer
visualizer = scheme.getConfig('visualizer');
visualizer.setDensity(1);
visualizer.setTimeline(scheme.getConfig('tFinal'), 1);
visualizer.setComponents(struct('D', [], 'Psi', [], 'Q', []));
visualizer.addDataset(scheme.getConfig('sId'), 'numeric');
visualizer.addPlotter('profile', '1d');

% Configure analyzer
analyzer = scheme.getConfig('analyzer');
analyzer.setNLevels(5);
analyzer.setDensity(4^config.nDims);
analyzer.setComponents(struct('D', 1:config.nh, 'Psi', 1, 'Q', 1));

vElement = approx.element.L2OrthotopeElement.hermite(1, config.nh, T = config.T0);
vDisc = approx.space.SpectralSpace(vElement);

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
state = physics.vlasov.ApHermiteState(xDisc, vDisc);
[state, result] = scheme.converge(state);
end

%% HELPER FUNCTIONS

function h = DInit(config)
h = @(x, varargin) DInitImpl(config, x, varargin{:});
end

function D = DInitImpl(config, x, varargin)
if nargin >= 1, nh = varargin{1}; end
kappa = config.kappa;
delta = config.delta;
rhoInf = rhoInfImpl(config, x);
D0 = (rhoInf + delta * cos(kappa*x)) ./ sqrt(rhoInf);
D = [D0, zeros(1, size(x, 2)*(nh-1))];
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