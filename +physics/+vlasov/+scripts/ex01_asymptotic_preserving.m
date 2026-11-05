%==========================================================================
% FileName: ex01_asymptotic_preserving
% Author: Yi Cai
% Date: 2025-08-03
% Description:
%   Diffusive scaling with periodic boundary conditions.
%==========================================================================

clc;
clear;
close all;

%% CONFIGURATION
configPath = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = core.io.OptionParser().parse(configPath);
config.nDims = 1; % Dimension
config.L = 6; % Length
config.xBBox = [-config.L, config.L]; % Spatial bounding box
config.tFinal = 0.1; % Final time
config.k = pi / config.L; % Wave length
config.T0 = 1; % Temperature parameter
config.delta = 0.01; % Perturbation
config.phiInf = phiInf(config.k); % Potential steady state
config.EInf = EInf(config.k); % Electrical field steady state
config.cInf = 2 * config.L / integral(@(x) exp(-config.phiInf(x)), -config.L, config.L); % Normalization parameter
config.rhoInf = rhoInf(config.k, config.cInf); % Density steady state
config.nx = repmat(8, 1, config.nDims); % Grid resolution
config.nh = 2; % Number of Hermite modes
config.ic.D = DInit(config.k, config.delta, config.cInf, config.nh); % Initial hermite coefficients
config.ic.omega = omegaInit(config.k, config.delta, config.cInf); % Initial potential
config.bc = []; % Boundary condition
config.exact = []; 
config.tau0 = 1; % 
config.epsilon = 1e-3; % Scaling parameter
config.timeScale = 1; % Time scale
config.collisionScale = 0; % Collision scale
config.verbose = 1; % Verbose flag
config.experimentId = 'AP1';
config.vDiscId = sprintf('H%d', config.nh);
config.tDiscId = upper(config.tDiscId);
config.xDiscId = sprintf('%s%d', upper(config.xDiscName), config.xBasisOrder);
config.schemeId = strjoin({config.xDiscId, config.vDiscId}, '-');

%% SIMULATION SETUP

scheme = physics.vlasov.HermiteScheme(config) ...
    .setTimeDiscretization(config.tDiscId, config.tFinal) ...
    .setTimer(config.verbose) ...
    .setVisualizer( ...
    'experimentId', config.experimentId, ...
    'schemeId', config.schemeId, ...
    'nDims', config.nDims, ...
    'nTimeNodes', 10, ...
    'final', config.tFinal, ...
    'density', ones(1, config.nDims), ...
    'components', struct('D', [], 'omega', [])) ...
    .addDataset('DG', 2) ...
    .addStyle('DG', {'Color', 'b', 'Marker', 'o', 'LineStyle', 'none'});

%% SIMULATION EXECUTION

xDisc = approx.discretization.FiniteElementDiscretization() ...
    .setDgOrthotopeElement( ...
    'nDims', config.nDims, ...
    'order', config.xBasisOrder, ...
    'basisType', config.xBasisType, ...
    'basisPattern', config.xBasisPattern ...
    ) ...
    .setElementOperator() ...
    .setDerivativeOrder(1) ...
    .setUniformGrid(config.nx, config.xBBox) ...
    .setMeshSpace();

state = physics.vlasov.HermiteState(xDisc, config.nh);
state = scheme.initialize(state);
[state, results] = scheme.run(state);
state = scheme.finalize(state);

figure(2)
semilogy(results.t, results.mass, 'o-')

figure(3)
semilogy(results.t, results.potential, 'o-')

%% HELPER FUNCTIONS

function h = fInit(k, delta)
h = @(x, v) fInitImpl(x, v, k, delta);
end

function f = fInitImpl(x, v, k, delta)
M = exp(-abs(v).^2/(2)) / sqrt(2*pi);
f = (rhoInf(x) + delta*cos(k*x)) .* M;
end

function h = DInit(k, delta, cInf, nh)
h = @(x) DInitImpl(x, k, delta, cInf, nh);
end

function D = DInitImpl(x, k, delta, cInf, nh)
D0 = (rhoInfImpl(x, k, cInf) + delta * cos(k*x)) ./ sqrt(rhoInfImpl(x, k, cInf));
D = [D0, zeros(1, size(x, 2) * nh)];
end

function h = omegaInit(k, delta, cInf)
h = @(x) omegaInitImpl(x, k, delta, cInf);
end

function omega = omegaInitImpl(x, k, delta, cInf)
omega = delta / k^2 * cos(k*x) .* sqrt(rhoInfImpl(x, k, cInf));
end

function h = phiInf(k)
h = @(x) phiInfImpl(x, k); 
end

function phi = phiInfImpl(x, k)
% phi = -0.2*sin(k*x);
phi = ones(1, size(x, 2));
end

function h = EInf(k)
h = @(x) EInfImpl(x, k);
end

function E = EInfImpl(x, k)
% E = 0.2*k*cos(k*x);
E = zeros(1, size(x, 2));
end

function h = rhoInf(k, cInf)
h = @(x) rhoInfImpl(x, k, cInf); 
end

function rho = rhoInfImpl(x, k, cInf)
rho = cInf * exp(-phiInfImpl(x, k));
end