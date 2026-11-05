%==========================================================================
% FileName: ex01_landau_damping_convergence
% Author: Yi Cai
% Date: 2025-08-31
% Description:
%   Landau damping convergence test.
%==========================================================================

%% EXECUTION
clc, clear;
dataPath = fullfile(fileparts(mfilename('fullpath')), 'ex01');

for k = 0
    [config, state, table] = simulate(10^(-k));
    fieldNames = fieldnames(table);
    for iName = 1:length(fieldNames)
        fieldName = fieldNames{iName};
        fileName = sprintf('%s_%s_%d.csv', config.SID, fieldName, k);
        writetable(table.(fieldName), fullfile(dataPath, fileName));
    end
    fileName = sprintf('%s_%d.mat', config.SID, k);
    state.save(fullfile(dataPath, fileName));
end

%% SIMULATION

function [config, state, table] = simulate(epsilon)
close all;
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
config.sqrtRhoInf = @(x) sqrt(config.rhoInf(x));
config.invRhoInf = @(x) 1 ./ config.rhoInf(x);
config.invSqrtRhoInf = @(x) 1 ./ sqrt(config.rhoInf(x));
config.nx = repmat(16, 1, config.nDims); % Grid resolution
config.nh = 63; % Number of Hermite modes - 1
config.ic.D = DInit(config.k, config.delta, config.cInf, config.nh); % Initial hermite coefficients
config.ic.omega = omegaInit(config.k, config.delta, config.cInf); % Initial potential
config.bc = []; % Boundary condition
config.exact = [];
config.tau0 = 1; %
config.epsilon = epsilon; % Scaling parameter
config.tau = config.tau0;
config.verbose = 1; % Verbose flag
config.TID = upper(config.tOdeIntName);
config.XID = sprintf('%s%d', upper(config.xDiscName), config.xBasisOrder);
config.VID = sprintf('H%d', config.nh);
config.EID = sprintf('AP%d', config.nDims); % Experiment ID
config.SID = strjoin({config.XID, config.VID}, '-'); % Scheme ID

config.records = {'Integration', 'Visualization', 'Analysis'};
config.timer = core.chrono.Timer(config.records, config.verbose);

config.visualizer = profilers.Visualizer(config.nDims);
config.visualizer.setDensity(ones(1, config.nDims));
config.visualizer.setTimeline(config.tFinal, 1);
config.visualizer.setComponents(struct('D', 1, 'omega', 1));
config.visualizer.setTitleHead(sprintf('%s with %s', config.EID, config.SID));
config.visualizer.addDataset(config.XID, 2);

switch config.nDims
    case 1
        config.visualizer.addPlotter('profile', '1d');
    case 2
        %         config.visualizer.addPlotter('profile', '2d');
        config.visualizer.addPlotter('slice', 'slice1d');
    case 3
        %         config.visualizer.addPlotter('profile', '3d');
        config.visualizer.addPlotter('slice', 'slice1d');
end

config.analyzer = profilers.Analyzer(6);
config.analyzer.setReductions({'L1', 'L2', 'Lx'});
config.analyzer.setDensity(repmat(4, 1, config.nDims));
config.analyzer.setComponents(struct('D', 1:(config.nh+1), 'omega', 1));

scheme = physics.vlasov.HermiteScheme(config);

xElement = approx.element.L2OrthotopeElement(config.nDims, ...
    config.xBasisOrder, config.xBasisType, config.xBasisPattern);
xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx/2);
xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
state = physics.vlasov.HermiteState(xDisc, config.nh);
[state, table] = scheme.converge(state);
end

%% HELPER FUNCTIONS

function h = fInit(k, delta)
h = @(x, v) fInitImpl(x, v, k, delta);
end

function f = fInitImpl(x, v, k, delta)
M = exp(-abs(v).^2/(2)) / sqrt(2*pi);
f = (rhoInf(x) + delta * cos(k*x)) .* M;
end

function h = DInit(k, delta, cInf, nh)
h = @(x) DInitImpl(x, k, delta, cInf, nh);
end

function D = DInitImpl(x, k, delta, cInf, nh)
rhoInf = rhoInfImpl(x, k, cInf);
D0 = (rhoInf + delta * cos(k*x)) ./ rhoInf;
D = [D0, zeros(1, size(x, 2)*nh)];
end

function h = omegaInit(k, delta, cInf)
h = @(x) omegaInitImpl(x, k, delta, cInf);
end

function omega = omegaInitImpl(x, k, delta, cInf)
rhoInf = rhoInfImpl(x, k, cInf);
omega = delta / k^2 * cos(k*x) .* sqrt(rhoInf);
end

function h = phiInf(k)
h = @(x) phiInfImpl(x, k);
end

function phi = phiInfImpl(x, k)
% phi = -0.2*sin(k*x);
phi = zeros(size(x));
end

function h = EInf(k)
h = @(x) EInfImpl(x, k);
end

function E = EInfImpl(x, k)
% E = 0.2*k*cos(k*x);
E = zeros(size(x));
end

function h = rhoInf(k, cInf)
h = @(x) rhoInfImpl(x, k, cInf);
end

function rho = rhoInfImpl(x, k, cInf)
rho = cInf * exp(-phiInfImpl(x, k));
end