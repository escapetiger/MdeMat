%==========================================================================
% FileName: ex01_landau_damping_convergence
% Author: Yi Cai
% Date: 2025-08-31
% Description:
%   Landau damping convergence test.
%==========================================================================

%% EXECUTION
clc, clear;
dataPath = fullfile(fileparts(mfilename('fullpath')), 'ex03');

tFinal = [0.8, 1.6, 3.0, 6.0];
for i = 1:4
    for j = 2:4
        for k = 1
            [config, state] = simulate(10^(-k), 10^(j), tFinal(i));
            fileName = sprintf('%s_%d.mat', config.SID, k);
            state.save(fullfile(dataPath, fileName));
        end
        ax = gca;
        fileName = sprintf('figure_TS%d_%d_%d.pdf', config.xBasisOrder, j, i);
        delete(ax.Children(2:end))
        exportgraphics(ax,fullfile(dataPath, fileName),'ContentType','vector');
    end
end

%% SIMULATION

function [config, state] = simulate(epsilon, tau0, tFinal)
close all;
configPath = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = core.io.OptionParser().parse(configPath);
config.nDims = 1; % Dimension
config.L = 6; % Length
config.xBBox = [-config.L, config.L]; % Spatial bounding box
config.tFinal = tFinal; % Final time
config.k = pi / config.L; % Wave length
config.T0 = 1; % Temperature parameter
config.alpha1 = 1/6;
config.alpha2 = 5/6;
config.delta = 0.01; % Perturbation
config.phiInf = phiInf(config.k); % Potential steady state
config.EInf = EInf(config.k); % Electrical field steady state
config.cInf = 2 * config.L / integral(@(x) exp(-config.phiInf(x)), -config.L, config.L); % Normalization parameter
config.rhoInf = rhoInf(config.k, config.cInf); % Density steady state
config.sqrtRhoInf = @(x) sqrt(config.rhoInf(x));
config.invRhoInf = @(x) 1 ./ config.rhoInf(x);
config.invSqrtRhoInf = @(x) 1 ./ sqrt(config.rhoInf(x));
config.nx = repmat(128, 1, config.nDims); % Grid resolution
config.nh = 200; % Number of Hermite modes - 1
config.ic.D = DInit(config.k, config.alpha1, config.alpha2, config.delta, config.cInf, config.nh); % Initial hermite coefficients
config.ic.omega = omegaInit(config.k, config.alpha1, config.alpha2, config.delta, config.cInf); % Initial potential
config.bc = []; % Boundary condition
config.exact = [];
config.tau0 = tau0; %
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

config.visualizer = physics.Visualizer(config.nDims);
config.visualizer.setDensity(ones(1, config.nDims));
config.visualizer.setTimeline(config.tFinal, 10);
config.visualizer.setComponents(struct('D', [], 'omega', []));
config.visualizer.setTitleHead(sprintf('%s with %s', config.EID, config.SID));
config.visualizer.addDataset(config.XID, 2);
config.visualizer.addPlotter('profile', '1d');

scheme = physics.vlasov.HermiteScheme(config);

xElement = approx.element.L2OrthotopeElement(config.nDims, ...
    config.xBasisOrder, config.xBasisType, config.xBasisPattern);
xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx);
xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
state = physics.vlasov.HermiteState(xDisc, config.nh);
state = scheme.run(state);
end

%% HELPER FUNCTIONS

function h = fInit(k, delta)
h = @(x, v) fInitImpl(x, v, k, delta);
end

function f = fInitImpl(x, v, k, delta)
M = exp(-abs(v).^2/(2)) / sqrt(2*pi);
f = (rhoInf(x) + delta * cos(k*x)) .* M;
end

function h = DInit(k, alpha1, alpha2, delta, cInf, nh)
h = @(x) DInitImpl(x, k, alpha1, alpha2, delta, cInf, nh);
end

function D = DInitImpl(x, k, alpha1, alpha2, delta, cInf, nh)
rhoInf = rhoInfImpl(x, k, cInf);
O = zeros(1, size(x, 2));
if nh >= 0
D0 = (alpha1 + alpha2) * (rhoInf + delta * cos(k*x)) ./ rhoInf;
end

if nh >= 2
D2 = sqrt(2) * alpha2 * (rhoInf + delta * cos(k*x)) ./ rhoInf;
end

if nh < 2
    D = [D0, repmat(O, 1, nh)];
else
    D = [D0, O, D2, repmat(O, 1, nh-2)];
end
end

function h = omegaInit(k, alpha1, alpha2, delta, cInf)
h = @(x) omegaInitImpl(x, k, alpha1, alpha2, delta, cInf);
end

function omega = omegaInitImpl(x, k, alpha1, alpha2, delta, cInf)
rhoInf = rhoInfImpl(x, k, cInf);
omega = (alpha1 + alpha2) * delta / k^2 * cos(k*x) .* sqrt(rhoInf);
end

function h = phiInf(k)
h = @(x) phiInfImpl(x, k);
end

function phi = phiInfImpl(x, k)
phi = 0.1*(1-cos(k*x));
end

function h = EInf(k)
h = @(x) EInfImpl(x, k);
end

function E = EInfImpl(x, k)
E = -0.1*k*sin(k*x);
end

function h = rhoInf(k, cInf)
h = @(x) rhoInfImpl(x, k, cInf);
end

function rho = rhoInfImpl(x, k, cInf)
rho = cInf * exp(-phiInfImpl(x, k));
end