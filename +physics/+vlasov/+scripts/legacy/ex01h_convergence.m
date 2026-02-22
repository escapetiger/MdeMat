%==========================================================================
% FileName: ex01h_convergence
% Author: Yi CAI
% Date: 13/09/2025
% Description:
%   Homogeneous Landau damping convergence test.
%==========================================================================

clc, clear, close all;

%% EXECUTION
% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = physics.scheme.ConfigParser.parseFile(common);

% Checkpoint information
config.ckptDir = fullfile(fileparts(mfilename('fullpath')), 'ex01h');
config.ckptPrefix = [upper(config.schemeName), num2str(config.xBasisOrder)];
if config.useFilter
    config.ckptPrefix = ['F', config.ckptPrefix];
end
config.ckptPrefix = strjoin({'LDH', config.ckptPrefix}, '_');
config.tFinal = 0.1;
config.ckptTimeStamps = [];

% Test the convergence rate from weakly to strongly collisional regime
tau0 = 10.^[6];
for k = 1:length(tau0)
    config.epsilon = 1;
    config.tau0 = tau0(k);
    config.ckptPostfix = sprintf('tau0%.0e', tau0(k));
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

% Test the convergence rate from short to long time scale
% epsilon = 10.^[0, -2, -4];
% for k = 1:length(epsilon)
%     config.epsilon = epsilon(k);
%     config.tau0 = 1;
%     config.tFinal = epsilon(k);
%     config.ckptPostfix = sprintf('epsilon%.0e', epsilon(k));
%     fileName = sprintf('%s_state_%s.mat', config.ckptPrefix, config.ckptPostfix);
%     fileName = fullfile(config.ckptDir, fileName);
% 
%     [config, scheme, state, result] = run(config);
%     state.save(fileName);
% 
%     fieldNames = fieldnames(result);
%     for iName = 1:length(fieldNames)
%         fieldName = fieldNames{iName};
%         fileName = sprintf('%s_%s_%s.csv', config.ckptPrefix, fieldName, config.ckptPostfix);
%         fileName = fullfile(config.ckptDir, fileName);
%         writetable(result.(fieldName), fileName);
%     end
% end

%% SIMULATION

function [config, scheme, state, result] = run(config)

config.nDims = 1; % Dimension
config.kappa = 0.5; % Wave length
config.L = 2 * pi / config.kappa; % Length of domain
config.xBBox = [0, config.L]; % Spatial bounding box
config.vBBox = [-6, 6]; % Velocity bounding box
config.T0 = 1; % Temperature parameter
config.delta = 0.05; % Perturbation
config.nx = 16; % Grid resolution
config.nh = 8; % Number of Hermite modes - 1
config.ic = DInit(config.kappa, config.delta); % Initial hermite coefficients
config.bc = []; % Boundary condition
config.exact = []; % Exact solution
config.verbose = 2; % Verbose flag

% Construct and configure scheme
switch config.schemeName
    case 'dghs'
        scheme = physics.vlasov.DghsScheme(config = config);
    case 'cdghs'
        scheme = physics.vlasov.CdghsScheme(config = config);
end

% Configure visualizer
visualizer = scheme.getConfig('visualizer');
visualizer.setDensity(1);
visualizer.setTimeline(scheme.getConfig('tFinal'), 1);
visualizer.setComponents(struct('D', [], 'P', [], 'E', []));
visualizer.addDataset(scheme.getConfig('sId'), 'numeric');
visualizer.addPlotter('profile', '1d');

% Configure analyzer
analyzer = scheme.getConfig('analyzer');
analyzer.setNLevels(5);
analyzer.setDensity(4^config.nDims);
analyzer.setComponents(struct('D', 1:config.nh, 'P', 1, 'E', 1));

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
state = physics.vlasov.HermiteState(xDisc, vDisc);
[state, result] = scheme.converge(state);
end

%% HELPER FUNCTIONS

function h = DInit(kappa, delta)
h = @(x, varargin) DInitImpl(kappa, delta, x, varargin{:});
end

function D = DInitImpl(kappa, delta, x, varargin)
if nargin >= 1, nh = varargin{1}; end
D0 = (1 + delta * cos(kappa*x));
D = [D0, zeros(1, size(x, 2)*(nh-1))];
end