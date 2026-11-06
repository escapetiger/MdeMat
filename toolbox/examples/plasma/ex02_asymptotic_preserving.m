function [state, result] = ex02_asymptotic_preserving(options)
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
% FileName: ex02_asymptotic_preserving
% Author: Yi CAI
% Date: 27/01/2026
% Description:
%   Asymptotic preserving test.
%==========================================================================


%% EXECUTION
% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = mdemat.tools.environment.readConfig(common);

% Checkpoint information
config.id = 'ex02i';
config.ckptDir = fullfile(fileparts(mfilename('fullpath')), config.id);
config.ckptPrefix = [upper(config.schemeName), num2str(config.xBasisOrder)];
if config.useFilter
    config.ckptPrefix = ['F', config.ckptPrefix];
end
switch config.id
    case 'ex02h'
        config.ckptPrefix = strjoin({'LDH', config.ckptPrefix}, '_');
    case 'ex02i'
        config.ckptPrefix = strjoin({'LDI', config.ckptPrefix}, '_');
end
config.tau0 = 10^5;
config.vBBox = [-3.2, 3.2];
config.ckptFigureName = 'f';

if options.Smoke
    config.epsilon = 1;
    config.tau0 = 1;
    [config, scheme, state, result] = run(config, options); %#ok<ASGLU>
    return;
end

override = 0;
lineStyleOptions = mdemat.tools.visual.plotStyle();
fieldStyleOptions = mdemat.tools.visual.plotStyle();
epsilonLevel = [0,-1,-2,-3,-4,-5,-6];
epsilon = 10.^epsilonLevel;
peLegs = cell(1, length(epsilon));
k0 = inf;
kmax = length(epsilon);
for k = 1:kmax
    config.epsilon = epsilon(k);
    config.tFinal = 100;
    config.ckptTimeStamps = [];
    config.ckptPostfix = sprintf('epsilon%.0e', epsilon(k));
    fileName = sprintf('%s_state_%s.mat', config.ckptPrefix, config.ckptPostfix);
    fileName = fullfile(config.ckptDir, fileName);

    if exist(fileName, 'file') && ~override
        loadOptions = options;
        loadOptions.LoadFile = fileName;
        [config, ~, state, result] = run(config, loadOptions);
    else
        [config, ~, state, result] = run(config, options);
        state.save(fileName, omitCallbacks=true, createDirectory=true, overwrite=logical(override));
    end

    % Save distribution plot
    figure(2);
    ax = gca;
    xRef = 0;
    V = linspace(config.vBBox(1), config.vBBox(2), 512);
    X = state.XDisc.Mesh.collocate(xRef);
    F = state.distribution(xRef, V);
    imagesc(ax, X, V, F.', 'Interpolation', 'bilinear');
    axis(ax, 'xy', 'equal', 'tight');
    xl = xlabel(ax, 'x');
    yl = ylabel(ax, 'v');
    set(xl, 'FontSize', fieldStyleOptions.AxisLabelFontSize);
    set(yl, 'FontSize', fieldStyleOptions.AxisLabelFontSize);
    colormap(ax, fieldStyleOptions.ColorMap);
    colorbar(ax);
    fileName = sprintf('%s_f_%s.pdf', config.ckptPrefix, config.ckptPostfix);
    fileName = fullfile(config.ckptDir, fileName);
    exportgraphics(ax, fileName, 'ContentType', 'vector');

    % Save potential energy plot
    figure(3);
    ax = gca;
    hold(ax, 'on');
    t = state.History.time;
    potentialEnergy = max(1e-35, state.History.potentialEnergy);
    % mask = state.History.time <= config.tFinal;
    % t = state.History.time(mask) / epsilon(k);
    % potentialEnergy = max(1e-35, state.History.potentialEnergy(mask));
    if k <= k0
        style = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, ...
            k+1, k, kind="line"));
    else
        style = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, ...
            k+1, k-k0, kind="scatter"));
    end
    line = plot(t, potentialEnergy, style{:});
    set(line, 'LineWidth', 3, 'MarkerSize', 8);
    set(ax, 'YScale', 'log');
    hold(ax, 'off');
    peLegs{k} = sprintf('$\\varepsilon=10^{%d}$', epsilonLevel(k));
    leg = legend(ax, peLegs(1:k));
    set(leg, 'Location', lineStyleOptions.LegendLocation);
    set(leg, 'FontSize', lineStyleOptions.LegendFontSize);
    set(leg, 'Interpreter', 'latex');
    xl = xlabel('$t$');
    % xl = xlabel('$\varepsilon t$');
    set(xl, 'Interpreter','latex', 'FontSize', 16);
    ax.XAxis.Exponent = 0;
    ax.XAxis.TickLabelFormat = '%.0f';
    xlim([0, config.tFinal]);
    % xlim([0, config.tFinal/epsilon(k)]);
    ylim([1e-32, 1e2]);
    yticks(logspace(-35, 0, 8));
    yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-35, 0, 8), 'Un', 0));
    title('Potential energy', 'FontSize', 16);
    fileName = sprintf('%s_potential_energy.pdf', config.ckptPrefix);
    % fileName = sprintf('%s_rescaled_potential_energy.pdf', config.ckptPrefix);
    fileName = fullfile(config.ckptDir, fileName);
    exportgraphics(ax, fileName, 'ContentType', 'vector');

    % Save relative entropy plot
    figure(4);
    ax = gca;
    hold(ax, 'on');
    t = state.History.time;
    L2FDistance = max(1e-35, state.History.L2FDistance);
    % mask = state.History.time <= config.tFinal;
    % t = state.History.time(mask) / epsilon(k);
    % L2FDistance = max(1e-35, state.History.L2FDistance(mask));
    if k <= k0
        style = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, ...
            k+1, k, kind="line"));
    else
        style = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, ...
            k+1, k-k0, kind="scatter"));
    end
    line = plot(t, L2FDistance, style{:});
    set(line, 'LineWidth', 3, 'MarkerSize', 8);
    peLegs{k} = sprintf('$\\varepsilon=10^{%d}$', epsilonLevel(k));
    leg = legend(ax, peLegs(1:k));
    set(leg, 'Location', lineStyleOptions.LegendLocation);
    set(leg, 'FontSize', lineStyleOptions.LegendFontSize);
    set(leg, 'Interpreter', 'latex');
    set(ax, 'YScale', 'log');
    hold(ax, 'off');
    xlim([0, config.tFinal]);
    % xlim([0, config.tFinal/epsilon(k)]);
    xl = xlabel('$t$');
    % xl = xlabel('$\varepsilon t$');
    set(xl, 'Interpreter','latex', 'FontSize', 16);
    ax.XAxis.Exponent = 0;
    ax.XAxis.TickLabelFormat = '%.0f';
    ylim([5e-12, 2e-1]);
    yticks(logspace(-16, 0, 17));
    yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-16, 0, 17), 'Un', 0));
    ti = title('$\|f-f_\infty\|_{L^2(f_\infty^{-1})}$');
    set(ti, 'Interpreter','latex', 'FontSize', 16);
    fileName = sprintf('%s_relative_entropy.pdf', config.ckptPrefix);
    % fileName = sprintf('%s_rescaled_relative_entropy.pdf', config.ckptPrefix);
    fileName = fullfile(config.ckptDir, fileName);
    exportgraphics(ax, fileName, 'ContentType', 'vector');

    % Save free energy plot
    figure(5);
    ax = gca;
    hold(ax, 'on');
    t = state.History.time;
    freeEnergy = max(1e-35, L2FDistance.^2 / 2 + potentialEnergy);
    % mask = state.History.time <= config.tFinal;
    % t = state.History.time(mask) / epsilon(k);
    % freeEnergy = max(1e-35, freeEnergy(mask));
    if k <= k0
        style = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, ...
            k+1, k, kind="line"));
    else
        style = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, ...
            k+1, k-k0, kind="scatter"));
    end
    line = plot(t, freeEnergy, style{:});
    set(line, 'LineWidth', 3, 'MarkerSize', 8);
    peLegs{k} = sprintf('$\\varepsilon=10^{%d}$', epsilonLevel(k));
    leg = legend(ax, peLegs(1:k));
    set(leg, 'Location', lineStyleOptions.LegendLocation);
    set(leg, 'FontSize', lineStyleOptions.LegendFontSize);
    set(leg, 'Interpreter', 'latex');
    set(ax, 'YScale', 'log');
    hold(ax, 'off');
    xlim([0, config.tFinal]);
    % xlim([0, config.tFinal/epsilon(k)]);
    xl = xlabel('$t$');
    % xl = xlabel('$\varepsilon t$');
    set(xl, 'Interpreter','latex', 'FontSize', 16);
    ax.XAxis.Exponent = 0;
    ax.XAxis.TickLabelFormat = '%.0f';
    ylim([1e-27, 1e2]);
    yticks(logspace(-35, 0, 8));
    yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-35, 0, 8), 'Un', 0));
    ti = title('$\mathcal{E}_h$');
    set(ti, 'Interpreter','latex', 'FontSize', 16);
    fileName = sprintf('%s_free_energy.pdf', config.ckptPrefix);
    % fileName = sprintf('%s_rescaled_free_energy.pdf', config.ckptPrefix);
    fileName = fullfile(config.ckptDir, fileName);
    exportgraphics(ax, fileName, 'ContentType', 'vector');
end



%% SIMULATION


if isempty(result) && ~isempty(state)
    result = struct('schemaVersion', 1, 'completed', true, 'diagnostics', state.History);
end

end

function [config, scheme, state, result] = run(config, options)

config.nDims = 1; % Dimension
switch config.id
    case 'ex02h'
        config.L = 4*pi; % Length of domain
        config.kappa = 2*pi / config.L; % Wave length
        config.xBBox = [0, config.L]; % Spatial bounding box
    case 'ex02i'
        config.L = 6; % Length of domain
        config.kappa = pi / config.L; % Wave length
        config.xBBox = [-config.L, config.L]; % Spatial bounding box
end
config.T0 = 1; % Temperature parameter
config.delta = 0.01; % Perturbation
config.nx = 128; % Grid resolution
config.nh = 640; % Number of Hermite modes - 1
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
scheme.setObservation('plot', plotOptions);


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
if isfield(options, 'LoadFile')
    prototype = scheme.initialize(state);
    state = mdemat.internal.execution.loadStudyState(options.LoadFile, prototype);
    result = struct('schemaVersion', 1, 'completed', true, 'diagnostics', state.History);
elseif options.Smoke
    result = mdemat.internal.execution.compute(scheme, state);
    state = result.state;
else
    state = scheme.run(state);
    result = mdemat.internal.execution.studyResult(scheme, state);
end
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
    case 'ex02h'
        phi = zeros(1, size(x, 2));
    case 'ex02i'
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
