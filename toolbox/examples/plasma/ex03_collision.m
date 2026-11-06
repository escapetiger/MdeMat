function [state, result] = ex03_collision(options)
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
% FileName: ex03_collision
% Author: Yi CAI
% Date: 27/01/2026
% Description:
%   Collision comparison.
%==========================================================================


%% EXECUTION
% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = mdemat.tools.environment.readConfig(common);

% Checkpoint information
config.id = 'ex03ldh';
config.ckptDir = fullfile(fileparts(mfilename('fullpath')), config.id);
config.ckptPrefix = [upper(config.schemeName), num2str(config.xBasisOrder)];
if config.useFilter
    config.ckptPrefix = ['F', config.ckptPrefix];
end
switch config.id
    case 'ex03ldh'
        config.ckptPrefix = strjoin({'LDH', config.ckptPrefix}, '_');
        config.tFinal = 50;
        config.ckptTimeStamps = [4, 8, 16, 32, 40, 50];
    case 'ex03ldi'
        config.ckptPrefix = strjoin({'LDI', config.ckptPrefix}, '_');
        config.tFinal = 70;
        config.ckptTimeStamps = [4, 8, 16, 30, 50, 70];
        % config.tFinal = 400;
        % config.ckptTimeStamps = [];
    case 'ex03tsh'
        config.ckptPrefix = strjoin({'TSH', config.ckptPrefix}, '_');
        config.tFinal = 50;
        config.ckptTimeStamps = [4, 8, 16, 30, 50];
    case 'ex03tsi'
        config.ckptPrefix = strjoin({'TSI', config.ckptPrefix}, '_');
        config.tFinal = 60;
        config.ckptTimeStamps = [8, 16, 30, 60];
end
config.epsilon = 1;
config.vBBox = [-3.2, 3.2];
config.ckptFigureName = 'f';
config.verbose = 2;

if options.Smoke
    config.epsilon = 1;
    config.tau0 = 1;
    [config, scheme, state, result] = run(config, options); %#ok<ASGLU>
    return;
end

override = 1;
lineStyleOptions = mdemat.tools.visual.plotStyle();
fieldStyleOptions = mdemat.tools.visual.plotStyle();
tau0Level = [inf];
tau0 = 10.^tau0Level;


kmax = length(tau0);
for k = 1:kmax
    config.tau0 = tau0(k);
    config.ckptPostfix = sprintf('tau0%.0e', tau0(k));
    % fileName = sprintf('%s_state_%s.mat', config.ckptPrefix, config.ckptPostfix);
    tckpt = 40;
    fileName = sprintf('%s_state_%s_t%dp0.mat', config.ckptPrefix, config.ckptPostfix, tckpt);
    fileName = fullfile(config.ckptDir, fileName);

    if exist(fileName, 'file') && ~override
        loadOptions = options;
        loadOptions.LoadFile = fileName;
        [config, scheme, state, result] = run(config, loadOptions);
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
    if config.verbose == 2
        imagesc(ax, X, V, F.', 'Interpolation', 'bilinear');
    else
        FInf = state.equilibrium(xRef, V);
        imagesc(ax, X, V, (F-FInf).', 'Interpolation', 'bilinear');
    end
    axis(ax, 'xy', 'equal', 'tight');
    xl = xlabel(ax, 'x');
    yl = ylabel(ax, 'v');
    set(ax, 'FontSize', fieldStyleOptions.FontSize);
    set(xl, 'FontSize', fieldStyleOptions.AxisLabelFontSize);
    set(yl, 'FontSize', fieldStyleOptions.AxisLabelFontSize);
    colormap(ax, fieldStyleOptions.ColorMap);
    colorbar(ax);
    ti = title(sprintf('t = %d', tckpt));
    set(ti, 'FontSize', fieldStyleOptions.TitleFontSize);
    caxis([-0.05, 0.65]);
    % fileName = sprintf('%s_f_%s_t%dp0.pdf', config.ckptPrefix, config.ckptPostfix, tckpt);
    fileName = sprintf('%s_f_%s.pdf', config.ckptPrefix, config.ckptPostfix);
    fileName = fullfile(config.ckptDir, fileName);
    exportgraphics(ax, fileName, 'ContentType', 'vector');

    % % Save free energy plot
    % figure(3);
    % ax = gca;
    % hold(ax, 'on');
    % t = state.History.time;
    % potentialEnergy = sqrt(2*state.History.potentialEnergy);
    % freeEnergy = max(1e-35, potentialEnergy);
    % if k <= k0
    %     style = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, k+1, k, kind="line"));
    % else
    %     style = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, k+1, k-k0, kind="scatter"));
    % end
    % line = plot(t, freeEnergy, style{:});
    % set(line, 'LineWidth', 3, 'MarkerSize', 8);
    % peLegs{k} = sprintf('$\\tau_0=10^{%d}$', tau0Level(k));
    % leg = legend(ax, peLegs(1:k));
    % set(leg, 'Location', lineStyleOptions.LegendLocation);
    % set(leg, 'FontSize', lineStyleOptions.LegendFontSize);
    % set(leg, 'Interpreter', 'latex');
    % set(ax, 'YScale', 'log');
    % hold(ax, 'off');
    % xlim([0, config.tFinal]);
    % xl = xlabel('t');
    % ax.XAxis.Exponent = 0;
    % ax.XAxis.TickLabelFormat = '%.0f';
    % ylim([5e-8, 5e2]);
    % yticks(logspace(-18, 2, 21));
    % yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-18, 2, 21), 'Un', 0));
    % set(ax, 'FontSize', fieldStyleOptions.FontSize);
    % set(xl, 'FontSize', fieldStyleOptions.AxisLabelFontSize);
    % ti = title('Potential energy');
    % set(ti, 'FontSize', fieldStyleOptions.TitleFontSize);
    % fileName = sprintf('%s_potential_energy.pdf', config.ckptPrefix);
    % fileName = fullfile(config.ckptDir, fileName);
    % exportgraphics(ax, fileName, 'ContentType', 'vector');
    % 
    % % Save relative entropy plot
    % figure(4);
    % ax = gca;
    % hold(ax, 'on');
    % t = state.History.time;
    % L2FDistance = max(1e-35, state.History.L2FDistance);
    % if k <= k0
    %     style = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, k+1, k, kind="line"));
    % else
    %     style = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, k+1, k-k0, kind="scatter"));
    % end
    % line = plot(t, L2FDistance, style{:});
    % set(line, 'LineWidth', 3, 'MarkerSize', 8);
    % peLegs{k} = sprintf('$\\tau_0=10^{%d}$', tau0Level(k));
    % leg = legend(ax, peLegs(1:k));
    % set(leg, 'Location', lineStyleOptions.LegendLocation);
    % set(leg, 'FontSize', lineStyleOptions.LegendFontSize);
    % set(leg, 'Interpreter', 'latex');
    % set(ax, 'YScale', 'log');
    % hold(ax, 'off');
    % xlim([0, config.tFinal]);
    % xl = xlabel('t');
    % ax.XAxis.Exponent = 0;
    % ax.XAxis.TickLabelFormat = '%.0f';
    % ylim([5e-5, 5e2]);
    % yticks(logspace(-16, 4, 21));
    % yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-16, 4, 21), 'Un', 0));
    % set(ax, 'FontSize', fieldStyleOptions.FontSize);
    % set(xl, 'FontSize', fieldStyleOptions.AxisLabelFontSize);
    % % ti = title('$\|f-f_\infty\|_{L^2(f_\infty^{-1})}$');
    % % set(ti, 'Interpreter','latex', 'FontSize', 16);
    % ti = title('Distance to kinetic equilibrium');
    % set(ti, 'FontSize', fieldStyleOptions.TitleFontSize);
    % fileName = sprintf('%s_relative_entropy.pdf', config.ckptPrefix);
    % fileName = fullfile(config.ckptDir, fileName);
    % exportgraphics(ax, fileName, 'ContentType', 'vector');
    % 
    % % Save trend to macro equilibrium plot
    % figure(5);
    % ax = gca;
    % hold(ax, 'on');
    % t = state.History.time;
    % L2RhoDistance = max(1e-35, state.History.L2RhoDistance);
    % if k <= k0
    %     style = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, k+1, k, kind="line"));
    % else
    %     style = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, k+1, k-k0, kind="scatter"));
    % end
    % line = plot(t, L2RhoDistance, style{:});
    % set(line, 'LineWidth', 3, 'MarkerSize', 8);
    % peLegs{k} = sprintf('$\\tau_0=10^{%d}$', tau0Level(k));
    % leg = legend(ax, peLegs(1:k));
    % set(leg, 'Location', lineStyleOptions.LegendLocation);
    % set(leg, 'FontSize', lineStyleOptions.LegendFontSize);
    % set(leg, 'Interpreter', 'latex');
    % set(ax, 'YScale', 'log');
    % hold(ax, 'off');
    % xlim([0, config.tFinal]);
    % xl = xlabel('t');
    % ax.XAxis.Exponent = 0;
    % ax.XAxis.TickLabelFormat = '%.0f';
    % ylim([5e-8, 5e2]);
    % yticks(logspace(-16, 4, 21));
    % yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-16, 4, 21), 'Un', 0));
    % set(ax, 'FontSize', fieldStyleOptions.FontSize);
    % set(xl, 'FontSize', fieldStyleOptions.AxisLabelFontSize);
    % ti = title('Distance to macroscopic equilibrium');
    % set(ti, 'FontSize', fieldStyleOptions.TitleFontSize);
    % fileName = sprintf('%s_macro_relative_entropy.pdf', config.ckptPrefix);
    % fileName = fullfile(config.ckptDir, fileName);
    % exportgraphics(ax, fileName, 'ContentType', 'vector');
    % 
    % % Save trend to local equilibrium plot
    % figure(6);
    % ax = gca;
    % hold(ax, 'on');
    % t = state.History.time;
    % L2Entropy = max(1e-35, state.History.L2Entropy);
    % if k <= k0
    %     style = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, k+1, k, kind="line"));
    % else
    %     style = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, k+1, k-k0, kind="scatter"));
    % end
    % line = plot(t, L2Entropy, style{:});
    % set(line, 'LineWidth', 3, 'MarkerSize', 8);
    % peLegs{k} = sprintf('$\\tau_0=10^{%d}$', tau0Level(k));
    % leg = legend(ax, peLegs(1:k));
    % set(leg, 'Location', lineStyleOptions.LegendLocation);
    % set(leg, 'FontSize', lineStyleOptions.LegendFontSize);
    % set(leg, 'Interpreter', 'latex');
    % set(ax, 'YScale', 'log');
    % hold(ax, 'off');
    % xlim([0, config.tFinal]);
    % xl = xlabel('t');
    % ax.XAxis.Exponent = 0;
    % ax.XAxis.TickLabelFormat = '%.0f';
    % ylim([5e-5, 5e2]);
    % yticks(logspace(-16, 4, 21));
    % yticklabels(arrayfun(@(x) sprintf('%.0e', x), logspace(-16, 4, 21), 'Un', 0));
    % set(ax, 'FontSize', fieldStyleOptions.FontSize);
    % set(xl, 'FontSize', fieldStyleOptions.AxisLabelFontSize);
    % ti = title('Distance to local equilibrium');
    % set(ti, 'FontSize', fieldStyleOptions.TitleFontSize);
    % fileName = sprintf('%s_micro_relative_entropy.pdf', config.ckptPrefix);
    % fileName = fullfile(config.ckptDir, fileName);
    % exportgraphics(ax, fileName, 'ContentType', 'vector');
end


% Save conservation plot
figure(8);
ax = gca;
hold(ax, 'on');
style1 = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, 1, 1, kind="line"));
style2 = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, 2, 2, kind="line"));
style3 = namedargs2cell(mdemat.tools.visual.curveStyle(lineStyleOptions, 3, 3, kind="line"));
t = state.History.time;
mass = state.History.mass;
momentum = state.History.momentum;
totalEnergy = state.History.totalEnergy;
mass = mass - mass(1);
momentum = momentum - momentum(1);
totalEnergy = totalEnergy - totalEnergy(1);
line = plot(t, mass, style1{:});
set(line, 'LineWidth', 2);
line = plot(t, momentum, style2{:});
set(line, 'LineWidth', 2);
line = plot(t, totalEnergy, style3{:});
set(line, 'LineWidth', 2);
distLegs = cell(1, 3);
distLegs{1} = 'mass';
distLegs{2} = 'momentum';
distLegs{3} = 'total energy';
leg = legend(ax, distLegs);
set(leg, 'Location', lineStyleOptions.LegendLocation);
set(leg, 'FontSize', lineStyleOptions.LegendFontSize);
set(leg, 'Interpreter', 'latex');
xlim([0, config.tFinal]);
xl = xlabel('$t$');
set(xl, 'Interpreter', 'latex');
ylim([-1e-5, 1e-5]);
yticks(-1e-5:5e-6:1e-5);
yticklabels(arrayfun(@(x) sprintf('%.0e', x), -1e-5:5e-6:1e-5, 'Un', 0));
title('Conservation');
fileName = sprintf('%s_conservation.pdf', config.ckptPrefix);
fileName = fullfile(config.ckptDir, fileName);
exportgraphics(ax, fileName, 'ContentType', 'vector');

%% SIMULATION


if isempty(result) && ~isempty(state)
    result = struct('schemaVersion', 1, 'completed', true, 'diagnostics', state.History);
end

end

function [config, scheme, state, result] = run(config, options)

config.nDims = 1; % Dimension
switch config.id
    case 'ex03ldh'
        config.L = 4*pi; % Length of domain
        config.kappa = 2*pi / config.L; % Wave length
        config.xBBox = [0, config.L]; % Spatial bounding box
        config.delta = 0.01; % Perturbation
    case 'ex03ldi'
        config.L = 6; % Length of domain
        config.kappa = pi / config.L; % Wave length
        config.xBBox = [-config.L, config.L]; % Spatial bounding box
        config.delta = 0.01; % Perturbation
    case 'ex03tsh'
        config.L = 4*pi; % Length of domain
        config.kappa = 2*pi / config.L; % Wave length
        config.xBBox = [0, config.L]; % Spatial bounding box
        config.delta = 0.01; % Perturbation
        config.alpha1 = 1 / 6;
        config.alpha2 = 5 / 6;
    case 'ex03tsi'
        config.L = 6; % Length of domain
        config.kappa = pi / config.L; % Wave length
        config.xBBox = [-config.L, config.L]; % Spatial bounding box
        config.T0 = 1; % Temperature parameter
        config.delta = 0.01; % Perturbation
        config.alpha1 = 1 / 6;
        config.alpha2 = 5 / 6;
end
config.T0 = 1; % Temperature parameter
config.nx = 128; % Grid resolution
config.nh = 50; % Number of Hermite modes - 1
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
switch config.id
    case {'ex03ldh', 'ex03ldi'}
        D0 = (rhoInf + delta * cos(kappa*x)) ./ sqrt(rhoInf);
        D = [D0; zeros(nh-1, size(x, 2))];
    case 'ex03tsh'
        alpha1 = config.alpha1;
        alpha2 = config.alpha2;
        tmp = (cos(2*kappa*x) + cos(3*kappa*x)) * alpha2 + cos(kappa*x);
        D0 = (alpha1 + alpha2) * (1 + delta * tmp) ;
        D2 = sqrt(2) * alpha2 * (1 + delta * tmp);
        if nh > 2
            D = [D0; zeros(1, size(x, 2)); D2; zeros(nh-3, size(x, 2))];
        else
            error('nh must be larger than 2.');
        end
    case 'ex03tsi'
        alpha1 = config.alpha1;
        alpha2 = config.alpha2;
        tmp = cos(kappa*x);
        D0 = (alpha1 + alpha2) * (1 + delta * tmp) ./ sqrt(rhoInf);
        D2 = sqrt(2) * alpha2 * (1 + delta * tmp) ./ sqrt(rhoInf);
        if nh > 2
            D = [D0; zeros(1, size(x, 2)); D2; zeros(nh-3, size(x, 2))];
        else
            error('nh must be larger than 2.');
        end
end
end

function h = phiInf(config)
h = @(x) phiInfImpl(config, x);
end

function phi = phiInfImpl(config, x)
switch config.id
    case {'ex03ldh', 'ex03tsh'}
        phi = zeros(1, size(x, 2));
    case 'ex03ldi'
        kappa = config.kappa;
        phi = 0.2*sin(kappa*x);
    case 'ex03tsi'
        kappa = config.kappa;
        phi = 0.1*(1-cos(kappa*x));
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
