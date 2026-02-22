%==========================================================================
% FileName: ex03_gaussian
% Author: Yi CAI
% Description:
%   Gaussian pulse test for radiation transport. An isotropic Gaussian
%   density pulse propagates outward through a pure scattering medium.
%==========================================================================

clc, clear, close all;

%% EXECUTION

% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = physics.scheme.ConfigParser.parseFile(common);

% Problem setup
config.id = 'ex03_2dt';

% Parse problem ID
tok = regexp(config.id, '_(\d+)d([ts])$', 'tokens');
assert(~isempty(tok), 'String format not recognized.');
config.nDims = str2double(tok{1}{1});
switch tok{1}{2}
    case 't'
        config.vDimReduction = 'topology';
    case 's'
        config.vDimReduction = 'symmetry';
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
config.L = 1;
config.xBBox = repmat([-config.L, config.L], 1, config.nDims);
config.nx = repmat(64, 1, config.nDims);

% Physical parameters
config.decomposition = '';
config.timeScale = 1;
config.scatteringScale = -1;
config.absorptionScale = 1;

% Scattering and absorption
config.scattering = fScattering(config);
config.absorption = [];
config.source = [];

% Initial and boundary conditions
config.ic = fInit(config);
config.bc = fBc(config);
config.exact = [];
config.filtering = true;

% Velocity discretization parameters based on dimensions
switch config.nDims
    case 1
        if strcmp(config.vDimReduction, 'topology')
            config.nu = 1;
            config.nv = 2;
        else
            config.nu = 2;
            config.nv = 16;
        end
    case 2
        if strcmp(config.vDimReduction, 'topology')
            config.nu = 5;
            config.nv = 8;
        else
            config.nu = 2;
            config.nv = 16;
        end
end

% Checkpoint information
config.ckptDir = fullfile(fileparts(mfilename('fullpath')), config.id);
if ~isfolder(config.ckptDir), mkdir(config.ckptDir); end
config.ckptTimeStamps = [0.1, 0.2, 0.3, 0.4];
config.verbose = 2;
config.tFinal = 0.4;

% Parameter sweep
override = 1;
order = [2];
epsilon = 10.^[0];
for k = 1:length(epsilon)
    plt = struct();
    plt.legend = {};
    plt.time = 0.2;
    config.epsilon = epsilon(k);
    config.ckptPostfix = sprintf('epsilon%.0e', epsilon(k));

    %< Reference
    config.ckptPrefix = [upper(config.schemeName), num2str(1)];
    switch config.id
        case 'ex03_2dt'
            config.ckptPrefix = strjoin({'GS2DT', config.ckptPrefix}, '_');
    end
    if ~isempty(plt.time) && ~override
        timeStr = sprintf('t%.1f', plt.time);
        timeStr = strrep(timeStr, '.', 'p');
        fileName = sprintf('%s_state_%s_%s.mat', config.ckptPrefix, config.ckptPostfix, timeStr);
    else
        fileName = sprintf('%s_state_%s.mat', config.ckptPrefix, config.ckptPostfix);
    end
    fileName = fullfile(config.ckptDir, fileName);

    if exist(fileName, 'file') && ~override
        state = physics.radiation.MacroMicroState.load(fileName);
        plt.style = 'line';
        plt.legend{end+1} = 'REF';
        plt.colorIdx = 1;
        plt.markerIdx = 1;
        switch state.NXDims
            case 1
                plotDensity1d(config, state, plt);
            case 2
                plotDensity2d(config, state, plt);
        end
    end

    for j = 1:length(order)
        %< Numeric
        switch order(j)
            case 1
                config.tOdeIntName = 'be';
                config.xBasisOrder = 1;
            case 2
                config.tOdeIntName = 'ssprk2';
                config.xBasisOrder = 2;
            case 3
                config.tOdeIntName = 'ssprk3';
                config.xBasisOrder = 3;
        end
        config.ckptPrefix = [upper(config.schemeName), num2str(config.xBasisOrder)];
        switch config.id
            case 'ex03_2dt'
                config.ckptPrefix = strjoin({'GS2DT', config.ckptPrefix}, '_');
        end

        if ~isempty(plt.time) && ~override
            timeStr = sprintf('t%.1f', plt.time);
            timeStr = strrep(timeStr, '.', 'p');
            fileName = sprintf('%s_state_%s_%s.mat', config.ckptPrefix, config.ckptPostfix, timeStr);
        else
            fileName = sprintf('%s_state_%s.mat', config.ckptPrefix, config.ckptPostfix);
        end
        fileName = fullfile(config.ckptDir, fileName);

        if exist(fileName, 'file') && ~override
            state = physics.radiation.MacroMicroState.load(fileName);
        else
            [config, scheme, state] = run(config);
            state.save(fileName);
        end

        if config.verbose >= 2
            plt.style = 'scatter';
            plt.legend{end+1} = sprintf('MMDG%d', order(j));
            plt.colorIdx = j+1;
            plt.markerIdx = j;
            switch state.NXDims
                case 1
                    plotDensity1d(config, state, plt);
                case 2
                    plotDensity2d(config, state, plt);
            end
        end
    end
end

%% SIMULATION

function [config, scheme, state] = run(config)

% Construct and configure scheme
switch config.schemeName
    case 'mmdg'
        scheme = physics.radiation.MmDgScheme(config=config);
end

% Configure visualizer
visualizer = scheme.getConfig('visualizer');
visualizer.setDensity(1);
visualizer.setTimeline(scheme.getConfig('tFinal'), 10);
visualizer.setComponents(struct('U', 1, 'G', []));
visualizer.addDataset(scheme.getConfig('sId'), 'numeric');
switch config.nDims
    case 1
        visualizer.addPlotter('profile', '1d');
    case 2
        visualizer.addPlotter('profile', '2d');
        visualizer.addPlotter('slice', 'slice1d');
end

% Build velocity discretization (SumSpace)
if config.nu > 0
    vMacroElement = approx.element.L2SphereElement.modal(config.nDims, ...
        config.nu, reduction=config.vDimReduction);
    vMacroDisc = approx.space.SpectralSpace(vMacroElement);
else
    vMacroDisc = [];
end

if all(config.nv > 0)
    if config.nDims == 2 && strcmp(config.vDimReduction, 'symmetry')
        nVDims = 3;
    else
        nVDims = config.nDims;
    end
    vMicroElement = approx.element.L2SphereElement.nodal(nVDims, ...
        config.nv, reduction=config.vDimReduction);
    vMicroDisc = approx.space.SpectralSpace(vMicroElement);
else
    vMicroDisc = [];
end

vDisc = approx.space.SumSpace(vMacroDisc, vMicroDisc);

% Build spatial discretization
switch scheme.getConfig('xBasisType')
    case 'modal'
        xElement = approx.element.BH1OrthotopeElement.modal(config.nDims, ...
            scheme.getConfig('xBasisOrder'), ...
            pattern=scheme.getConfig('xBasisPattern'));
    case 'nodal'
        xElement = approx.element.BH1OrthotopeElement.nodal(config.nDims, ...
            scheme.getConfig('xBasisOrder'));
end
xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx);
xDisc = approx.space.FiniteElementSpace(xMesh, xElement);

% Create state and run simulation
state = physics.radiation.MacroMicroState(xDisc, vDisc);
state = scheme.run(state);
end

%% HELPER FUNCTIONS

function h = fInit(config)
h = @(x, v) fInitImpl(x, v, config);
end

function f = fInitImpl(x, v, config)
E = config.E;
sigma = sqrt(1e-2);
f = exp(-sum(x.^2, 1) / (2 * sigma^2)) * E;
end

function h = fBc(config)
h = @(x, v, t) fBcImpl(x, v, t, config);
end

function f = fBcImpl(x, v, t, config)
f = zeros(1, size(x, 2));
end

function h = fScattering(config)
h = @(x) fScatteringImpl(x, config);
end

function C = fScatteringImpl(x, config)
C = ones(1, size(x, 2));
end


function plotDensity1d(config, state, plt)
strategy1d = physics.visual.Strategy1d();

figure(2);
ax = gca;

x = state.XDisc.Mesh.collocate(0);
density = state.XDisc.eval(0, state.Dofs.U(:, 1));

hold(ax, 'on');
switch plt.style
    case 'line'
        style = strategy1d.getDefaultLineStyle(plt.colorIdx, plt.markerIdx);
    case 'scatter'
        style = strategy1d.getDefaultScatterStyle(plt.colorIdx, plt.markerIdx);
end
line = plot(x, density, style{:});
set(line, 'LineWidth', 3, 'MarkerSize', 8);
leg = legend(ax, plt.legend);
set(leg, 'Location', strategy1d.DefaultLegendPosition);
set(leg, 'FontSize', strategy1d.DefaultLegendFontSize);
set(leg, 'Interpreter', 'latex');
hold(ax, 'off');
xlim([0, config.L]);
xl = xlabel('x');
set(ax, 'FontSize', strategy1d.DefaultFontSize);
set(xl, 'FontSize', strategy1d.DefaultAxisLabelFontSize);
ti = title('Density');
set(ti, 'FontSize', strategy1d.DefaultTitleFontSize);
if ~isempty(plt.time)
    timeStr = sprintf('t%.1f', plt.time);
    timeStr = strrep(timeStr, '.', 'p');
    fileName = sprintf('%s_density_%s_%s.pdf', config.ckptPrefix, config.ckptPostfix, timeStr);
else
    fileName = sprintf('%s_density_%s.pdf', config.ckptPrefix, config.ckptPostfix);
end
fileName = fullfile(config.ckptDir, fileName);
exportgraphics(ax, fileName, 'ContentType', 'vector');
end

function plotDensity2d(config, state, plt)
strategy2d = physics.visual.Strategy2d();

figure(2);
ax = gca;
density = state.XDisc.eval(zeros(state.NXDims, 1), state.Dofs.U(:, 1));
x = state.XDisc.Mesh.collocate(repmat({0}, state.NXDims, 1));
density = reshape(density, length(x{1}), length(x{2}));
imagesc(ax, x{1}, x{2}, density.', 'Interpolation', 'bilinear');
axis(ax, 'xy', 'equal', 'tight');
xl = xlabel(ax, 'x');
yl = ylabel(ax, 'y');
set(ax, 'FontSize', strategy2d.DefaultFontSize);
set(xl, 'FontSize', strategy2d.DefaultAxisLabelFontSize);
set(yl, 'FontSize', strategy2d.DefaultAxisLabelFontSize);
colormap(ax, strategy2d.ColorMap);
colorbar(ax);
fileName = sprintf('%s_density_%s.pdf', config.ckptPrefix, config.ckptPostfix);
fileName = fullfile(config.ckptDir, fileName);
exportgraphics(ax, fileName, 'ContentType', 'vector');
end
