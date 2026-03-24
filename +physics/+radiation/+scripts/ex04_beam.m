%==========================================================================
% FileName: ex04_beam
% Author: Yi CAI
% Description:
%   Beam problem for radiation transport. A directional Gaussian source
%   at the origin emits particles in the +x1 direction through a pure
%   scattering medium with periodic boundary conditions.
%==========================================================================

clc, clear, close all;

%% EXECUTION

% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = physics.scheme.ConfigParser.parseFile(common);

% Problem setup
config.id = 'ex04_2dt';

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

% Physical parameters
config.epsilon = 1;
config.decomposition = '';
config.timeScale = 1;
config.scatteringScale = -1;
config.absorptionScale = 1;

% Scattering and absorption
config.scattering = fScattering(config);
config.absorption = [];
config.source = fSource(config);

% Initial and boundary conditions
config.ic = fInit(config);
config.bc = [];
config.exact = [];

% Temporal domain
config.tFinal = 1.0;

% Checkpoint information
config.ckptDir = fullfile(fileparts(mfilename('fullpath')), config.id);
if ~isfolder(config.ckptDir), mkdir(config.ckptDir); end
config.ckptTimeStamps = 0.1:0.1:config.tFinal;
config.verbose = 2;

% Scheme control
config.useFilter = false;
config.filterType = 'exp';
config.filterStrength = 0.01;
config.filterOrder = 4;
config.filterCutOff = 2/3;
config.filterCollisionRate = 1;
config.usePositivityLimiter = true;
config.positivityLimiterType = 'zhang_shu';

% Simulation control
enableRef = 1;
overrideRef = 0;
enableNum = 1;
overrideNum = 0;

% Plotting control
plt = struct();
plt.legend = {};
switch config.nDims
    case 2
        plt.strategy = physics.visual.Strategy2d();
        plt.figIdx = 2;
        plt.colorIdx = 1;
        plt.markerIdx = 1;
end
plt.time = 1.0;
plt.slicePostfix = sprintf('epsilon%.0e', config.epsilon);

% Reference
if enableRef
    config.tOdeIntName = 'ssprk2';
    config.xBasisOrder = 2;
    config.cfl = [];
    config.dt = 0.0005;
    config.nx = repmat(50, 1, config.nDims);
    config.nu = 0;
    config.nv = 100;
    config.ckptPostfix = sprintf('epsilon%.0e', config.epsilon);
    config.ckptPrefix = 'REF';
    switch config.id
        case 'ex04_2dt'
            config.ckptPrefix = strjoin({'BM2DT', config.ckptPrefix}, '_');
    end
    if ~isempty(plt.time) && ~overrideRef
        timeStr = sprintf('t%.1f', plt.time);
        timeStr = strrep(timeStr, '.', 'p');
        fileName = sprintf('%s_state_%s_%s.mat', ...
            config.ckptPrefix, config.ckptPostfix, timeStr);
    else
        fileName = sprintf('%s_state_%s.mat', ...
            config.ckptPrefix, config.ckptPostfix);
    end
    fileName = fullfile(config.ckptDir, fileName);

    if exist(fileName, 'file') && ~overrideRef
        state = physics.radiation.MacroMicroState.load(fileName);
    else
        [config, scheme, state] = run(config);
        state.save(fileName);
    end
    plt.style = 'line';
    plt.legend{end+1} = 'REF';
    switch state.NXDims
        case 2
            plotDensitySlice1d(config, state, plt);
            plotDensity2d(config, state, plt);
            plt.figIdx = plt.figIdx + 1;
            plt.colorIdx = plt.colorIdx + 1;
    end
end

%< Numeric
if enableNum
    config.cfl = [];
    config.dt = 0.05;
    config.nx = repmat(50, 1, config.nDims);
    order = [2];
    % nuv = [0, 19;6,8;10, 0]; % same DoFs
    % nuv = [0, 33;2,30;8,18;17, 0;8, 0;0, 18]; % same DoFs
    % nuv = [6, 0; 6, 4;6, 16];
    % nuv = [1, 8;8, 8;16, 8];
    nuv = [12, 0;12, 4;12, 16];
    for k = 1:size(nuv, 1)
        config.nu = nuv(k, 1);
        config.nv = nuv(k, 2:end);
        config.ckptPostfix = sprintf('epsilon%.0e_nu%d_nv%d', ...
            config.epsilon, config.nu, config.nv);
        for j = 1:length(order)
            switch order(j)
                case 1
                    config.tOdeIntName = 'be';
                    config.xBasisOrder = 1;
                case 2
                    config.tOdeIntName = 'sdirk3';
                    config.xBasisOrder = 2;
                case 3
                    config.tOdeIntName = 'sdirk3';
                    config.xBasisOrder = 3;
            end
            config.ckptPrefix = [upper(config.schemeName), num2str(config.xBasisOrder)];
            switch config.id
                case 'ex04_2dt'
                    config.ckptPrefix = strjoin({'BM2DT', config.ckptPrefix}, '_');
            end

            if ~isempty(plt.time) && ~overrideNum
                timeStr = sprintf('t%.1f', plt.time);
                timeStr = strrep(timeStr, '.', 'p');
                fileName = sprintf('%s_state_%s_%s.mat', ...
                    config.ckptPrefix, config.ckptPostfix, timeStr);
            else
                fileName = sprintf('%s_state_%s.mat', ...
                    config.ckptPrefix, config.ckptPostfix);
            end
            fileName = fullfile(config.ckptDir, fileName);

            if exist(fileName, 'file') && ~overrideNum
                state = physics.radiation.MacroMicroState.load(fileName);
            else
                [config, scheme, state] = run(config);
                state.save(fileName);
            end

            plt.style = 'scatter';
            plt.legend{end+1} = sprintf('MMDG%d', order(j));
            switch state.NXDims
                case 2
                    plotDensitySlice1d(config, state, plt);
                    plotDensity2d(config, state, plt);
            end
            plt.figIdx = plt.figIdx + 1;
            plt.colorIdx = plt.colorIdx + 1;
            plt.markerIdx = plt.markerIdx + 1;
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
visualizer.setComponents(struct('U', [], 'G', []));
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
f = zeros(1, size(x, 2));
end

function h = fSource(config)
h = @(x, v, t) fSourceImpl(x, v, t, config);
end

function f = fSourceImpl(x, v, t, config)
v0 = [sqrt(3)/2; 1/2];
kappa = 100;
sigma = sqrt(1e-2);
nd = size(x, 1);
x0 = zeros(nd, 1);
Gx = exp(-sum((x - x0).^2, 1) / (2 * sigma^2)) / (sqrt(2*pi)*sigma)^nd;
Gv = exp(kappa * v0.' * v) / (2*pi*besseli(0, kappa));
f = Gx .* Gv;
end

function h = fScattering(config)
h = @(x) fScatteringImpl(x, config);
end

function C = fScatteringImpl(x, config)
C = double(x(1, :) > 0.5) * 10;
end


function plotDensity2d(config, state, plt)

figure(plt.figIdx);
ax = gca;
density = state.density(zeros(state.NXDims, 1));
x = state.XDisc.Mesh.collocate(repmat({0}, state.NXDims, 1));
density = reshape(density, length(x{1}), length(x{2}));

imagesc(ax, x{1}, x{2}, density.', 'Interpolation', 'bilinear');
axis(ax, 'xy', 'equal', 'tight');
xl = xlabel(ax, 'x');
yl = ylabel(ax, 'y');
set(ax, 'FontSize', plt.strategy.DefaultFontSize);
set(xl, 'FontSize', plt.strategy.DefaultAxisLabelFontSize);
set(yl, 'FontSize', plt.strategy.DefaultAxisLabelFontSize);
colormap(ax, plt.strategy.ColorMap);
colorbar(ax);
clim([-0.05, 5.05]);
xline(ax, 0.5, '--w', 'LineWidth', 3, 'HandleVisibility', 'off');
% if ~isempty(plt.time)
%     t = plt.time;
% else
%     t = config.tFinal;
% end
% ti = title(sprintf('t = %.2f', t));
% set(ti, 'FontSize', plt.strategy.DefaultTitleFontSize);
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

function plotDensitySlice1d(config, state, plt)

figure(1);
ax = gca;
hold(ax, 'on');
density = state.density(zeros(state.NXDims, 1));
x = state.XDisc.Mesh.collocate(repmat({0}, state.NXDims, 1));
density = reshape(density, length(x{1}), length(x{2}));

r = linspace(min(x{1}), max(x{1}), 256);
s = interp2(x{1}, x{2}, density.', r*sqrt(3)/2, r/2);

strategy = physics.visual.Strategy1d();
if strcmpi(plt.style, 'line')
    style = strategy.getDefaultLineStyle(plt.colorIdx, plt.markerIdx);
else
    style = strategy.getDefaultScatterStyle(plt.colorIdx, plt.markerIdx);
end
if contains(plt.legend{end}, 'REF')
    label = sprintf('%s', plt.legend{end});
else
    if config.nu == 0
        label = sprintf('S_{%d}', config.nv);
    elseif config.nv == 0
        label = sprintf('P_{%d}', config.nu-1);
    else
        label = sprintf('P_{%d}S_{%d}', config.nu-1, config.nv);
    end
end
plot(ax, r, s, style{:}, 'DisplayName', label);
xline(ax, 0.5, '--r', 'LineWidth', 3, 'HandleVisibility', 'off');
legend(ax, 'show', 'FontSize', strategy.DefaultLegendFontSize, ...
    'Location', strategy.DefaultLegendPosition);
xlim([-config.L, config.L]);
xticks(-config.L:0.25:config.L);
ylim([-0.05, 5.05]);
xl = xlabel(ax, 'x');
yl = ylabel(ax, '\rho');
set(ax, 'FontSize', strategy.DefaultFontSize);
set(xl, 'FontSize', strategy.DefaultAxisLabelFontSize);
set(yl, 'FontSize', strategy.DefaultAxisLabelFontSize);
if ~isempty(plt.time)
    t = plt.time;
else
    t = config.tFinal;
end
ti = title(sprintf('t = %.2f', t));
set(ti, 'FontSize', plt.strategy.DefaultTitleFontSize);
if ~isempty(plt.time)
    timeStr = sprintf('t%.1f', plt.time);
    timeStr = strrep(timeStr, '.', 'p');
    fileName = sprintf('%s_slice_%s_%s.pdf', config.ckptPrefix, plt.slicePostfix, timeStr);
else
    fileName = sprintf('%s_slice_%s.pdf', config.ckptPrefix, plt.slicePostfix);
end
fileName = fullfile(config.ckptDir, fileName);
exportgraphics(ax, fileName, 'ContentType', 'vector');
end
