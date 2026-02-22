%==========================================================================
% FileName: ex04_lattice
% Author: Yi CAI
% Description:
%   Lattice (checkerboard) problem for radiation transport. A 7x7 grid of
%   alternating scattering and absorbing material cells with a fixed
%   isotropic source in the center cell.
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

% Spatial domain [0, 7]^d
config.L = 7;
config.xBBox = repmat([0, config.L], 1, config.nDims);
config.nx = repmat(35, 1, config.nDims);

% Physical parameters
config.decomposition = '';
config.timeScale = 1;
config.scatteringScale = -1;
config.absorptionScale = 1;

% Scattering, absorption, and source
config.scattering = fScattering(config);
config.absorption = fAbsorption(config);
config.source = fSource(config);

% Initial and boundary conditions
config.ic = fInit(config);
config.bc = fBc(config);
config.exact = [];

% Velocity discretization parameters based on dimensions
switch config.nDims
    case 1
        if strcmp(config.vDimReduction, 'topology')
            config.nu = 1;
            config.nv = 2;
            config.nup = config.nu;
        else
            config.nu = 2;
            config.nv = 16;
            config.nup = config.nu + 2;
        end
    case 2
        if strcmp(config.vDimReduction, 'topology')
            config.nu = 2;
            config.nv = 16;
            config.nup = 2*config.nu + 2;
        else
            config.nu = 2;
            config.nv = 16;
            config.nup = [config.nu + 1, 2*(config.nu-1) + 2];
        end
end

% Checkpoint information
config.ckptDir = fullfile(fileparts(mfilename('fullpath')), config.id);
if ~isfolder(config.ckptDir), mkdir(config.ckptDir); end
config.ckptTimeStamps = [0.5, 1.0, 2.0, 3.0, 3.5];
config.verbose = 2;
config.tFinal = 3.5;

% Parameter sweep
override = 0;
order = [2];
epsilon = 10.^[0];
for k = 1:length(epsilon)
    plt = struct();
    plt.legend = {};
    plt.time = 3.5;
    config.epsilon = epsilon(k);
    config.ckptPostfix = sprintf('epsilon%.0e', epsilon(k));

    %< Reference
    config.ckptPrefix = [upper(config.schemeName), num2str(1)];
    switch config.id
        case 'ex04_2dt'
            config.ckptPrefix = strjoin({'LT2DT', config.ckptPrefix}, '_');
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
        plotDensity2d(config, state, plt);
    end

    for j = 1:length(order)
        %< Numeric
        switch order(j)
            case 1
                config.tOdeIntName = 'be';
                config.xBasisOrder = 1;
            case 2
                config.tOdeIntName = 'sdirk2';
                config.xBasisOrder = 2;
            case 3
                config.tOdeIntName = 'sdirk3';
                config.xBasisOrder = 3;
        end
        config.ckptPrefix = [upper(config.schemeName), num2str(config.xBasisOrder)];
        switch config.id
            case 'ex04_2dt'
                config.ckptPrefix = strjoin({'LT2DT', config.ckptPrefix}, '_');
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
            plotDensity2d(config, state, plt);
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
ci = min(6, max(0, floor(x(1, :))));
cj = min(6, max(0, floor(x(2, :))));
isScattering = (mod(ci + cj, 2) == 0);
C = isScattering * 1.0;
end

function h = fAbsorption(config)
h = @(x) fAbsorptionImpl(x, config);
end

function C = fAbsorptionImpl(x, config)
ci = min(6, max(0, floor(x(1, :))));
cj = min(6, max(0, floor(x(2, :))));
isAbsorbing = (mod(ci + cj, 2) == 1);
C = isAbsorbing * 10.0;
end

function h = fSource(config)
h = @(x, v, t) fSourceImpl(x, v, t, config);
end

function S = fSourceImpl(x, v, t, config)
inCenter = (x(1,:) >= 3 & x(1,:) <= 4 & x(2,:) >= 3 & x(2,:) <= 4);
S = double(inCenter);
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
