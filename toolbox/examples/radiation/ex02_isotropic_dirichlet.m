function [state, result] = ex02_isotropic_dirichlet(options)
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
% FileName: ex02_isotropic_dirichlet
% Author: Yi CAI
% Description:
%   Isotropic Dirichlet boundary condition test for radiation transport.
%==========================================================================


%% EXECUTION

% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = mdemat.tools.environment.readConfig(common);

% Problem setup
config.id = 'ex02_2dt';

% Parse settings from problem ID
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

% Spatial domain [0, L]^d (not symmetric)
config.L = 1;
config.xBBox = repmat([0, config.L], 1, config.nDims);
config.nx = repmat(64, 1, config.nDims);

% Physical parameters
config.decomposition = '';
config.timeScale = 1;
config.scatteringScale = -1;
config.absorptionScale = 1;

% Velocity discretization parameters based on dimensions
switch config.nDims
    case 1
        if strcmp(config.vDimReduction, 'topology')
            config.nu = 1; % order: 0 or 1
            config.nv = 2; % micro coupling
        else
            config.nu = 2; % order: 0 for none; nu - 1 for degree
            config.nv = 16; % micro coupling
        end
    case 2
        if strcmp(config.vDimReduction, 'topology')
            config.nu = 2; % order: 0 for none; nu - 1 for degree
            config.nv = 16; % micro coupling
        else
            config.nu = 2; % order: 0 for none; nu - 1 for degree
            config.nv = [4, 4]; % micro coupling
        end
end

% Checkpoint information
config.ckptDir = fullfile(fileparts(mfilename('fullpath')), config.id);
if options.Write && ~isfolder(config.ckptDir), mkdir(config.ckptDir); end
config.ckptTimeStamps = [0.1, 0.4, 1.0, 1.6, 4.0];
config.verbose = 2;
config.tFinal = 4.0;
config.cfl = [];
config.dt = 0.1;

if options.Smoke
    config.epsilon = 1;
    config.tau0 = 1;
    [config, scheme, state, result] = run(config, options); %#ok<ASGLU>
    return;
end

% Parameter sweep
override = 0;
order = [1];
epsilon = 10.^0;
for k = 1:length(epsilon)
    plt = struct();
    plt.legend = {};
    if override 
        plt.time = config.tFinal;
    else
        plt.time = 4.0;
    end
    switch config.nDims
        case 1
            plt.styleOptions = mdemat.tools.visual.plotStyle();
        case 2
            plt.styleOptions = mdemat.tools.visual.plotStyle();
    end

    config.epsilon = epsilon(k);
    config.ckptPostfix = sprintf('epsilon%.0e', epsilon(k));

    %< Reference 
    config.ckptPrefix = [upper(config.schemeName), num2str(1)];
    switch config.id
        case 'ex02_1ds'
            config.ckptPrefix = strjoin({'ID1DS', config.ckptPrefix}, '_');
        case 'ex02_2dt'
            config.ckptPrefix = strjoin({'ID2DT', config.ckptPrefix}, '_');
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
        loadOptions = options;
        loadOptions.LoadFile = fileName;
        [config, ~, state, result] = run(config, loadOptions); 
        plt.style = 'line';
        plt.legend{end+1} = 'REF';
        plt.colorIdx = 1;
        plt.markerIdx = 1;
        switch state.NXDims
            case 1
                plt.figIdx = 2;
                plotDensity1d(config, state, plt);
            case 2
                plt.figIdx = 1;
                plotDensity2d(config, state, plt);
        end
    end

    %< Numeric
    for j = 1:length(order)
        switch config.nDims
            case 1
                plt.figIdx = 2;
            case 2
                plt.figIdx = j + 1;
        end
        switch order(j)
            case 1
                config.tOdeIntName = 'fe';
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
            case 'ex02_1ds'
                config.ckptPrefix = strjoin({'ID1DS', config.ckptPrefix}, '_');
            case 'ex02_2dt'
                config.ckptPrefix = strjoin({'ID2DT', config.ckptPrefix}, '_');
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
            loadOptions = options;
        loadOptions.LoadFile = fileName;
        [config, ~, state, result] = run(config, loadOptions);
        else
            [config, ~, state, result] = run(config, options);
            state.save(fileName, omitCallbacks=true, createDirectory=true, overwrite=logical(override));
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


if isempty(result) && ~isempty(state)
    result = struct('schemaVersion', 1, 'completed', true, 'diagnostics', state.History);
end

end

function [config, scheme, state, result] = run(config, options)

% Scattering and absorption
config.scattering = fScattering(config);
config.absorption = [];
config.source = [];

% Initial and boundary conditions
config.ic = fInit(config);
config.bc = fBc(config);
config.exact = [];

config = mdemat.internal.execution.studyConfig(config, options);

% Construct and configure scheme
switch config.schemeName
    case 'mmdg'
        scheme = mdemat.application.radiation.MmDgScheme(config=config);
end

% Configure visualizer
plotOptions = scheme.Observation.Plot;
plotOptions.PointsPerDimension = 1;
plotOptions.Times = linspace(0, scheme.getConfig('tFinal'), 10+1);
plotOptions.Components = struct('U', [], 'G', [], 'F', 1);
plotOptions.DatasetTypes.(scheme.getConfig('sId')) = 'numeric';
switch config.nDims
    case 1
        if options.Plot, plotOptions.Modes{end+1} = '1d'; end
    case 2
        if options.Plot, plotOptions.Modes{end+1} = '2d'; end
        if options.Plot, plotOptions.Modes{end+1} = 'slice1d'; end
end
scheme.setObservation('plot', plotOptions);

% Build velocity discretization (SumSpace)
if config.nu > 0
    vMacroElement = mdemat.common.fem.L2SphereElement.modal(config.nDims, ...
        config.nu, reduction=config.vDimReduction, np=config.nv);
    vMacroDisc = mdemat.common.fem.SpectralSpace(vMacroElement);
else
    vMacroDisc = [];
end

if all(config.nv > 0)
    if config.nDims == 2 && strcmp(config.vDimReduction, 'symmetry')
        nVDims = 3;
    else
        nVDims = config.nDims;
    end
    vMicroElement = mdemat.common.fem.L2SphereElement.nodal(nVDims, ...
        config.nv, reduction=config.vDimReduction);
    vMicroDisc = mdemat.common.fem.SpectralSpace(vMicroElement);
else
    vMicroDisc = [];
end

vDisc = mdemat.common.fem.SumSpace(vMacroDisc, vMicroDisc);

% Build spatial discretization
switch scheme.getConfig('xBasisType')
    case 'modal'
        xElement = mdemat.common.fem.BH1OrthotopeElement.modal(config.nDims, ...
            scheme.getConfig('xBasisOrder'), ...
            pattern=scheme.getConfig('xBasisPattern'));
    case 'nodal'
        xElement = mdemat.common.fem.BH1OrthotopeElement.nodal(config.nDims, ...
            scheme.getConfig('xBasisOrder'));
end
xMesh = mdemat.data.mesh.UniformGrid(config.xBBox, config.nx);
xDisc = mdemat.common.fem.FiniteElementSpace(xMesh, xElement);

% Create state and run simulation
state = mdemat.application.radiation.MacroMicroState(xDisc, vDisc);
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
E = config.E;
L = config.L;
fL = 1; % Value at x = 0
fR = 0; % Value at x = L
f = (any(x == 0, 1) .* fL + any(x == L, 1) .* fR) .* E;
end

function h = fScattering(config)
h = @(x) fScatteringImpl(x, config);
end

function C = fScatteringImpl(x, config)
C = ones(1, size(x, 2));
end


function plotDensity1d(config, state, plt)
figure(plt.figIdx);
ax = gca;

x = state.XDisc.Mesh.collocate(0);
if config.nu > 0
    density = state.XDisc.eval(0, state.Dofs.U(:, 1));
else
    w = state.VDisc.Rhs.Element.Volume.Weights;
    v = state.VDisc.Rhs.Element.Volume.Nodes;
    F = state.XDisc.eval(0, state.Dofs.F);
    F = state.VDisc.Rhs.eval(v, F.');
    density = F * w(:);
end
hold(ax, 'on');
switch plt.style
    case 'line'
        style = namedargs2cell(mdemat.tools.visual.curveStyle(plt.styleOptions, ...
            plt.colorIdx, plt.markerIdx, kind="line"));
    case 'scatter'
        style = namedargs2cell(mdemat.tools.visual.curveStyle(plt.styleOptions, ...
            plt.colorIdx, plt.markerIdx, kind="scatter"));
end
line = plot(x, density, style{:});
set(line, 'LineWidth', 3, 'MarkerSize', 8);
leg = legend(ax, plt.legend);
set(leg, 'Location', plt.styleOptions.LegendLocation);
set(leg, 'FontSize', plt.styleOptions.LegendFontSize);
set(leg, 'Interpreter', 'latex');
hold(ax, 'off');
xlim([0, 1]);
xl = xlabel('x');
set(ax, 'FontSize', plt.styleOptions.FontSize);
set(xl, 'FontSize', plt.styleOptions.AxisLabelFontSize);
ti = title('Density');
set(ti, 'FontSize', plt.styleOptions.TitleFontSize);
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

figure(plt.figIdx);
ax = gca;
density = state.density(zeros(state.NXDims, 1));
x = state.XDisc.Mesh.collocate(repmat({0}, state.NXDims, 1));
density = reshape(density, length(x{1}), length(x{2}));

imagesc(ax, x{1}, x{2}, density.', 'Interpolation', 'bilinear');
axis(ax, 'xy', 'equal', 'tight');
xl = xlabel(ax, 'x');
yl = ylabel(ax, 'y');
set(ax, 'FontSize', plt.styleOptions.FontSize);
set(xl, 'FontSize', plt.styleOptions.AxisLabelFontSize);
set(yl, 'FontSize', plt.styleOptions.AxisLabelFontSize);
colormap(ax, plt.styleOptions.ColorMap);
colorbar(ax);
ti = title(sprintf('t = %.2f', plt.time));
set(ti, 'FontSize', plt.styleOptions.TitleFontSize);
if ~isempty(plt.time)
    timeStr = sprintf('t%.1f', plt.time);
    timeStr = strrep(timeStr, '.', 'p');
    fileName = sprintf('%s_density_%s_%s.pdf', config.ckptPrefix, config.ckptPostfix, timeStr);
else
    fileName = sprintf('%s_density_%s.pdf', config.ckptPrefix, config.ckptPostfix);
end
caxis([-0.05, 0.85]);
fileName = fullfile(config.ckptDir, fileName);
exportgraphics(ax, fileName, 'ContentType', 'vector');
end
