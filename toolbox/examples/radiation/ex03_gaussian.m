function [state, result] = ex03_gaussian(options)
%   Smoke=true executes this original model at bounded resolution and time,
%   returns actual results, and skips figures/files. Defaults retain the full
%   research parameters and sweeps below.
% Run the retained research study; return its final workspace state and assessment.
arguments
    options.Smoke (1, 1) logical = false
end
options.Plot = ~options.Smoke;
options.Write = ~options.Smoke;


%==========================================================================
% FileName: ex03_gaussian
% Author: Yi CAI
% Description:
%   Gaussian pulse test for radiation transport. An isotropic Gaussian
%   density pulse propagates outward through a pure scattering medium.
%==========================================================================


%% EXECUTION

% Common configuration file
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config = mdemat.tools.environment.readConfig(common);

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

% Physical parameters
config.epsilon = 1;
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
config.bc = [];
config.exact = [];

% Temporal domain 
config.tFinal = 0.6;

% Checkpoint information
config.ckptDir = fullfile(fileparts(mfilename('fullpath')), config.id);
if options.Write && ~isfolder(config.ckptDir), mkdir(config.ckptDir); end
config.ckptTimeStamps = 0:0.05:config.tFinal;
config.verbose = 2;

% Scheme control
config.useFilter = false;
config.filterType = 'rot';
config.filterStrength = 0.01;
config.filterOrder = 4;
config.filterCutoff = 2/3;
config.filterCollisionRate = 1;
config.usePositivityLimiter = true;
config.positivityLimiterType = 'zhang_shu';

% Simulation control 
enableRef = 1;
overrideRef = 0;
enableNum = 1;
overrideNum = 0;

if options.Smoke
    config.epsilon = 1;
    config.tau0 = 1;
    [config, scheme, state, result] = run(config, options); %#ok<ASGLU>
    return;
end

% Plotting control
plt = struct();
plt.angleFigIdx = 1;
plt.sliceFigIdx = 2;
plt.legend = {};
switch config.nDims
    case 2
        plt.styleOptions = mdemat.tools.visual.plotStyle();
        plt.figIdx = plt.sliceFigIdx + 1;
        plt.colorIdx = 1;
        plt.markerIdx = 1;
end
plt.time = 0.6;
plt.slicePostfix = sprintf('epsilon%.0e', config.epsilon);

% Reference
if enableRef
    config.tOdeIntName = 'ssprk2';
    config.xBasisOrder = 2;
    config.cfl = [];
    config.dt = 0.0005;
    config.nx = repmat(50, 1, config.nDims);
    config.nu = 1;
    config.nv = 100;
    config.ckptPostfix = sprintf('epsilon%.0e', config.epsilon);
    config.ckptPrefix = 'REF';
    switch config.id
        case 'ex03_2dt'
            config.ckptPrefix = strjoin({'GS2DT', config.ckptPrefix}, '_');
        case 'ex03_2ds'
            config.ckptPrefix = strjoin({'GS2DS', config.ckptPrefix}, '_');
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
        loadOptions = options;
        loadOptions.LoadFile = fileName;
        [config, ~, state, result] = run(config, loadOptions);
    else
        [config, ~, state, result] = run(config, options);
        state.save(fileName, omitCallbacks=true, createDirectory=true, overwrite=logical(overrideRef));
    end
    plt.style = 'line';
    plt.legend{end+1} = 'REF';
    switch state.NXDims
        case 2
            plotAngleDist2d(config, state, plt);
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
    % nuv = [4,0;6,0;8,0;10,0;12,0;14,0;16,0;18,0;20,0]; % varying
    % nuv = [0,7;0,11;0,15;0,19;0,23;0,27;0,31;0,35;0,39]; % varying
    % nuv = [2,36;8,24;16,8;20,0];
    % nuv = [0, 19;1, 18;2, 16;6, 8;8, 4;10, 0]; % same DoFs
    % nuv = [0, 33;1, 32; 2, 30; 9, 16; 12, 8; 17, 0];
    % Select "varyNu" to run the alternate fixed-nv sweep.
    nuv = gaussianSweep("varyNv");
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
                case 'ex03_2dt'
                    config.ckptPrefix = strjoin({'GS2DT', config.ckptPrefix}, '_');
                case 'ex03_2ds'
                    config.ckptPrefix = strjoin({'GS2DS', config.ckptPrefix}, '_');
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
                loadOptions = options;
        loadOptions.LoadFile = fileName;
        [config, ~, state, result] = run(config, loadOptions);
            else
                [config, ~, state, result] = run(config, options);
                state.save(fileName, omitCallbacks=true, createDirectory=true, overwrite=logical(overrideNum));
            end
    
            plt.legend{end+1} = sprintf('MMDG%d', order(j));
            switch state.NXDims
                case 2
                    plt.style = 'line';
                    plotAngleDist2d(config, state, plt);
                    plt.style = 'scatter';
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


if isempty(result) && ~isempty(state)
    result = struct('schemaVersion', 1, 'completed', true, 'diagnostics', state.History);
end

end

function [config, scheme, state, result] = run(config, options)

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
        config.nu, reduction=config.vDimReduction);
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
E = config.E;
sigma = sqrt(1e-2);
nd = size(x, 1);
f = exp(-sum(x.^2, 1) / (2 * sigma^2)) / (sqrt(2*pi)*sigma)^nd * E;
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


function plotDensity2d(config, state, plt)

figure(plt.figIdx);
ax = gca;
xRefCell = repmat({0}, state.NXDims, 1);
xRef = cell(size(xRefCell));
[xRef{:}] = ndgrid(xRefCell{:});
xRef = reshape(cat(ndims(xRef{1})+1, xRef{:}), [], numel(xRefCell)).';
density = state.density(xRef);
x = state.XDisc.Mesh.collocate(xRefCell);
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

if plt.time == 0.3
    clim([-0.05, 8.05]);
elseif plt.time == 0.4
    clim([-0.05, 5.55]);
elseif plt.time == 0.6
    clim([-0.05, 0.85]);
end

% if ~isempty(plt.time)
%     t = plt.time;
% else
%     t = config.tFinal;
% end
% ti = title(sprintf('t = %.2f', t));
% set(ti, 'FontSize', plt.styleOptions.TitleFontSize);
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

figure(plt.sliceFigIdx);
ax = gca;
hold(ax, 'on');
density = state.density(zeros(state.NXDims, 1));
x = state.XDisc.Mesh.collocate(repmat({0}, state.NXDims, 1));
density = reshape(density, length(x{1}), length(x{2}));

R = sqrt(state.NXDims);
r = linspace(min(x{1}), max(x{1}), 256);
s = interp2(x{1}, x{2}, density.', r/R, r/R);

styleOptions = mdemat.tools.visual.plotStyle();
if strcmpi(plt.style, 'line')
    style = namedargs2cell(mdemat.tools.visual.curveStyle(styleOptions, ...
        plt.colorIdx, plt.markerIdx, kind="line"));
else
    style = namedargs2cell(mdemat.tools.visual.curveStyle(styleOptions, ...
        plt.colorIdx, plt.markerIdx, kind="scatter"));
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
legend(ax, 'show', 'FontSize', styleOptions.LegendFontSize, ...
    'Location', styleOptions.LegendLocation);
xlim([0, config.L]);
xticks(0:0.2:config.L);

if plt.time == 0.3
    ylim([-0.05, 8.05]);
elseif plt.time == 0.4
    ylim([-0.05, 5.55]);
elseif plt.time == 0.6
    ylim([-0.05, 1.05]);
end

xl = xlabel(ax, 'x');
yl = ylabel(ax, '\rho');
set(ax, 'FontSize', styleOptions.FontSize);
set(xl, 'FontSize', styleOptions.AxisLabelFontSize);
set(yl, 'FontSize', styleOptions.AxisLabelFontSize);
if ~isempty(plt.time)
    t = plt.time;
else
    t = config.tFinal;
end
ti = title(ax, sprintf('t = %.2f', t));
set(ti, 'FontSize', styleOptions.TitleFontSize);

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

function plotAngleDist2d(config, state, plt)

ne = state.XDisc.NMeshElements;
nvq = 512;
theta = linspace(0, 2*pi, nvq+1);
nxq = state.XDisc.Element.Volume.NPoints;
wxq = state.XDisc.Element.Volume.Weights;
xRef = state.XDisc.Element.Volume.Nodes;
detJac = state.XDisc.Mesh.computeElementJacobianDeterminants();

M = 1/(2*pi);
F = state.distribution(xRef, theta);
E = wxq * reshape((F-M).^2/M/2, nxq, []);
E = detJac(:).' * reshape(E, ne, []);

figure(plt.angleFigIdx);
ax = gca;
hold(ax, 'on');
styleOptions = mdemat.tools.visual.plotStyle();
if strcmpi(plt.style, 'line')
    style = namedargs2cell(mdemat.tools.visual.curveStyle(styleOptions, ...
        plt.colorIdx, plt.markerIdx, kind="line"));
else
    style = namedargs2cell(mdemat.tools.visual.curveStyle(styleOptions, ...
        plt.colorIdx, plt.markerIdx, kind="scatter"));
end
if contains(plt.legend{end}, 'REF')
    label = plt.legend{end};
else
    if config.nu == 0
        label = sprintf('S_{%d}', config.nv);
    elseif config.nv == 0
        label = sprintf('P_{%d}', config.nu-1);
    else
        label = sprintf('P_{%d}S_{%d}', config.nu-1, config.nv);
    end
end
plot(ax, theta, E, style{:}, 'DisplayName', label);
legend(ax, 'show', 'FontSize', styleOptions.LegendFontSize, ...
    'Location', styleOptions.LegendLocation);
xlim([0, 2*pi]);
xticks(0:pi/4:2*pi);
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi','5\pi/4','3\pi/2','7\pi/4','2\pi'});
xl = xlabel(ax, '\theta');
% yl = ylabel(ax, '\rho(\theta)');
set(ax, 'FontSize', styleOptions.FontSize);
set(xl, 'FontSize', styleOptions.AxisLabelFontSize);
% set(yl, 'FontSize', styleOptions.AxisLabelFontSize);

if ~isempty(plt.time)
    timeStr = sprintf('t%.1f', plt.time);
    timeStr = strrep(timeStr, '.', 'p');
    fileName = sprintf('%s_angle_%s_%s.pdf', config.ckptPrefix, plt.slicePostfix, timeStr);
else
    fileName = sprintf('%s_angle_%s.pdf', config.ckptPrefix, plt.slicePostfix);
end
fileName = fullfile(config.ckptDir, fileName);
exportgraphics(ax, fileName, 'ContentType', 'vector');
end
