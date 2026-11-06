function verifyInstalledPackage()
%VERIFYINSTALLEDPACKAGE Standalone child-process check; needs only the .mltbx.
work = pwd;
resetSearchPath();
fprintf('INSTALL_CHECK native folder: %s\n', work);
report = struct('passed', false, 'release', version('-release'), ...
    'computer', computer, 'workingDirectory', work, 'checks', {{}});
setting = settings;
folderSetting = setting.matlab.addons.InstallationFolder;
hadTemporary = hasTemporaryValue(folderSetting);
if hadTemporary, previous = folderSetting.TemporaryValue; else, previous = ''; end
installed = [];
temporaryFolderSelected = false;
settingsGuard = onCleanup(@() restoreSetting(folderSetting, hadTemporary, previous));
uninstallGuard = [];
try
    % Query the caller's installation before switching to the disposable root.
    prior = matlab.addons.toolbox.installedToolboxes;
    if ~isempty(prior) && any(strcmp({prior.Guid}, '41e45d88-ea42-4d8f-bcac-46759e59ed38'))
        error('mdemat:build:ExistingInstallation', 'An existing MdeMat installation is registered.');
    end
    folderSetting.TemporaryValue = fullfile(work, 'installed-addons');
    temporaryFolderSelected = true;
    installed = matlab.addons.toolbox.installToolbox(fullfile(work, 'MdeMat.mltbx'), true);
    % Capture the returned registration by value. No nested mutable workspace
    % is referenced when the guard is destroyed during success or unwinding.
    uninstallGuard = onCleanup(@() uninstallOwned(installed));
    fprintf('INSTALL_CHECK installed: %s\n', installed.Name);
    source = which('mdemat.common.fdm.secondDerivative');
    assert(contains(source, fullfile(work, 'installed-addons')), ...
        'mdemat:build:WrongInstalledPath', 'Code did not resolve from temporary installation.');
    assert(isempty(which('approx.integrate.LebedevRule')) && ...
        isempty(which('core.function.CompiledFunction')), ...
        'mdemat:build:OldNamespaceVisible', 'An old source namespace remains visible.');
    basis = mdemat.common.approximation.CompiledFunction(1, 4, 'monic_unit_legendre');
    x = [-0.3, 0.1, 0.4];
    assert(max(abs(basis.eval(x) - [ones(size(x)); x; x.^2-1/12; x.^3-3*x/20]), ...
        [], 'all') < 2e-14, 'mdemat:build:ResourceInvariant', 'Installed basis resource failed.');
    sphere = mdemat.data.geometry.Sphere([0, 0, 0], 1);
    rule = mdemat.common.quadrature.LebedevRule();
    [nodes, weights] = rule.generate(sphere, 14);
    assert(abs(sum(weights)-4*pi) < 1e-13 && size(nodes, 2) == 14, ...
        'mdemat:build:VendorInvariant', 'Installed Lebedev quadrature failed.');
    [L, nodes, lift] = mdemat.common.fdm.secondDerivative(7, [0, 1]);
    exact = nodes .* (1-nodes);
    computed = (-L) \ (2*ones(size(nodes)) + lift*[0; 0]);
    assert(norm(computed-exact, Inf) < 1e-13, ...
        'mdemat:build:SolverInvariant', 'Installed finite-difference solve failed.');
    assert(contains(help('mdemat.common.fdm.secondDerivative'), 'Dirichlet'), ...
        'mdemat:build:HelpMissing', 'Native help is unavailable.');
    toolboxRoot = fileparts(fileparts(fileparts(fileparts(source))));
    assert(isfile(fullfile(toolboxRoot, 'doc', 'GettingStarted.m')), ...
        'mdemat:build:HelpMissing', 'Getting Started guide is missing.');
    introduction = runGuide(toolboxRoot);
    guideError = max(abs(introduction.Values-introduction.Reference), [], 'all');
    assert(guideError < 1e-13 && introduction.Assessment.L2 < 1e-13, ...
        'mdemat:build:GuideInvariant', 'The callable guide failed its continuous quadratic solution.');
    space = mdemat.common.fem.FiniteElementSpace( ...
        mdemat.data.mesh.UniformGrid([0, 1], 4), ...
        mdemat.common.fem.BH1OrthotopeElement.modal(1, 3));
    problem = mdemat.solver.elliptic.poissonProblem(@(x) -2*ones(1, size(x, 2)), ...
        boundaryType="dirichlet", boundaryValue=@(x) x(1, :).^2);
    [solution, execution] = mdemat.solver.elliptic.solvePoisson(problem, space);
    referenceNodes = space.Element.Volume.Nodes;
    physicalNodes = space.Mesh.collocate(referenceNodes);
    potentialError = norm(space.eval(referenceNodes, solution.U)-physicalNodes(:).^2, Inf);
    gradientError = norm(space.eval(referenceNodes, solution.Q)-2*physicalNodes(:), Inf);
    assert(execution.Converged && potentialError < 2e-10 && gradientError < 2e-10, ...
        'mdemat:build:PoissonInvariant', 'Installed pure Poisson potential/gradient failed.');
    state = mdemat.common.fem.SpatialState(space);
    state.setDof('U', solution.U); state.setDof('Q', solution.Q);
    stateFile = fullfile(work, 'state.mat');
    state.save(stateFile);
    stateCopy = mdemat.common.fem.SpatialState.load(stateFile);
    assert(isequal(stateCopy.Dofs, state.Dofs) && ...
        norm(stateCopy.XDisc.eval(referenceNodes, stateCopy.Dofs.U)-physicalNodes(:).^2, Inf) < 2e-10, ...
        'mdemat:build:StateInvariant', 'Installed numeric state/space roundtrip failed.');
    families = verifyEvolutionFamilies();
    meshing = verifyDistMesh(toolboxRoot);
    verifyTools(toolboxRoot);
    identity = mdemat.tools.environment.versionInfo();
    manifest = jsondecode(fileread(fullfile(toolboxRoot, 'manifest.json')));
    payload = string(manifest.files);
    normalizedPayload = replace(payload, '\', '/');
    forbidden = regexp(cellstr(normalizedPayload), ...
        '(^|/)(docs|tests|benchmarks|buildUtilities|outputs|release|legacy|\.git)(/|$)', 'once');
    assert(strcmp(manifest.version, identity.ToolboxVersion) && ...
        all(isfile(fullfile(toolboxRoot, payload))) && ...
        all(cellfun(@isempty, forbidden)) && ~any(endsWith(lower(payload), '.pdf')), ...
        'mdemat:build:ManifestInvariant', 'Installed manifest/payload/version check failed.');
    dependencies = string(manifest.dependencyFiles);
    assert(all(startsWith(dependencies, 'toolbox/') | startsWith(dependencies, 'MATLABROOT/')), ...
        'mdemat:build:ManifestInvariant', 'The dependency manifest contains development paths.');
    run = mdemat.tools.io.createRunRecord(struct('U', computed), struct('Status', 'COMPLETED'));
    recordFile = fullfile(work, 'verification.json');
    mdemat.tools.io.writeRecord(recordFile, run);
    [restored, recordInfo] = mdemat.tools.io.readRecord(recordFile);
    assert(strcmp(restored.ValidationReport.Status, 'UNVERIFIED') && ...
        strcmp(recordInfo.ToolboxVersion, identity.ToolboxVersion), ...
        'mdemat:build:RunRecordInvariant', 'Installed run-record metadata failed.');
    try
        mdemat.tools.environment.projectRoot();
        error('mdemat:build:UnexpectedCheckout', 'An installed add-on claimed to be a checkout.');
    catch exception
        assert(strcmp(exception.identifier, 'mdemat:environment:NotCheckout'), ...
            'mdemat:build:UnexpectedCheckout', 'Installed project-root policy failed.');
    end
    report.installedSource = source;
    report.manifest = struct('version', manifest.version, 'files', numel(payload), ...
        'dependencyProducts', {manifest.dependencyProducts});
    report.measurements = struct('GuideMaxError', guideError, ...
        'PoissonPotentialMaxError', potentialError, 'PoissonGradientMaxError', gradientError, ...
        'EvolutionFamilies', families, 'DistMesh', meshing);
    report.checks = {'temporary native install', 'foreign directory/factory path', ...
        'old namespace absent', 'numeric basis resources', 'Lebedev area invariant', ...
        'independent quadratic FDM solution', 'native help/guide', ...
        'callable guide continuous solution', 'pure Poisson potential/gradient', ...
        'numeric State/space roundtrip', ...
        'periodic advection/diffusion constant and mass', ...
        'radiation scattering/absorption/source equilibrium without observers', ...
        'AP and DG Hermite warm-Maxwellian equilibrium without observers', ...
        'low-rank uniform transport and mass', ...
        'bundled DistMesh circle/refinement with scoped path/RNG/figures', ...
        'installed manifest/version/payload/dependency paths', ...
        'versioned RunRecord/IO metadata', 'installed checkout discovery rejected', ...
        'value Config isolation and validation', 'packaged external Graph plot method', ...
        'value field rendering and original component IDs', 'removed tool classes absent', ...
        'consolidated analyze namespace and Timer'};
    clear uninstallGuard;
    assertInstallationRemoved(source);
    clear settingsGuard;
    assertSettingRestored(folderSetting, hadTemporary, previous);
    report.cleanupComplete = true;
    report.passed = true;
catch exception
    report.errorIdentifier = exception.identifier;
    report.reason = getReport(exception, 'extended', 'hyperlinks', 'off');
    try
        if exist('uninstallGuard', 'var'), clear uninstallGuard; end
        if temporaryFolderSelected
            removeOwnedRegistration();
        end
    catch cleanupException
        report.cleanupError = getReport(cleanupException, 'extended', 'hyperlinks', 'off');
    end
    if exist('settingsGuard', 'var'), clear settingsGuard; end
    try
        assertSettingRestored(folderSetting, hadTemporary, previous);
    catch cleanupException
        report.settingsCleanupError = getReport(cleanupException, 'extended', 'hyperlinks', 'off');
    end
end

stream = fopen(fullfile(work, 'installed-package.json'), 'w', 'n', 'UTF-8');
fileGuard = onCleanup(@() fclose(stream)); %#ok<NASGU>
fprintf(stream, '%s\n', jsonencode(report, PrettyPrint=true));
assert(report.passed, 'mdemat:build:InstalledArtifactFailed', 'Installed artifact verification failed.');
end

function verifyTools(toolboxRoot)
timer = mdemat.tools.analyze.Timer(verbose=false);
[~, interval] = timer.start(record='installed');
timer.stop(interval=interval);
snapshot = timer.snapshot();
assert(~snapshot.IsActive && numel(snapshot.Records)==1, ...
    'mdemat:build:TimerInvariant','The consolidated Timer failed its installed lifecycle.');
for package = {'+analysis','+timing','+except'}
    assert(~isfolder(fullfile(toolboxRoot,'+mdemat','+tools',package{1})), ...
        'mdemat:build:OldToolVisible','An obsolete tools package remains installed.');
end
config = mdemat.tools.environment.Config();
config = config.define('order',2,@(x) isnumeric(x)&&isscalar(x)&&x>=1&&x==fix(x));
changed = config.set('order',3);
assert(config.get('order')==2 && changed.get('order')==3, ...
    'mdemat:build:ConfigInvariant','Value Config copies must update independently.');
rejected = false;
try
    changed.set('order',0);
catch
    rejected = true;
end
assert(rejected,'mdemat:build:ConfigInvariant','Config accepted an invalid value.');
for name = ["mdemat.tools.environment.ConfigParser", "mdemat.tools.analysis.Analyzer", ...
        "mdemat.tools.analysis.Metric", "mdemat.tools.visual.Dataset", ...
        "mdemat.tools.visual.Database", "mdemat.tools.visual.GraphPlotter", ...
        "mdemat.tools.analysis.weightedError", "mdemat.tools.timing.Timer"]
    assert(isempty(which(name)),'mdemat:build:OldToolVisible','A removed tool class is packaged.');
end
fig = figure('Visible','off');
guard = onCleanup(@() delete(fig)); %#ok<NASGU>
ax = axes(fig);
graph = mdemat.data.mesh.Graph([0,0;1,1],[1,2]);
handles = plot(graph,axes=ax);
assert(isgraphics(handles.Edges) && isequal(handles.Edges.Parent,ax), ...
    'mdemat:build:GraphInvariant','The packaged external plot method failed.');
view = struct('Values',[2,4,6;3,6,9],'ComponentIds',[2,4],'Axes',{{[0,1,2]}});
samples = struct('RUN',struct('Type','numeric','Data',struct('u',view)));
selected = mdemat.tools.visual.selectComponent(samples,'u',4);
assert(isequal(selected.RUN.Data.u.Values,[3,6,9]), ...
    'mdemat:build:SamplesInvariant','Original component IDs were lost.');
artists = mdemat.tools.visual.plotField(ax,samples,'u',4,mode="1d");
assert(~isempty(artists) && all(isgraphics(artists)), ...
    'mdemat:build:SamplesInvariant','Installed field rendering failed.');
for name = {'Configuration.m','Graphics.m','Assessment.m'}
    assert(isfile(fullfile(toolboxRoot,'doc',name{1})), ...
        'mdemat:build:HelpMissing','The tools migration guide is missing.');
end
end

function uninstallOwned(installed)
matlab.addons.toolbox.uninstallToolbox(installed);
end

function removeOwnedRegistration()
% Called only after selecting our disposable installation root. The same GUID
% was ruled out in the caller's original registration before changing it.
registrations = matlab.addons.toolbox.installedToolboxes;
for index = 1:numel(registrations)
    if strcmp(registrations(index).Guid, '41e45d88-ea42-4d8f-bcac-46759e59ed38')
        matlab.addons.toolbox.uninstallToolbox(registrations(index));
    end
end
assertInstallationRemoved('');
end

function assertInstallationRemoved(source)
registrations = matlab.addons.toolbox.installedToolboxes;
assert(isempty(registrations) || ...
    ~any(strcmp({registrations.Guid}, '41e45d88-ea42-4d8f-bcac-46759e59ed38')), ...
    'mdemat:build:UninstallFailed', 'The disposable toolbox registration remains.');
assert(isempty(source) || ~isfile(source), 'mdemat:build:UninstallFailed', ...
    'The disposable installed source remains after uninstall.');
end

function assertSettingRestored(setting, hadTemporary, previous)
if hadTemporary
    restored = hasTemporaryValue(setting) && isequal(setting.TemporaryValue, previous);
else
    restored = ~hasTemporaryValue(setting);
end
assert(restored, 'mdemat:build:SettingRestoreFailed', ...
    'The caller add-on installation setting was not restored.');
end

function measurements = verifyDistMesh(toolboxRoot)
originalPath = path;
restorePath = onCleanup(@() path(originalPath)); %#ok<NASGU>
folder = fullfile(toolboxRoot, 'examples', 'geometry');
addpath(folder);
assert(strcmp(which('ex04_nonuniform_trimesh'), ...
    fullfile(folder, 'ex04_nonuniform_trimesh.m')), ...
    'mdemat:build:DistMeshPath', 'The mesher example must resolve inside the installed artifact.');
assert(isfile(fullfile(toolboxRoot, 'third_party', 'distmesh', 'LICENSE')) && ...
    isfile(fullfile(toolboxRoot, 'third_party', 'distmesh', 'src', 'distmesh2d.m')), ...
    'mdemat:build:DistMeshPayload', 'The bundled mesher/license is missing.');
scopedPath = path;
originalRng = rng;
figures = findall(groot, 'Type', 'figure');
assert(isempty(which('distmesh2d')), 'mdemat:build:DistMeshPath', ...
    'Vendor functions must stay outside the default installed MATLAB path.');
result = ex04_nonuniform_trimesh(ElementSize=0.3, MakePlots=false, Seed=7);
assert(abs(result.physicalArea-pi)<0.1 && result.quality.minQuality>0.2 && ...
    result.refined.NElements==4*result.mesh.NElements, ...
    'mdemat:build:DistMeshInvariant', 'Installed circle mesh/refinement failed.');
assert(strcmp(path, scopedPath) && isequal(rng, originalRng) && ...
    isequal(findall(groot, 'Type', 'figure'), figures) && isempty(which('distmesh2d')), ...
    'mdemat:build:DistMeshCleanup', 'The mesher changed caller path/RNG/figures.');
measurements = struct('CircleAreaError', abs(result.physicalArea-pi), ...
    'MinimumQuality', result.quality.minQuality, 'Triangles', result.mesh.NElements);
end

function restoreSetting(setting, hadTemporary, previous)
if hadTemporary, setting.TemporaryValue = previous;
elseif hasTemporaryValue(setting), clearTemporaryValue(setting); end
end

function resetSearchPath()
% restoredefaultpath is a script; run it outside the nested static workspace.
restoredefaultpath;
end

function result = runGuide(toolboxRoot)
original = path;
guard = onCleanup(@() path(original)); %#ok<NASGU>
addpath(fullfile(toolboxRoot, 'doc'));
assert(contains(help('GettingStarted'), 'exact continuous'), ...
    'mdemat:build:HelpMissing', 'The guide native help is unavailable.');
result = GettingStarted();
end

function measurements = verifyEvolutionFamilies()
% Independent stationary solutions exercise installed family entry points.
figures = findall(groot, 'Type', 'figure');
files = folderEntries(pwd);
space = mdemat.common.fem.FiniteElementSpace( ...
    mdemat.data.mesh.UniformGrid([0, 1], 2), ...
    mdemat.common.fem.BH1OrthotopeElement.modal(1, 2));
nodes = space.Element.Volume.Nodes;
[advection, transportReport] = mdemat.solver.hyperbolic.solveAdvection( ...
    mdemat.solver.hyperbolic.advectionProblem(1, 1), space, ...
    mdemat.common.ode.Exrk4Integrator(0.01), dt=0.005);
[diffusion, diffusionReport] = mdemat.solver.parabolic.solveDiffusion( ...
    mdemat.solver.parabolic.diffusionProblem(1, 1), space, ...
    mdemat.common.ode.Sdirk2Integrator(0.01), dt=0.005);
transportError = max(abs(space.eval(nodes, advection.U)-1), [], 'all');
diffusionError = max(abs(space.eval(nodes, diffusion.U)-1), [], 'all');
assert(transportError < 2e-12 && diffusionError < 2e-12 && ...
    abs(transportReport.FinalMass-1) < 2e-12 && abs(diffusionReport.FinalMass-1) < 2e-12, ...
    'mdemat:build:EvolutionInvariant', 'Installed periodic advection/diffusion equilibrium failed.');

spatial = mdemat.common.fem.FiniteElementSpace( ...
    mdemat.data.mesh.UniformGrid([0, 1], 2), ...
    mdemat.common.fem.BH1OrthotopeElement.nodal(1, 2));
macro = mdemat.common.fem.SpectralSpace( ...
    mdemat.common.fem.L2SphereElement.modal(1, 1, reduction='topology'));
micro = mdemat.common.fem.SpectralSpace( ...
    mdemat.common.fem.L2SphereElement.nodal(1, 2, reduction='topology'));
angular = mdemat.common.fem.SumSpace(macro, micro);
model = mdemat.application.radiation.problem(spatial, angular, 1, 0.02, 0.01, ...
    config=struct('scattering', 0.7, 'absorption', 0.2, 'source', 0.2, ...
    'decomposition', 'partial', 'tOdeIntName', 'sdirk2'));
[radiation, radiationScheme] = mdemat.application.radiation.solve(model);
radiationError = max(abs(radiation.state.distribution(0, micro.Element.Volume.Nodes)-1), [], 'all');
assert(radiationError < 2e-10 && abs(radiation.diagnostics.mass(end)-2) < 3e-10, ...
    'mdemat:build:RadiationInvariant', 'Installed radiation equilibrium or angular normalization failed.');
assertObserversEmpty(radiationScheme);

temperature = 0.7;
spatial = mdemat.common.fem.FiniteElementSpace( ...
    mdemat.data.mesh.UniformGrid([0, 2*pi], 2), ...
    mdemat.common.fem.BH1OrthotopeElement.nodal(1, 2));
velocity = mdemat.common.fem.SpectralSpace( ...
    mdemat.common.fem.L2OrthotopeElement.hermite(1, 3, T=temperature));
initial = @(x, nh) [ones(1, size(x, 2)); zeros(nh-1, size(x, 2))];
velocities = [-1, 0, 1]*sqrt(temperature);
exact = exp(-velocities.^2/(2*temperature))/sqrt(2*pi*temperature);
plasmaErrors = zeros(1, 2);
methods = {'ap', 'dghs'};
for index = 1:numel(methods)
    model = mdemat.application.plasma.problem(spatial, velocity, initial, 0.02, 0.01, ...
        method=methods{index}, config=struct('tau', 0.4, 'tau0', 0.4));
    [plasma, plasmaScheme] = mdemat.application.plasma.solve(model);
    plasmaErrors(index) = max(abs(plasma.state.distribution(0, velocities)-exact), [], 'all');
    assert(plasmaErrors(index) < 5e-11 && ...
        abs(plasma.diagnostics.mass(end)-2*pi) < 3e-11 && ...
        abs(plasma.diagnostics.kineticEnergy(end)-pi*temperature) < 3e-11, ...
        'mdemat:build:PlasmaInvariant', 'Installed warm-Hermite equilibrium/normalization failed.');
    assertObserversEmpty(plasmaScheme);
end

weights = [0.5; 0.5];
factors = mdemat.data.tensor.MatrixFactors(ones(4, 1), 1, ones(2, 1), ...
    LeftMetric=speye(4)/4, RightMetric=spdiags(weights, 0, 2, 2));
[transport, lowRankReport] = mdemat.application.radiation.solveLowRankTransport( ...
    factors, [-1, 1], weights, [0, 0.02], 0.01, Scattering=0.7, MaxRank=1);
lowRankError = max(abs(transport.Field.materialize()-1), [], 'all');
assert(lowRankError < 2e-12 && abs(lowRankReport.FinalMass-1) < 2e-12, ...
    'mdemat:build:LowRankInvariant', 'Installed low-rank transport equilibrium/mass failed.');
assert(isequal(findall(groot, 'Type', 'figure'), figures), ...
    'mdemat:build:UnexpectedFigures', 'Pure installed family APIs created figures.');
assert(isequal(folderEntries(pwd), files), 'mdemat:build:UnexpectedOutput', ...
    'Pure installed family APIs created files or directories.');
measurements = struct('AdvectionMaxError', transportError, 'DiffusionMaxError', diffusionError, ...
    'RadiationMaxError', radiationError, 'WarmPlasmaMaxErrors', plasmaErrors, ...
    'LowRankMaxError', lowRankError, 'LowRankMass', lowRankReport.FinalMass);
end

function entries = folderEntries(folder)
items = dir(fullfile(folder, '**', '*'));
entries = sort(string(fullfile({items.folder}, {items.name})));
end

function assertObserversEmpty(scheme)
assert(isempty(scheme.Observers.timer) && isempty(scheme.Observers.visualizer) && ...
    scheme.Observation.Analysis.PointsPerDimension == 0, 'mdemat:build:UnexpectedObserver', ...
    'A pure installed application constructed an observer.');
end

