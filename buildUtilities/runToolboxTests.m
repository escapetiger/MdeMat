function report = runToolboxTests(options)
%RUNTOOLBOXTESTS Run explicit suites, with real optional statement coverage.
%   Suite="all" includes retained/migrated tests, unit, regression and numerical.
%   Optional symbolic tests have a separate mandatory-capability job.
%   Coverage excludes only third-party vendor code and optional symbolic source;
%   the excluded paths are written into the report. The project target is 80%.
arguments
    options.Suite (1, 1) string {mustBeMember(options.Suite, ["all", "numerical", "symbolic"])} = "all"
    options.Coverage (1, 1) logical = false
    options.MinimumCoverage (1, 1) double {mustBeNonnegative, mustBeLessThanOrEqual(options.MinimumCoverage, 1)} = 0.8
end
root = fileparts(fileparts(mfilename('fullpath')));
original = path;
guard = onCleanup(@() path(original)); %#ok<NASGU>
addpath(fullfile(root, 'toolbox'));
switch options.Suite
    case "all"
        folders = {'migrated', 'unit', 'regression', 'numerical'};
    case "numerical"
        folders = {'numerical', fullfile('regression', 'solvers')};
    case "symbolic"
        if ~license('test', 'Symbolic_Toolbox') || exist('sym', 'file') ~= 2
            error('mdemat:build:SymbolicRequired', 'Symbolic Math Toolbox is required by this job.');
        end
        folders = {fullfile('optional', 'symbolic')};
end
suite = [];
for index = 1:numel(folders)
    folder = fullfile(root, 'tests', folders{index});
    if isfolder(folder)
        part = matlab.unittest.TestSuite.fromFolder(folder, IncludingSubfolders=true);
        if isempty(suite)
            suite = reshape(part, 1, []);
        else
            suite = [suite, reshape(part, 1, [])]; %#ok<AGROW>
        end
    end
end
if isempty(suite), error('mdemat:build:EmptySuite', 'No tests found for %s.', options.Suite); end
output = fullfile(root, 'outputs', 'checks');
if ~isfolder(output), mkdir(output); end
runner = matlab.unittest.TestRunner.withTextOutput;
runner.addPlugin(matlab.unittest.plugins.DiagnosticsRecordingPlugin);
runner.addPlugin(matlab.unittest.plugins.XMLPlugin.producingJUnitFormat( ...
    fullfile(output, options.Suite + "-tests.xml")));
excluded = strings(0, 1);
coverageFile = fullfile(output, options.Suite + "-coverage.xml");
if options.Coverage
    files = dir(fullfile(root, 'toolbox', '+mdemat', '**', '*.m'));
    paths = string(fullfile({files.folder}, {files.name}));
    omit = contains(paths, fullfile('+internal', '+lebedev')) | ...
        contains(paths, [filesep, '+symbolic', filesep]);
    upstream = dir(fullfile(root, 'third_party', 'distmesh', 'src', '*.m'));
    excluded = [paths(omit), string(fullfile({upstream.folder}, {upstream.name}))];
    format = matlab.unittest.plugins.codecoverage.CoberturaFormat(coverageFile);
    plugin = matlab.unittest.plugins.CodeCoveragePlugin.forFile(paths(~omit), ...
        MetricLevel="statement", Producing=format);
    runner.addPlugin(plugin);
end
visibility = get(groot, 'DefaultFigureVisible');
graphicsGuard = onCleanup(@() set(groot, 'DefaultFigureVisible', visibility)); %#ok<NASGU>
set(groot, 'DefaultFigureVisible', 'off');
results = runner.run(suite);
report = struct('suite', options.Suite, 'release', version('-release'), ...
    'capturedUTC', char(datetime('now',TimeZone='UTC')), ...
    'passed', sum([results.Passed]), 'failed', sum([results.Failed]), ...
    'incomplete', sum([results.Incomplete]), 'coverageMeasured', options.Coverage, ...
    'statementCoverage', NaN, 'minimumCoverage', options.MinimumCoverage, ...
    'coverageExclusions', excluded, 'tests', string({results.Name}));
if options.Coverage
    document = xmlread(coverageFile);
    report.statementCoverage = str2double(char( ...
        document.getDocumentElement().getAttribute('line-rate')));
end
writeBuildEvidence(fullfile(output, options.Suite + "-tests.json"), report);
if any([results.Failed] | [results.Incomplete])
    error('mdemat:build:TestsFailed', 'The %s suite has failed/incomplete tests.', options.Suite);
end
if options.Coverage && (~isfinite(report.statementCoverage) || ...
        report.statementCoverage < options.MinimumCoverage)
    error('mdemat:build:CoverageBelowTarget', ...
        'Measured statement coverage %.2f%% is below %.2f%%.', ...
        100*report.statementCoverage, 100*options.MinimumCoverage);
end
end
