function results = runTests(options)
%RUNTESTS Execute required or explicitly optional toolbox test suites.
%   RESULTS = runTests() runs tests/ against toolbox/.
%   RESULTS = runTests(Suite="symbolic") tests the optional symbolic feature.
%   OutputFile writes a compact evidence JSON without changing suite outcomes.
arguments
    options.Suite (1, 1) string {mustBeMember(options.Suite, ["required", "symbolic", "all"])} = "required"
    options.OutputFile (1, 1) string = ""
    options.Pattern (1, 1) string = ""
    options.FailOnFailure (1, 1) logical = true
end
root = fileparts(mfilename('fullpath'));
oldPath = path;
restorePath = onCleanup(@() path(oldPath)); %#ok<NASGU>
setupMdeMat();
oldVisibility = get(groot, 'DefaultFigureVisible');
restoreGraphics = onCleanup(@() set(groot, 'DefaultFigureVisible', oldVisibility)); %#ok<NASGU>
set(groot, 'DefaultFigureVisible', 'off');
suite = [];
testRoot = fullfile(root, 'tests');
folders = ["unit", "regression", "migrated", "numerical"];
if options.Suite == "symbolic"
    folders = "optional/symbolic";
elseif options.Suite == "all"
    folders(end+1) = "optional/symbolic";
end
for folder = folders
    location = fullfile(testRoot, folder);
    if isfolder(location)
        suite = appendSuite(suite, matlab.unittest.TestSuite.fromFolder( ...
            location, IncludingSubfolders=true));
    end
end
if strlength(options.Pattern) > 0 && ~isempty(suite)
    selected = ~cellfun(@isempty, regexp({suite.Name}, char(options.Pattern), 'once'));
    suite = suite(selected);
end
if isempty(suite)
    error('mdemat:test:EmptySuite', 'No tests were discovered for suite %s.', options.Suite);
end
runner = matlab.unittest.TestRunner.withTextOutput;
runner.addPlugin(matlab.unittest.plugins.DiagnosticsRecordingPlugin);
results = runner.run(suite);
if strlength(options.OutputFile) > 0
    writeEvidence(options.OutputFile, results, options.Suite);
end
if options.FailOnFailure && any([results.Failed] | [results.Incomplete])
    error('mdemat:test:Failed', 'Required suite has failed or incomplete tests.');
end
end

function suite = appendSuite(suite, part)
if isempty(suite)
    suite = reshape(part, 1, []);
else
    suite = [suite, reshape(part, 1, [])];
end
end

function writeEvidence(file, results, mode)
records = repmat(struct('name', '', 'passed', false, 'failed', false, ...
    'incomplete', false, 'duration', 0, 'diagnostics', ''), 1, numel(results));
for index = 1:numel(results)
    result = results(index);
    diagnostics = '';
    if result.Failed || result.Incomplete
        recordsForTest = result.Details.DiagnosticRecord;
        parts = strings(1, numel(recordsForTest));
        for item = 1:numel(recordsForTest)
            parts(item) = string(recordsForTest(item).Report);
        end
        diagnostics = char(join(parts, newline));
    end
    records(index) = struct('name', char(result.Name), 'passed', result.Passed, ...
        'failed', result.Failed, 'incomplete', result.Incomplete, ...
        'duration', result.Duration, 'diagnostics', diagnostics);
end
report = struct('mode', char(mode), 'release', version('-release'), ...
    'capturedUTC', char(datetime('now', TimeZone='UTC')), ...
    'passed', sum([results.Passed]), 'failed', sum([results.Failed]), ...
    'incomplete', sum([results.Incomplete]), 'results', records);
parent = fileparts(file);
if strlength(parent) > 0 && ~isfolder(parent)
    mkdir(parent);
end
handle = fopen(file, 'w', 'n', 'UTF-8');
if handle < 0
    error('mdemat:test:EvidenceWrite', 'Cannot write %s.', file);
end
closeFile = onCleanup(@() fclose(handle)); %#ok<NASGU>
fprintf(handle, '%s\n', jsonencode(report, PrettyPrint=true));
end
