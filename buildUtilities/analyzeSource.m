function report = analyzeSource()
%ANALYZESOURCE Check every distributable MATLAB file; retain warning evidence.
root = fileparts(fileparts(mfilename('fullpath')));
files = dir(fullfile(root, 'toolbox', '**', '*.m'));
files = [files; dir(fullfile(root, 'third_party', 'distmesh', 'src', '*.m'))];
records = cell(numel(files), 1);
errors = 0;
warnings = 0;
for index = 1:numel(files)
    file = fullfile(files(index).folder, files(index).name);
    result = codeIssues(file, CodeAnalyzerConfiguration="factory");
    issues = result.Issues;
    errors = errors + sum(strcmpi(string(issues.Severity), "error"));
    warnings = warnings + sum(strcmpi(string(issues.Severity), "warning"));
    records{index} = struct('file', erase(file, [root, filesep]), ...
        'vendor', contains(file, fullfile('+internal', '+lebedev')) || ...
        contains(file, fullfile('third_party', 'distmesh', 'src')), ...
        'issues', table2struct(issues));
end
report = struct('release', version('-release'), ...
    'capturedUTC', char(datetime('now', TimeZone='UTC')), 'errors', errors, ...
    'warnings', warnings, 'warningReviewRequired', warnings > 0, 'files', {records});
writeBuildEvidence(fullfile(root, 'outputs', 'checks', 'analyzer.json'), report);
if errors > 0
    error('mdemat:build:AnalyzerErrors', '%d analyzer errors remain; inspect outputs/checks.', errors);
end
end
