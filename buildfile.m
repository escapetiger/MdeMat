function plan = buildfile
%BUILDFILE Analyze, test, assess, package and verify with MATLAB buildtool.
plan = buildplan(localfunctions);
plan.DefaultTasks = ["analyzer", "test"];
plan("package").Dependencies = ["analyzer", "test"];
end

function analyzerTask(~)
% Check every distributable MATLAB file and save all analyzer findings.
utility("analyzeSource");
end

function checkTask(~)
% Compatibility spelling for the analyzer task.
utility("analyzeSource");
end

function testTask(~)
% Include migrated v4 tests and enforce real coverage.
utility("runToolboxTests", Coverage=true);
end

function numericalTask(~)
% Run independently selected numerical and solver tests.
utility("runToolboxTests", Suite="numerical");
end

function symbolicTask(~)
% Require the optional symbolic capability explicitly.
utility("runToolboxTests", Suite="symbolic");
end

function packageTask(~)
% Stage reviewed payload and create an installable .mltbx.
utility("packageToolboxArtifact", Overwrite=true);
end

function verifyPackageTask(~)
% Verify an existing artifact in a fresh external MATLAB process.
root = fileparts(mfilename('fullpath'));
metadata = jsondecode(fileread(fullfile(root, 'toolbox', 'resources', 'toolbox.json')));
utility("verifyPackageArtifact", fullfile(root, 'release', ...
    "MdeMat-" + metadata.version + ".mltbx"));
end

function utility(name, varargin)
root = fileparts(mfilename('fullpath'));
original = path;
guard = onCleanup(@() path(original)); %#ok<NASGU>
addpath(fullfile(root, 'buildUtilities'));
feval(name, varargin{:});
end
