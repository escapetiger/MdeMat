function report = verifyPackageArtifact(artifact)
%VERIFYPACKAGEARTIFACT Install/exercise an artifact in this fresh MATLAB process.
%   Invoke from a NEW external MATLAB batch process after packaging exits.
%   No nested MATLAB process is launched. Only the artifact and standalone
%   verifier enter a native temporary folder. Cached checkout code is refused.
arguments
    artifact {mustBeTextScalar, mustBeFile}
end
root = fileparts(fileparts(mfilename('fullpath')));
loaded = string(inmem('-completenames'));
if any(contains(loaded, [filesep, '+mdemat', filesep]))
    error('mdemat:build:FreshProcessRequired', ...
        'Start a new MATLAB batch process for installed-artifact verification.');
end
work = tempname;
mkdir(work);
originalFolder = pwd;
originalPath = path;
guard = onCleanup(@() cleanWork(work, originalFolder, originalPath)); %#ok<NASGU>
copyfile(artifact, fullfile(work, 'MdeMat.mltbx'));
copyfile(fullfile(root, 'buildUtilities', 'verifyInstalledPackage.m'), work);
cd(work);
addpath(work);
failure = [];
try
    verifyInstalledPackage();
catch exception
    failure = exception;
end
cd(originalFolder);
path(originalPath);
file = fullfile(work, 'installed-package.json');
if isfile(file), report = jsondecode(fileread(file));
else, report = struct('passed', false, 'reason', 'Verifier did not produce evidence.'); end
clear guard;
report.temporaryWorkRemoved = ~isfolder(work);
writeBuildEvidence(fullfile(root, 'release', 'installed-package.json'), report);
if ~isempty(failure), rethrow(failure); end
if ~report.passed
    error('mdemat:build:InstalledArtifactFailed', 'Installed check failed; inspect release evidence.');
end
end

function cleanWork(work, originalFolder, originalPath)
cd(originalFolder);
path(originalPath);
if isfolder(work), rmdir(work, 's'); end
end
