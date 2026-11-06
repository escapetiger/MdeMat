function report = verifyNativeProject()
%VERIFYNATIVEPROJECT Open the real project after copying to a native temp drive.
%   This checks portable project metadata without relying on a WSL UNC checkout.
root = fileparts(fileparts(mfilename('fullpath')));
work = tempname;
mkdir(work);
original = path;
originalFolder = pwd;
cleanup = onCleanup(@() cleanProject(work, original, originalFolder)); %#ok<NASGU>
folders = {'toolbox', 'resources', 'tests', 'benchmarks', 'docs', ...
    'buildUtilities', 'third_party', '.github'};
for index = 1:numel(folders)
    copyfile(fullfile(root, folders{index}), fullfile(work, folders{index}));
end
files = dir(root);
files = files(~[files.isdir]);
for index = 1:numel(files), copyfile(fullfile(root, files(index).name), work); end
project = matlab.project.loadProject(work);
for retired = {'+analysis', '+timing', '+except'}
    assert(~isfolder(fullfile(work, 'toolbox', '+mdemat', '+tools', retired{1})), ...
        'mdemat:build:RetiredPackageRestored', ...
        'Opening the project recreated a retired tools package.');
end
assert(strcmpi(char(project.RootFolder), work), ...
    'mdemat:build:ProjectRoot', 'Project root did not follow its copied location.');
resolved = which('mdemat.common.fdm.secondDerivative');
assert(startsWith(resolved, work, IgnoreCase=true), ...
    'mdemat:build:ProjectPath', 'Project source path is not portable.');
members = string({project.Files.Path});
assert(all(isfile(members) | isfolder(members)), ...
    'mdemat:build:ProjectMembership', 'The project contains stale or uncopied entries.');
projectPaths = project.ProjectPath;
assert(numel(projectPaths) == 1 && ...
    strcmpi(char(projectPaths(1).File), fullfile(work, 'toolbox')), ...
    'mdemat:build:ProjectPath', 'Only the toolbox root may enter the project path.');
report = struct('passed', true, 'release', version('-release'), ...
    'nativeRoot', work, 'resolvedSource', resolved, 'members', numel(members), ...
    'pathFolders', 1, 'retiredToolPackagesAbsent', true, 'cleanupComplete', false);
cleanProject(work, original, originalFolder);
report.cleanupComplete = ~isfolder(work);
writeBuildEvidence(fullfile(root, 'release', 'native-project.json'), report);
end

function cleanProject(work, originalPath, originalFolder)
project = matlab.project.currentProject;
if ~isempty(project) && strcmpi(char(project.RootFolder), work), close(project); end
cd(originalFolder);
path(originalPath);
if isfolder(work), rmdir(work, 's'); end
end
