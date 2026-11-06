function initializeProject()
%INITIALIZEPROJECT Create/update the actual MATLAB project using native APIs.
root = fileparts(fileparts(mfilename('fullpath')));
projectFile = fullfile(root, 'MdeMat.prj');
if isfile(projectFile)
    project = matlab.project.loadProject(root);
else
    project = matlab.project.createProject(Folder=root, Name="MdeMat");
end
closeProject = onCleanup(@() close(project)); %#ok<NASGU>
% Migration removes source locations. Native project membership must describe
% the current checkout, including when this initializer is run repeatedly.
members = project.Files;
for index = 1:numel(members)
    filename = members(index).Path;
    if ~isfile(filename) && ~isfolder(filename)
        project.removeFile(filename);
    end
end
% Project metadata can recreate registered empty folders on open. Retire
% these former packages explicitly, after removing their stale file members.
retired = {fullfile('+except','private'), '+except', '+timing', '+analysis'};
for index = 1:numel(retired)
    folder = fullfile(root, 'toolbox', '+mdemat', '+tools', retired{index});
    if isfolder(folder)
        contents = dir(fullfile(folder, '**', '*'));
        if any(~[contents.isdir])
            error('mdemat:build:RetiredPackageNotEmpty', ...
                'Refusing to remove unexpected files from retired package %s.', folder);
        end
    end
    member = project.findFile(folder);
    if ~isempty(member), project.removeFile(folder); end
    if isfolder(folder), rmdir(folder, 's'); end
end
projectPaths = project.ProjectPath;
toolboxFolder = fullfile(root, 'toolbox');
for index = 1:numel(projectPaths)
    if ~strcmpi(char(projectPaths(index).File), toolboxFolder)
        project.removePath(projectPaths(index).File);
    end
end
project.addPath(fullfile(root, 'toolbox'));
folders = {'toolbox', 'tests', 'benchmarks', 'docs', 'buildUtilities', 'third_party', '.github'};
for index = 1:numel(folders)
    project.addFolderIncludingChildFiles(fullfile(root, folders{index}));
end
files = dir(root);
files = files(~[files.isdir]);
for index = 1:numel(files)
    if ~strcmp(files(index).name, 'MdeMat.prj')
        project.addFile(fullfile(root, files(index).name));
    end
end
end
