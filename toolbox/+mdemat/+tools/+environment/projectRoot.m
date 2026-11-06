function root = projectRoot(startDirectory)
%PROJECTROOT Find this checkout, or search ancestors of an explicit directory.
%   Without arguments only the function's own checkout layout is accepted.
%   An installed copy nested inside another checkout cannot impersonate it.
%   Installed runtime resources are relative to toolboxRoot, never cwd.
arguments
    startDirectory {mustBeTextScalar, mustBeFolder} = fileparts(mfilename('fullpath'))
end
if nargin == 0
    root = canonicalDirectory(fileparts(mdemat.tools.environment.toolboxRoot()));
    if isCheckout(root), return; end
else
    root = canonicalDirectory(char(startDirectory));
    while true
        if isCheckout(root), return; end
        parent = fileparts(root);
        if isempty(parent) || strcmp(root, parent), break; end
        root = parent;
    end
end
error('mdemat:environment:NotCheckout', ...
    'No MdeMat checkout was found. Installed resources are relative to toolboxRoot.');
end

function tf = isCheckout(root)
tf = isfile(fullfile(root, 'MdeMat.prj')) && ...
    isfile(fullfile(root, 'setupMdeMat.m')) && isfolder(fullfile(root, 'toolbox'));
end

function directory = canonicalDirectory(directory)
[ok, attributes] = fileattrib(directory);
if ~ok || ~isstruct(attributes) || ~isfield(attributes, 'Name')
    error('mdemat:environment:PathUnavailable', ...
        'Cannot resolve source directory: %s', directory);
end
directory = attributes.Name;
end
