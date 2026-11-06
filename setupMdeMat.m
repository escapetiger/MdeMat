function environment = setupMdeMat(options)
%SETUPMDEMAT Add the explicit source root for a checkout.
%   setupMdeMat() enables the distributable toolbox.
%   No recursive paths or persistent MATLAB settings are changed.
root = fileparts(mfilename('fullpath'));
toolboxRoot = fullfile(root, 'toolbox');
sourceRoot = toolboxRoot;
if ~isfolder(sourceRoot)
    error('mdemat:environment:MissingSource', 'Missing source root: %s', sourceRoot);
end
addpath(sourceRoot);
environment = struct('Root', root, 'Toolbox', toolboxRoot, ...
    'Source', sourceRoot, 'Release', version('-release'));
end
