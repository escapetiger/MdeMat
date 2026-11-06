function report = packageToolboxArtifact(options)
%PACKAGETOOLBOXARTIFACT Stage an explicit payload and create a real .mltbx.
%   Toolbox/+mdemat, resources, doc, examples and reviewed vendor files are eligible. Repository
%   docs/PDFs, tests, build tools, Git and old namespaces cannot enter staging.
%   LGPL notices and full editable runtime source remain in the distribution.
arguments
    options.Overwrite (1, 1) logical = false
end
root = fileparts(fileparts(mfilename('fullpath')));
metadata = jsondecode(fileread(fullfile(root, 'toolbox', 'resources', 'toolbox.json')));
outputFolder = fullfile(root, 'release');
if ~isfolder(outputFolder), mkdir(outputFolder); end
output = fullfile(outputFolder, "MdeMat-" + metadata.version + ".mltbx");
if isfile(output) && ~options.Overwrite
    error('mdemat:build:ArtifactExists', 'Refusing to overwrite %s.', output);
end
work = tempname;
stage = fullfile(work, 'MdeMat');
mkdir(stage);
originalPath = path;
cleanup = onCleanup(@() cleanStage(work, originalPath)); %#ok<NASGU>
source = fullfile(root, 'toolbox');
allowedRoots = {'+mdemat', 'resources', 'doc', 'examples'};
extensions = {'.m', '.mat', '.json', '.csv', '.txt', '.md', '.xml', '.html', '.css', '.png'};
files = dir(fullfile(source, '**', '*'));
relative = strings(0, 1);
for index = 1:numel(files)
    if files(index).isdir, continue; end
    absolute = fullfile(files(index).folder, files(index).name);
    name = erase(absolute, [source, filesep]);
    tokens = strsplit(name, filesep);
    [~, stem, extension] = fileparts(name);
    forbidden = any(ismember(tokens, {'docs', 'tests', '+tests', '.git', ...
        'benchmarks', 'buildUtilities', 'outputs', 'release', 'legacy'}));
    if ~ismember(tokens{1}, allowedRoots) || ~ismember(lower(extension), extensions) || ...
            forbidden || startsWith(stem, 'Test')
        error('mdemat:build:ForbiddenPayload', 'Unreviewed toolbox payload: %s.', name);
    end
    copyInto(absolute, fullfile(stage, name));
    relative(end+1, 1) = string(name); %#ok<AGROW>
end
provenanceFile = fullfile(root, 'buildUtilities', 'symbolic', 'basis-provenance.json');
provenance = jsondecode(fileread(provenanceFile));
resources = provenance.resources;
if iscell(resources), resources = [resources{:}]; end
required = strcat({resources.name}, '.mat');
if isempty(required), error('mdemat:build:MissingResource', 'The basis export manifest is empty.'); end
for index = 1:numel(required)
    resource = fullfile(stage, 'resources', 'basis', required{index});
    if ~isfile(resource)
        error('mdemat:build:MissingResource', 'Missing numeric resource %s.', required{index});
    end
    loaded = load(resource, 'metadata');
    if ~isfield(loaded, 'metadata') || ~plainResource(loaded.metadata)
        error('mdemat:build:OpaqueResource', 'Resource %s is not portable value data.', required{index});
    end
end
if ~isfolder(fullfile(stage, '+mdemat', '+internal', '+lebedev'))
    error('mdemat:build:VendorMissing', 'The licensed Lebedev runtime source is missing.');
end
copyInto(fullfile(root, 'LICENSE'), fullfile(stage, 'LICENSE'));
relative(end+1) = "LICENSE";
notices = {'NOTICE', 'README.md', 'COPYING', 'COPYING.LESSER', 'provenance.json', 'upstream.html'};
for index = 1:numel(notices)
    destination = fullfile('third_party', 'lebedev', notices{index});
    copyInto(fullfile(root, destination), fullfile(stage, destination));
    relative(end+1) = string(destination); %#ok<AGROW>
end
% DistMesh stays outside the default MATLAB path. Its complete runtime source
% and license are explicit payload; the integration example scopes src itself.
distmesh = fullfile(root, 'third_party', 'distmesh');
distmeshRecords = {'LICENSE', 'README.md', 'NOTICE', 'INTEGRATION.md', ...
    'provenance.json', 'CITATION.bib', 'CITATION.cff'};
distmeshSource = dir(fullfile(distmesh, 'src', '*.m'));
if isempty(distmeshSource)
    error('mdemat:build:VendorMissing', 'The bundled DistMesh runtime source is missing.');
end
distmeshRecords = [distmeshRecords, cellfun(@(name) fullfile('src', name), ...
    {distmeshSource.name}, UniformOutput=false)];
for index = 1:numel(distmeshRecords)
    destination = fullfile('third_party', 'distmesh', distmeshRecords{index});
    copyInto(fullfile(distmesh, distmeshRecords{index}), fullfile(stage, destination));
    relative(end+1) = string(destination); %#ok<AGROW>
end
guide = fullfile(stage, 'doc', 'GettingStarted.m');
if ~isfile(guide), error('mdemat:build:MissingHelp', 'GettingStarted.m is required.'); end
addpath(stage);
addpath(fullfile(stage, 'third_party', 'distmesh', 'src'), '-begin');
matlabFiles = fullfile(stage, relative(endsWith(relative, '.m')));
[dependencies, products] = matlab.codetools.requiredFilesAndProducts(matlabFiles);
[dependencies, aliases] = portableDependencies(dependencies, stage, root, relative);
productNames = string({products.Name});
if any(~ismember(productNames, ["MATLAB", "Symbolic Math Toolbox"]))
    error('mdemat:build:UnexpectedProductDependency', ...
        'Unreviewed product dependency: %s.', join(productNames, ', '));
end
dependencyReview = struct('CoreProduct', 'MATLAB', ...
    'OptionalSourceProduct', 'Symbolic Math Toolbox: explicit common.symbolic features only', ...
    'BundledMeshGenerator', 'DistMesh 1.2: MIT, pure MATLAB, scoped example path', ...
    'CachedSourceAliases', aliases, ...
    'Availability', 'Static discovery is not evidence of license availability');
manifest = struct('schemaVersion', 1, 'name', metadata.name, 'version', metadata.version, ...
    'identifier', metadata.identifier, 'files', sort(relative), ...
    'dependencyProducts', products, 'dependencyFiles', dependencies, ...
    'optionalCapabilities', {{'Symbolic Math Toolbox: symbolic inputs/generation only', ...
    'DistMesh 1.2: included source and license; loaded explicitly by geometry example'}}, ...
    'dependencyReview', dependencyReview, ...
    'platformEvidence', 'Windows MATLAB R2025b; other platforms unverified');
writeBuildEvidence(fullfile(stage, 'manifest.json'), manifest);
relative(end+1) = "manifest.json";
opts = matlab.addons.toolbox.ToolboxOptions(stage, metadata.identifier);
opts.ToolboxName = metadata.name;
opts.ToolboxVersion = metadata.version;
opts.AuthorName = metadata.author;
opts.Summary = metadata.summary;
opts.Description = 'Numerical methods, DG PDE solvers and kinetic applications; see GettingStarted.';
opts.MinimumMatlabRelease = metadata.minimumRelease;
opts.SupportedPlatforms = struct('Win64', true, 'Glnxa64', false, ...
    'Maci64', false, 'MatlabOnline', false);
opts.ToolboxFiles = fullfile(stage, relative);
opts.ToolboxMatlabPath = string(stage);
opts.ToolboxGettingStartedGuide = string(guide);
opts.OutputFile = output;
matlab.addons.toolbox.packageToolbox(opts);
report = struct('artifact', output, 'release', version('-release'), ...
    'artifactSHA256', fileSHA256(output), 'manifest', manifest, ...
    'payloadFiles', numel(relative), 'packaged', isfile(output));
writeBuildEvidence(fullfile(outputFolder, 'package-report.json'), report);
end

function [paths, aliases] = portableDependencies(files, stage, root, payload)
% No development checkout or deleted staging location becomes public metadata.
% Loaded class definitions can make static analysis report their checkout
% paths. Accept only explicit payload members with byte-identical staged
% copies; never clear the caller's classes or accept arbitrary external code.
paths = string(files);
stagePrefix = replace(string(stage), '\', '/') + "/";
matlabPrefix = replace(string(matlabroot), '\', '/') + "/";
sourcePrefix = replace(string(fullfile(root, 'toolbox')), '\', '/') + "/";
rootPrefix = replace(string(root), '\', '/') + "/";
eligible = replace(string(payload), '\', '/');
aliases = struct('File', {}, 'SHA256', {});
for index = 1:numel(paths)
    path = replace(paths(index), '\', '/');
    if startsWith(path, stagePrefix, IgnoreCase=ispc)
        paths(index) = "toolbox/" + extractAfter(path, strlength(stagePrefix));
    elseif startsWith(path, matlabPrefix, IgnoreCase=ispc)
        paths(index) = "MATLABROOT/" + extractAfter(path, strlength(matlabPrefix));
    else
        relative = "";
        if startsWith(path, sourcePrefix, IgnoreCase=ispc)
            relative = extractAfter(path, strlength(sourcePrefix));
        elseif startsWith(path, rootPrefix, IgnoreCase=ispc)
            relative = extractAfter(path, strlength(rootPrefix));
        end
        if ~ismember(relative, eligible) || ~isfile(fullfile(stage, relative))
            error('mdemat:build:UnstagedDependency', ...
                'A required file is outside the explicit staged payload: %s.', files{index});
        end
        hash = fileSHA256(files{index});
        if ~strcmp(hash, fileSHA256(fullfile(stage, relative)))
            error('mdemat:build:ChangedDependency', ...
                'A cached source dependency differs from its staged copy: %s.', files{index});
        end
        paths(index) = "toolbox/" + relative;
        aliases(end+1) = struct('File', relative, 'SHA256', hash); %#ok<AGROW>
    end
end
paths = sort(unique(paths));
end

function hash = fileSHA256(filename)
stream = fopen(filename, 'rb');
if stream < 0, error('mdemat:build:ArtifactRead', 'Cannot read the packaged artifact.'); end
guard = onCleanup(@() fclose(stream)); %#ok<NASGU>
digest = java.security.MessageDigest.getInstance('SHA-256');
information = dir(filename);
remaining = information.bytes;
while remaining > 0
    requested = min(remaining, 1024*1024);
    bytes = fread(stream, requested, '*uint8');
    if numel(bytes) ~= requested
        error('mdemat:build:ArtifactRead', 'The packaged artifact could not be read completely.');
    end
    digest.update(typecast(bytes, 'int8'));
    remaining = remaining-numel(bytes);
end
[message, number] = ferror(stream);
if number ~= 0, error('mdemat:build:ArtifactRead', '%s', message); end
hash = lower(reshape(dec2hex(typecast(digest.digest(), 'uint8'), 2).', 1, []));
end

function copyInto(source, destination)
parent = fileparts(destination);
if ~isfolder(parent), mkdir(parent); end
[ok, message] = copyfile(source, destination);
if ~ok, error('mdemat:build:StageCopy', '%s', message); end
end

function cleanStage(work, originalPath)
path(originalPath);
if isfolder(work), rmdir(work, 's'); end
end

function valid = plainResource(value)
if isnumeric(value) || islogical(value) || ischar(value), valid = true;
elseif iscell(value), valid = all(cellfun(@plainResource, value), 'all');
elseif isstruct(value), valid = all(cellfun(@plainResource, struct2cell(value)), 'all');
else, valid = false; end
end
