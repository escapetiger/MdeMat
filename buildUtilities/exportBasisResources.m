function report = exportBasisResources(options)
%EXPORTBASISRESOURCES Export symbolic metadata as portable numeric resources.
%   report = exportBasisResources() reads the development MAT files and writes
%   toolbox/resources/basis/*.mat. Each runtime file contains only a versioned
%   structure of numbers, character vectors and cells: dimensions, derivative
%   orders and self-contained anonymous-function text. No sym objects or saved
%   function workspaces are distributed. The runtime loader rebuilds handles
%   from this trusted generated text using str2func.
%
%   Symbolic expressions, variable names, source SHA-256 and normalization
%   descriptions remain in the development provenance JSON. Existing files
%   require overwrite=true. All samples and derivative columns are checked
%   against their source handles before any output file is written.
%
%   Optional names: sourceDirectory, destinationDirectory, provenanceFile,
%   overwrite. Exporting the original symbolic metadata needs Symbolic Math
%   Toolbox; loading or evaluating the exported runtime resources does not.

arguments
    options.sourceDirectory (1, 1) string = ""
    options.destinationDirectory (1, 1) string = ""
    options.provenanceFile (1, 1) string = ""
    options.overwrite (1, 1) logical = false
end

root = fileparts(fileparts(mfilename('fullpath')));
if options.sourceDirectory == ""
    options.sourceDirectory = fullfile(root, 'buildUtilities', 'symbolic', 'metadata');
    if ~isfolder(options.sourceDirectory)
        options.sourceDirectory = fullfile(root, 'buildUtilities', 'symbolic', ...
            'original-resources');
    end
end
if options.destinationDirectory == ""
    options.destinationDirectory = fullfile(root, 'toolbox', 'resources', 'basis');
end
if options.provenanceFile == ""
    options.provenanceFile = fullfile(root, 'buildUtilities', 'symbolic', ...
        'basis-provenance.json');
end

files = dir(fullfile(options.sourceDirectory, '*.mat'));
assert(~isempty(files), 'MdeMat:BasisExport:MissingSource', ...
    'No symbolic metadata found in %s.', options.sourceDirectory);
destinations = fullfile(options.destinationDirectory, {files.name});
assert(options.overwrite || (~any(isfile(destinations)) && ...
    ~isfile(options.provenanceFile)), 'MdeMat:BasisExport:FileExists', ...
    'An export destination exists. Set overwrite=true to replace it.');

resources = cell(1, numel(files));
provenance = cell(1, numel(files));
report = repmat(struct('name', '', 'nDims', 0, 'nCodims', 0, ...
    'maxDerivativeOrder', 0, 'sourceSha256', '', 'maxRoundtripError', 0), ...
    1, numel(files));

for index = 1:numel(files)
    sourceFile = fullfile(files(index).folder, files(index).name);
    loaded = load(sourceFile, 'metadata');
    source = loaded.metadata;
    [~, name] = fileparts(files(index).name);
    assert(source.nDims == 1 || size(source.handles, 2) == 1, ...
        'MdeMat:BasisExport:UnsupportedDerivatives', ...
        'Multivariate derivative columns need explicit multi-index metadata.');
    assert(size(source.handles, 1) == source.nCodims, ...
        'MdeMat:BasisExport:InvalidDimensions', ...
        'Handle count does not match the declared output dimension in %s.', sourceFile);

    metadata = struct('schemaVersion', 1, 'basisName', name, ...
        'nDims', double(source.nDims), 'nCodims', double(source.nCodims), ...
        'derivativeOrders', 0:(size(source.handles, 2)-1), ...
        'functionText', {cell(size(source.handles))});
    for entry = 1:numel(source.handles)
        original = source.handles{entry};
        if isempty(original)
            metadata.functionText{entry} = '';
            continue;
        end
        text = func2str(original);
        assert(startsWith(text, '@('), 'MdeMat:BasisExport:UnsupportedHandle', ...
            'A resource handle must be a self-contained anonymous function.');
        metadata.functionText{entry} = text;
    end
    handles = mdemat.common.approximation.CompiledFunction.compileResource( ...
        metadata, source.nDims, source.nCodims);
    samples = sampleInputs(source.nDims);
    maxError = 0;
    for entry = 1:numel(source.handles)
        original = source.handles{entry};
        if isempty(original)
            continue;
        end
        rebuilt = handles{entry};
        expected = original(samples{:});
        actual = rebuilt(samples{:});
        assert(isnumeric(actual) && all(isfinite(actual), 'all'), ...
            'MdeMat:BasisExport:NonfiniteSample', ...
            'Resource sample evaluation must return finite numeric values.');
        difference = max(abs(actual-expected), [], 'all');
        scale = max(1, max(abs(expected), [], 'all'));
        assert(difference <= 2e-14*scale, 'MdeMat:BasisExport:NonEquivalent', ...
            'Rebuilt numeric function differs from its source in %s.', sourceFile);
        maxError = max(maxError, difference);
    end
    resources{index} = metadata;
    digest = sha256(sourceFile);
    report(index) = struct('name', name, 'nDims', metadata.nDims, ...
        'nCodims', metadata.nCodims, ...
        'maxDerivativeOrder', metadata.derivativeOrders(end), ...
        'sourceSha256', digest, 'maxRoundtripError', maxError);
    provenance{index} = report(index);
    provenance{index}.variables = arrayfun(@char, source.variables(:), ...
        'UniformOutput', false);
    provenance{index}.expressions = arrayfun(@char, source.expressions(:), ...
        'UniformOutput', false);
    provenance{index}.normalization = normalizationDescription(name);
    provenance{index}.sourceFile = files(index).name;
end

if ~isfolder(options.destinationDirectory)
    mkdir(options.destinationDirectory);
end
for index = 1:numel(files)
    metadata = resources{index}; %#ok<NASGU>
    save(destinations{index}, 'metadata', '-v7');
end
parent = fileparts(options.provenanceFile);
if ~isfolder(parent)
    mkdir(parent);
end
manifest = struct('schemaVersion', 1, 'generator', 'exportBasisResources', ...
    'resources', {provenance});
file = fopen(options.provenanceFile, 'w', 'n', 'UTF-8');
assert(file ~= -1, 'MdeMat:BasisExport:CannotWrite', ...
    'Cannot write development provenance: %s.', options.provenanceFile);
closeFile = onCleanup(@() fclose(file)); %#ok<NASGU>
fwrite(file, jsonencode(manifest, PrettyPrint=true), 'char');
end

function samples = sampleInputs(nDims)
assert(nDims == 1 || nDims == 2, 'MdeMat:BasisExport:UnsupportedDimensions', ...
    'Only the existing one-dimensional and spherical bases are supported.');
if nDims == 1
    samples = {linspace(-0.37, 0.41, 19)};
else
    samples = {linspace(-0.8, 0.8, 19), linspace(0.1, 5.9, 19)};
end
end

function description = normalizationDescription(name)
if startsWith(name, 'unit_gauss_')
    description = 'Nodal Lagrange basis on [-1/2,1/2], original node order preserved.';
else
    switch name
        case 'chebyshev'
            description = 'Classical first-kind Chebyshev basis on [-1,1].';
        case 'monic_unit_legendre'
            description = 'Monic Legendre basis on [-1/2,1/2].';
        case 'normal_legendre'
            description = 'Orthonormal Legendre basis for unweighted measure on [-1,1].';
        case 'normal_hermite'
            description = ['Probabilist Hermite basis He_n/sqrt(n!), orthonormal for ', ...
                'exp(-x^2/2)/sqrt(2*pi) dx.'];
        case 'normal_fourier'
            description = 'Real orthonormal Fourier basis on [0,2*pi].';
        case 'normal_spherical_harmonic'
            description = ['Real orthonormal spherical harmonics for dmu*dphi; ', ...
                'Condon--Shortley phase cancelled, order [0,cos(phi),sin(phi),...].'];
        otherwise
            description = 'Preserved from the named symbolic metadata and expressions.';
    end
end
end

function digest = sha256(fileName)
file = fopen(fileName, 'rb');
assert(file ~= -1, 'MdeMat:BasisExport:CannotRead', 'Cannot read %s.', fileName);
closeFile = onCleanup(@() fclose(file)); %#ok<NASGU>
bytes = fread(file, Inf, '*uint8');
engine = java.security.MessageDigest.getInstance('SHA-256');
engine.update(bytes);
digest = lower(reshape(dec2hex(typecast(engine.digest(), 'uint8'), 2).', 1, []));
end
