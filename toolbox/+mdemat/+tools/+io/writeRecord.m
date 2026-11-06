function writeRecord(file, data, options)
%WRITERECORD Save versioned value data without object serialization.
%   writeRecord(FILE, DATA) accepts numeric/logical arrays, text, cells and
%   structs. FILE ends in .mat or .json. Classes, shapes, sparsity, complex
%   values and non-finite values are retained. Sparse arrays support double
%   and logical. Objects/function handles fail.
%
%   Existing files are refused unless overwrite=true. Missing parent folders
%   require createDirectory=true. A temporary file is committed with Java NIO
%   exclusive move semantics; this requires the MATLAB JVM. No crash-atomic
%   storage guarantee is made. MAT stores primitive arrays and a UTF-8 schema.
%
%   See also mdemat.tools.io.readRecord.
%   limits overrides the same resource policy documented by readRecord.
arguments
    file {mustBeTextScalar}
    data
    options.overwrite (1, 1) logical = false
    options.createDirectory (1, 1) logical = false
    options.limits (1,1) struct = struct()
end
limits = recordLimits(options.limits);
file = char(file);
[directory, ~, extension] = fileparts(file);
extension = lower(extension);
if ~ismember(extension, {'.mat', '.json'})
    error('mdemat:io:InvalidFormat', 'Record filename must end in .mat or .json.');
end
if isempty(directory), directory = pwd; file = fullfile(directory, file); end
if isfile(file) && ~options.overwrite
    error('mdemat:io:FileExists', 'Refusing to overwrite %s.', file);
end
if ~usejava('jvm')
    error('mdemat:io:ExclusiveCommitUnavailable', 'Exclusive record creation requires the JVM.');
end
[tree, arrays] = packRecord(data,limits);
descriptions = struct();
names = fieldnames(arrays);
for index = 1:numel(names)
    value = arrays.(names{index});
    description = whos('value');
    descriptions.(names{index}) = struct('Type',class(value),'Bytes',description.bytes);
end
validateRecordTree(tree,descriptions,limits);
if ~isfolder(directory)
    if ~options.createDirectory
        error('mdemat:io:DirectoryNotFound', 'Output directory does not exist: %s', directory);
    end
    [created, message] = mkdir(directory);
    if ~created, error('mdemat:io:WriteFailed', '%s', message); end
end
temporary = [tempname(directory), extension];
guard = onCleanup(@() removeTemporary(temporary));
identity = mdemat.tools.environment.versionInfo();
if identity.NumericRecordSchemaVersion ~= 1
    error('mdemat:io:UnsupportedSchema','This writer implements numeric record schema v1.');
end
envelope = struct('format', 'mdemat.numericRecord', ...
    'schemaVersion', identity.NumericRecordSchemaVersion, ...
    'toolboxVersion', identity.ToolboxVersion, 'tree', tree);
try
    if strcmp(extension, '.mat')
        container = arrays;
        container.manifest = unicode2native(jsonencode(envelope), 'UTF-8');
        % -v7 skips unsupported >=2 GiB variables with a warning. Refuse
        % before save, then verify the temporary inventory before commit.
        fields = fieldnames(container);
        for index = 1:numel(fields)
            value = container.(fields{index}); %#ok<NASGU> whos inspects this variable by name.
            description = whos('value');
            if description.bytes >= 2^31
                error('mdemat:io:MatSizeLimit','MAT v7 variables must be smaller than 2 GiB.');
            end
        end
        save(temporary, '-struct', 'container', '-v7');
        inventory = whos('-file',temporary);
        if ~isequal(sort({inventory.name}.'),sort(fields))
            error('mdemat:io:IncompleteWrite','Temporary MAT variable inventory is incomplete.');
        end
        for index = 1:numel(inventory)
            value = container.(inventory(index).name);
            if ~strcmp(inventory(index).class,class(value)) || ...
                    ~isequal(inventory(index).size,size(value)) || ...
                    inventory(index).sparse ~= issparse(value) || inventory(index).complex ~= ~isreal(value)
                error('mdemat:io:IncompleteWrite','Temporary MAT variable descriptor differs.');
            end
        end
    else
        metadataBytes = numel(unicode2native(jsonencode(envelope),'UTF-8'));
        if metadataBytes+minimumJsonPayloadBytes(arrays) > limits.MaxFileBytes
            error('mdemat:io:ResourceLimit','JSON output exceeds the byte budget.');
        end
        envelope.arrays = encodeArrays(arrays);
        bytes = unicode2native(jsonencode(envelope), 'UTF-8');
        if numel(bytes) > limits.MaxFileBytes
            error('mdemat:io:ResourceLimit','JSON output exceeds the byte budget.');
        end
        stream = fopen(temporary, 'wb');
        if stream < 0, error('mdemat:io:WriteFailed', 'Cannot open temporary output.'); end
        streamGuard = onCleanup(@() closeStream(stream));
        written = fwrite(stream, bytes, 'uint8');
        closeStatus = fclose(stream);
        clear streamGuard
        if written ~= numel(bytes) || closeStatus ~= 0
            error('mdemat:io:WriteFailed', 'Incomplete record write.');
        end
    end
    fileInfo = dir(temporary);
    if fileInfo.bytes > limits.MaxFileBytes
        error('mdemat:io:ResourceLimit','Output file exceeds the byte budget.');
    end
    source = java.io.File(temporary);
    target = java.io.File(file);
    moveOptions = javaArray('java.nio.file.CopyOption', double(options.overwrite));
    if options.overwrite
        moveOptions(1) = javaMethod('valueOf', ...
            'java.nio.file.StandardCopyOption', 'REPLACE_EXISTING');
    end
    javaMethod('move', 'java.nio.file.Files', source.toPath(), target.toPath(), moveOptions);
catch exception
    if startsWith(exception.identifier,'mdemat:io:')
        rethrow(exception);
    end
    if ~options.overwrite && isfile(file)
        failure = MException('mdemat:io:FileExists', 'Refusing to overwrite %s.', file);
    else
        failure = MException('mdemat:io:WriteFailed', 'Record write failed: %s', file);
    end
    throwAsCaller(addCause(failure, exception));
end
end

function bytes = minimumJsonPayloadBytes(arrays)
% A lower bound suffices to refuse oversized payloads before base64 allocation.
bytes = 0;
names = fieldnames(arrays);
for index = 1:numel(names)
    value = arrays.(names{index});
    spec = primitiveSpec(class(value));
    count = numel(value);
    if issparse(value), count = nnz(value); end
    bytes = bytes+4*ceil(count*spec.Bytes/3)*(1+double(~isreal(value)));
end
end

function removeTemporary(file)
if isfile(file), delete(file); end
end

function closeStream(stream)
% fclose after a successful close is harmlessly ignored during cleanup.
if any(openedFiles == stream), fclose(stream); end
end
