function [data, info] = readRecord(file, options)
%READRECORD Read a bounded v1 primitive record without saved-class construction.
%   [DATA,INFO] = readRecord(FILE,limits=STRUCT) reads MAT/JSON written by
%   writeRecord. Limits may override MaxFileBytes (2 GiB), MaxArrayBytes
%   (2 GiB-1), MaxTotalBytes (4 GiB), MaxElements (1e8), MaxNodes (1e6),
%   MaxDepth (128) and MaxDimensions (32). They are finite positive integers.
%   Input/parser allocation is bounded by file/manifest limits; this is not
%   complete memory isolation. The input must remain unchanged during reading.
%   Only whitelisted primitive MAT variables are loaded. Objects, constructors,
%   function handles and eval are not part of this format.
arguments
    file {mustBeTextScalar, mustBeFile}
    options.limits (1,1) struct = struct()
end
limits = recordLimits(options.limits);
file = char(file);
[~,~,extension] = fileparts(file);
try
    listing = dir(file);
    if ~isscalar(listing) || listing.bytes > limits.MaxFileBytes
        error('mdemat:io:ResourceLimit','Input file exceeds the byte budget.');
    end
    switch lower(extension)
        case '.mat'
            variables = whos('-file',file);
            if numel(variables) > limits.MaxNodes+1
                error('mdemat:io:ResourceLimit','Too many MAT primitive variables.');
            end
            names = {variables.name};
            descriptions = struct();
            for index = 1:numel(variables)
                item = variables(index);
                if ~strcmp(item.name,'manifest') && isempty(regexp(item.name,'^n\d+$','once'))
                    error('mdemat:io:InvalidRecord','Invalid MAT primitive variable name.');
                end
                spec = primitiveSpec(item.class);
                recordShape(item.size,limits,item.sparse);
                if (item.sparse && ~spec.Sparse) || (item.complex && ~spec.Complex)
                    error('mdemat:io:InvalidRecord','Invalid MAT primitive representation.');
                end
                if item.bytes > limits.MaxArrayBytes
                    error('mdemat:io:ResourceLimit','MAT variable exceeds the byte budget.');
                end
                descriptions.(item.name) = struct('Type',item.class,'Bytes',item.bytes);
            end
            if ~ismember('manifest',names) || ~strcmp(descriptions.manifest.Type,'uint8')
                error('mdemat:io:InvalidRecord','MAT record requires a uint8 manifest.');
            end
            manifest = load(file,'manifest');
            envelope = jsondecode(native2unicode(manifest.manifest(:).','UTF-8'));
            descriptions = rmfield(descriptions,'manifest');
        case '.json'
            envelope = jsondecode(fileread(file,Encoding='UTF-8'));
        otherwise
            error('mdemat:io:InvalidFormat','Record filename must end in .mat or .json.');
    end
    if ~isstruct(envelope) || ~isscalar(envelope) || ...
            ~all(isfield(envelope,{'format','schemaVersion','tree'})) || ...
            ~isequal(envelope.format,'mdemat.numericRecord') || ~isequal(envelope.schemaVersion,1)
        error('mdemat:io:UnsupportedSchema','Unsupported numeric record schema.');
    end
    expected = {'format','schemaVersion','tree'};
    producer = '';
    if isfield(envelope,'toolboxVersion')
        expected{end+1} = 'toolboxVersion';
        producer = envelope.toolboxVersion;
        if ~ischar(producer) || ~isrow(producer) || ...
                isempty(regexp(producer,'^\d+\.\d+\.\d+(\.\d+)?$','once'))
            error('mdemat:io:InvalidRecord','Invalid producer toolbox version.');
        end
    end
    if strcmpi(extension,'.json')
        expected{end+1} = 'arrays';
    end
    if ~isequal(sort(fieldnames(envelope)),sort(expected(:)))
        error('mdemat:io:InvalidRecord','Unexpected numeric envelope fields.');
    end
    if strcmpi(extension,'.json')
        descriptions = describeJsonArrays(envelope.arrays,limits);
    end
    references = validateRecordTree(envelope.tree,descriptions,limits);
    if strcmpi(extension,'.mat')
        arrays = struct();
        if ~isempty(references), arrays = load(file,references{:}); end
    else
        arrays = decodeArrays(envelope.arrays,references);
    end
    data = unpackRecord(envelope.tree,arrays);
catch exception
    if startsWith(exception.identifier,'mdemat:io:'), rethrow(exception); end
    failure = MException('mdemat:io:InvalidRecord','Cannot read numeric record: %s',file);
    throwAsCaller(addCause(failure,exception));
end
info = struct('Format','mdemat.numericRecord','SchemaVersion',1,'ToolboxVersion',producer);
end
