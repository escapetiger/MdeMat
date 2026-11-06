function identity = versionInfo()
%VERSIONINFO Read trusted installed identity and supported persistence schemas.
%   Missing, malformed or unsupported metadata raises InvalidMetadata with
%   its original cause. Metadata is read afresh so development edits remain
%   visible; unknown schemas are never silently accepted.
file = fullfile(mdemat.tools.environment.toolboxRoot(), 'resources', 'toolbox.json');
try
    metadata = jsondecode(fileread(file));
catch exception
    failure = MException('mdemat:environment:InvalidMetadata', ...
        'Cannot read installed toolbox metadata: %s', file);
    throwAsCaller(addCause(failure, exception));
end
expected = {'name', 'version', 'identifier', 'minimumRelease', 'author', 'summary', 'schemas'};
valid = isstruct(metadata) && isscalar(metadata) && ...
    isequal(sort(fieldnames(metadata)), sort(expected(:)));
if valid
    textFields = expected(1:6);
    valid = all(cellfun(@(name) ischar(metadata.(name)) && ...
        isrow(metadata.(name)) && ~isempty(metadata.(name)), textFields));
end
if valid
    valid = ~isempty(regexp(metadata.version, '^\d+\.\d+\.\d+(\.\d+)?$', 'once')) && ...
        ~isempty(regexp(metadata.identifier, ...
        '^[0-9a-fA-F]{8}(-[0-9a-fA-F]{4}){3}-[0-9a-fA-F]{12}$', 'once')) && ...
        ~isempty(regexp(metadata.minimumRelease, '^R\d{4}[ab]$', 'once')) && ...
        isequal(metadata.schemas, struct('numericRecord', 1, 'runRecord', 1));
end
if ~valid
    error('mdemat:environment:InvalidMetadata', 'Unsupported toolbox identity/schema metadata.');
end
identity = struct('Name', metadata.name, 'ToolboxVersion', metadata.version, ...
    'Identifier', metadata.identifier, 'NumericRecordSchemaVersion', metadata.schemas.numericRecord, ...
    'RunRecordSchemaVersion', metadata.schemas.runRecord);
end
