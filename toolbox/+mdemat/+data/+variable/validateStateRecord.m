function record = validateStateRecord(record, kind)
%VALIDATESTATERECORD Check the shared v1 State value schema without binding code.
%   Kind/space/application parameters are selected by explicit fromRecord calls.
arguments
    record (1,1) struct
    kind {mustBeTextScalar}
end
expected = {'format','schemaVersion','kind','dofs','coefs','history', ...
    'callbackNames','spaces','parameters'};
if ~isequal(sort(fieldnames(record)),sort(expected(:))) || ...
        ~isequal(record.format,'mdemat.state') || ~isequal(record.schemaVersion,1) || ...
        ~isequal(record.kind,char(kind)) || ~numericStruct(record.dofs) || ...
        ~numericStruct(record.coefs) || ~numericStruct(record.history) || ...
        ~isstruct(record.spaces) || ~isscalar(record.spaces) || ...
        ~isstruct(record.parameters) || ~isscalar(record.parameters)
    error('mdemat:persistence:InvalidState','State schema/kind/data does not match this loader.');
end
names = record.callbackNames;
if ~iscell(names) || ~all(cellfun(@(name) ischar(name) && isvarname(name),names),'all') || ...
        numel(unique(names)) ~= numel(names) || any(ismember(names(:),fieldnames(record.coefs)))
    error('mdemat:persistence:InvalidState','Invalid required callback names.');
end
end

function valid = numericStruct(value)
valid = isstruct(value) && isscalar(value);
if valid, valid = all(structfun(@isnumeric,value)); end
end
