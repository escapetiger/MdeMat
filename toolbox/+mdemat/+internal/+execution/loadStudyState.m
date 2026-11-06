function state = loadStudyState(filename, prototype)
%LOADSTUDYSTATE Bind callbacks from an explicitly built research model.
%   PROTOTYPE is supplied by the entrypoint's fixed constructor and initialized
%   model. The record never selects a class constructor or executable callback.
record = mdemat.tools.io.readRecord(filename);
if isstruct(record) && isscalar(record) && isfield(record,'format') && ...
        isequal(record.format,'mdemat.molCheckpoint')
    checkpoint = mdemat.internal.execution.readCheckpoint(filename);
    record = checkpoint.state;
end
prototypeRecord = prototype.snapshot(omitCallbacks=true);
record = mdemat.data.variable.validateStateRecord(record,prototypeRecord.kind);
callbacks = struct();
for index = 1:numel(record.callbackNames)
    name = record.callbackNames{index};
    if ~isfield(prototype.Coefs, name) || ~isa(prototype.Coefs.(name), 'function_handle')
        error('mdemat:study:MissingCallback', ...
            'The explicit research model does not supply callback %s.', name);
    end
    callbacks.(name) = prototype.Coefs.(name);
end
state = prototype.fromRecord(record, callbacks=callbacks);
end
