function record = bindStateCallbacks(record, callbacks)
%BINDSTATECALLBACKS Bind exactly the callback names in a validated State record.
%   Saved metadata never supplies code or constructors. The caller supplies all
%   callbacks explicitly, and no numeric coefficient can be overwritten.
arguments
    record (1,1) struct
    callbacks (1,1) struct = struct()
end
names = record.callbackNames;
if ~isequal(sort(fieldnames(callbacks)),sort(names(:))) || ...
        ~all(structfun(@(value) isa(value,'function_handle'),callbacks)) || ...
        any(ismember(names(:),fieldnames(record.coefs)))
    error('mdemat:persistence:CallbacksRequired', ...
        'Supply exactly the recorded coefficient callbacks: %s.',strjoin(names,', '));
end
for index = 1:numel(names), record.coefs.(names{index}) = callbacks.(names{index}); end
end
