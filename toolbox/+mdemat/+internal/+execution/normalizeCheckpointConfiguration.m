function description = normalizeCheckpointConfiguration(description, options)
%NORMALIZECHECKPOINTCONFIGURATION Project v1/v2 saved identities before comparison.
%   Old v1 observer fields are removed in memory. Original files are unchanged.
%   New descriptors and v1 descriptors normalize to the same schemaVersion=2
%   representation; an explicit modelKeys subset follows the model owner.
arguments
    description (1,1) struct
    options.modelKeys cell = {}
end
expected = {'values','callbackNames'};
if isfield(description,'schemaVersion')
    expected{end+1} = 'schemaVersion';
    if ~isequal(description.schemaVersion,2)
        error('mdemat:persistence:InvalidCheckpoint','Unsupported configuration description.');
    end
end
if ~isequal(sort(fieldnames(description)),sort(expected(:))) || ...
        ~isstruct(description.values) || ~isscalar(description.values) || ...
        ~iscell(description.callbackNames) || ...
        ~all(cellfun(@(name) ischar(name) && isrow(name),description.callbackNames),'all') || ...
        numel(unique(description.callbackNames)) ~= numel(description.callbackNames)
    error('mdemat:persistence:InvalidCheckpoint','Invalid configuration description.');
end
values = description.values;
callbacks = description.callbackNames(:).';
% Old v1 omitted top-level callable fields entirely.
for index = 1:numel(callbacks)
    name = callbacks{index};
    if isvarname(name) && ~isfield(values,name), values.(name) = []; end
end
if isempty(options.modelKeys)
    keys = mdemat.internal.execution.checkpointModelKeys(fieldnames(values));
else
    keys = options.modelKeys(:);
    if ~all(cellfun(@(name) ischar(name) && isvarname(name),keys)) || ...
            ~all(isfield(values,keys))
        error('mdemat:persistence:InvalidCheckpoint','Saved model parameters do not match selected keys.');
    end
end
values = rmfield(values,setdiff(fieldnames(values),keys));
keep = false(size(callbacks));
for index = 1:numel(callbacks)
    root = regexp(callbacks{index},'^[A-Za-z]\w*','match','once');
    if isempty(regexp(callbacks{index}, ...
            '^[A-Za-z]\w*(\.\w+|\([1-9]\d*\)|\{[1-9]\d*\})*$','once'))
        error('mdemat:persistence:InvalidCheckpoint','Invalid saved model callback path.');
    end
    keep(index) = ismember(root,keys);
end
description = struct('schemaVersion',2,'values',orderfields(values), ...
    'callbackNames',{sort(callbacks(keep))});
end
