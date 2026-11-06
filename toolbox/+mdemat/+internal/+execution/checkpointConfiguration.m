function description = checkpointConfiguration(config, options)
%CHECKPOINTCONFIGURATION Describe validated model values and callback paths.
%   CONFIG is a resolved plain value struct, never a Config or validator table.
%   modelKeys explicitly selects model parameters; the default excludes the
%   execution observer/output keys. Callback code is not stored or compared.
%   Nested callback paths are recorded and their value position becomes [].
%   Matching paths requires the caller to supply the configured model; it does
%   not prove that two independently supplied callback implementations agree.
arguments
    config (1,1) struct
    options.modelKeys cell = mdemat.internal.execution.checkpointModelKeys(fieldnames(config))
end
names = options.modelKeys(:);
if ~all(cellfun(@(name) ischar(name) && isvarname(name),names)) || ...
        numel(unique(names)) ~= numel(names) || ~all(isfield(config,names))
    error('mdemat:persistence:UnsupportedConfiguration','Select explicit existing model keys.');
end
values = struct();
callbacks = {};
for index = 1:numel(names)
    key = names{index};
    values.(key) = describe(config.(key),key);
end
description = struct('schemaVersion',2,'values',orderfields(values), ...
    'callbackNames',{sort(callbacks(:).')});
    function value = describe(value, path)
        if isa(value,'function_handle')
            callbacks{end+1} = path;
            value = [];
        elseif isnumeric(value) || islogical(value) || ischar(value) || isstring(value)
            return;
        elseif isstruct(value)
            fields = fieldnames(value);
            for element = 1:numel(value)
                prefix = path;
                if ~isscalar(value), prefix = sprintf('%s(%d)',path,element); end
                for field = 1:numel(fields)
                    fieldName = fields{field};
                    value(element).(fieldName) = describe( ...
                        value(element).(fieldName),[prefix,'.',fieldName]);
                end
            end
            value = orderfields(value);
        elseif iscell(value)
            for element = 1:numel(value)
                value{element} = describe(value{element},sprintf('%s{%d}',path,element));
            end
        else
            error('mdemat:persistence:UnsupportedConfiguration', ...
                'Model parameter %s requires an explicit value/callback description.',path);
        end
    end
end
