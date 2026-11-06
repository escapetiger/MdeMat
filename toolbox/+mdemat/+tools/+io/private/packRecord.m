function [tree, arrays] = packRecord(value, limits)
% Fixed data-only grammar; no user classes or callable values enter storage.
arrays = struct();
index = 0;
nodes = 0;
totalBytes = 0;
if nargin < 2, limits = recordLimits(); end
primitives = {'double','single','logical','char','int8','uint8', ...
    'int16','uint16','int32','uint32','int64','uint64'};
tree = encode(value,1);
    function node = encode(item,depth)
        nodes = nodes+1;
        if depth > limits.MaxDepth || nodes > limits.MaxNodes
            error('mdemat:io:ResourceLimit','Record node/depth budget exceeded.');
        end
        recordShape(size(item),limits,issparse(item));
        if ismember(class(item), primitives)
            if issparse(item) && ~ismember(class(item), {'double', 'logical'})
                error('mdemat:io:UnsupportedValue', ...
                    'Sparse records support double and logical arrays.');
            end
            index = index + 1;
            name = sprintf('n%06d', index);
            kind = 'array';
            if ischar(item), kind = 'char'; item = uint16(item); end
            description = whos('item');
            totalBytes = totalBytes+description.bytes;
            if description.bytes > limits.MaxArrayBytes || totalBytes > limits.MaxTotalBytes
                error('mdemat:io:ResourceLimit','Primitive byte budget exceeded.');
            end
            arrays.(name) = item;
            node = struct('kind', kind, 'variable', name);
        elseif isstring(item)
            if sum(strlength(item),'all','omitnan')*2+8*numel(item) > limits.MaxArrayBytes
                error('mdemat:io:ResourceLimit','String text exceeds the byte budget.');
            end
            missing = ismissing(item);
            item(missing) = "";
            node = struct('kind', 'string', 'shape', size(item), ...
                'text', {cellstr(item(:))}, 'missing', missing(:).');
        elseif isstruct(item)
            fields = fieldnames(item);
            if numel(fields)*numel(item) > limits.MaxNodes
                error('mdemat:io:ResourceLimit','Struct children exceed the node budget.');
            end
            values = cell(numel(fields)*numel(item), 1);
            cursor = 0;
            for element = 1:numel(item)
                for field = 1:numel(fields)
                    cursor = cursor + 1;
                    values{cursor} = encode(item(element).(fields{field}),depth+1);
                end
            end
            node = struct('kind', 'struct', 'shape', size(item), ...
                'fields', {fields}, 'items', {values});
        elseif iscell(item)
            if numel(item) > limits.MaxNodes
                error('mdemat:io:ResourceLimit','Cell children exceed the node budget.');
            end
            values = cell(numel(item), 1);
            for element = 1:numel(item), values{element} = encode(item{element},depth+1); end
            node = struct('kind', 'cell', 'shape', size(item), 'items', {values});
        else
            error('mdemat:io:UnsupportedValue', ...
                'Only numeric/logical arrays, text, structs and cells are supported; received %s.', ...
                class(item));
        end
    end
end
