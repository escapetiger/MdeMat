function value = unpackRecord(node, arrays)
if ~isstruct(node) || ~isscalar(node) || ~isfield(node, 'kind')
    error('mdemat:io:InvalidRecord', 'Invalid record node.');
end
switch node.kind
    case {'array', 'char'}
        if ~isfield(node, 'variable') || ...
                isempty(regexp(node.variable, '^n\d+$', 'once')) || ...
                ~isfield(arrays, node.variable)
            error('mdemat:io:InvalidRecord', 'Missing primitive record array.');
        end
        value = arrays.(node.variable);
        if strcmp(node.kind, 'char')
            if ~isa(value, 'uint16')
                error('mdemat:io:InvalidRecord', 'Character data must be UTF-16 numeric units.');
            end
            value = char(value);
        end
    case 'string'
        shape = validShape(node.shape);
        text = textItems(node.text);
        missingMask = logical(node.missing(:));
        if numel(text) ~= prod(shape) || numel(missingMask) ~= prod(shape)
            error('mdemat:io:InvalidRecord', 'String shape does not match data.');
        end
        value = string(text);
        value(missingMask) = missing;
        value = reshape(value, shape);
    case {'cell', 'struct'}
        shape = validShape(node.shape);
        items = itemCells(node.items);
        if strcmp(node.kind, 'cell')
            if numel(items) ~= prod(shape)
                error('mdemat:io:InvalidRecord', 'Cell shape does not match data.');
            end
            value = cell(shape);
            for i = 1:numel(items), value{i} = unpackRecord(items{i}, arrays); end
        else
            fields = textItems(node.fields);
            if numel(unique(fields)) ~= numel(fields) || ...
                    ~all(cellfun(@isvarname, fields)) || numel(items) ~= prod(shape)*numel(fields)
                error('mdemat:io:InvalidRecord', 'Invalid struct fields or data shape.');
            end
            template = cell2struct(cell(numel(fields), 1), fields, 1);
            value = repmat(template, shape);
            cursor = 0;
            for i = 1:numel(value)
                for field = 1:numel(fields)
                    cursor = cursor + 1;
                    value(i).(fields{field}) = unpackRecord(items{cursor}, arrays);
                end
            end
        end
    otherwise
        error('mdemat:io:InvalidRecord', 'Unknown record node kind.');
end
end

function shape = validShape(shape)
shape = double(shape(:).');
if numel(shape) < 2 || any(~isfinite(shape) | shape < 0 | shape ~= fix(shape))
    error('mdemat:io:InvalidRecord', 'Invalid record dimensions.');
end
end

function items = itemCells(items)
if isempty(items), items = cell(0, 1);
elseif isstruct(items), items = num2cell(items(:));
elseif iscell(items), items = items(:);
else, error('mdemat:io:InvalidRecord', 'Invalid record children.');
end
end

function values = textItems(values)
if isempty(values), values = cell(0, 1);
elseif ischar(values), values = {values};
elseif isstring(values), values = cellstr(values(:));
elseif iscell(values) && all(cellfun(@ischar, values)), values = values(:);
else, error('mdemat:io:InvalidRecord', 'Invalid record text.');
end
end

