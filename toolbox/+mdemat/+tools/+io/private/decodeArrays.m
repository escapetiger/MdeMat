function output = decodeArrays(encoded, names)
if ~isstruct(encoded) || ~isscalar(encoded)
    error('mdemat:io:InvalidRecord', 'Invalid JSON array container.');
end
output = struct();
if nargin < 2, names = fieldnames(encoded); end
allowed = {'double', 'single', 'logical', 'int8', 'uint8', ...
    'int16', 'uint16', 'int32', 'uint32', 'int64', 'uint64'};
for i = 1:numel(names)
    item = encoded.(names{i});
    if isempty(regexp(names{i}, '^n\d+$', 'once')) || ~ismember(item.type, allowed)
        error('mdemat:io:InvalidRecord', 'Invalid JSON primitive type.');
    end
    shape = double(item.shape(:).');
    if numel(shape) < 2 || any(~isfinite(shape) | shape < 0 | shape ~= fix(shape))
        error('mdemat:io:InvalidRecord', 'Invalid JSON array dimensions.');
    end
    values = decodeBytes(item.real, item.type);
    if item.complex
        if ~ismember(item.type, {'double', 'single'})
            error('mdemat:io:InvalidRecord', 'Only floating arrays may be complex.');
        end
        imaginary = decodeBytes(item.imaginary, item.type);
        if numel(values) ~= numel(imaginary)
            error('mdemat:io:InvalidRecord','Real and imaginary payloads must have identical lengths.');
        end
        values = complex(values, imaginary);
    end
    if item.sparse
        if numel(shape) ~= 2 || ~ismember(item.type, {'double', 'logical'}) || ...
                numel(item.rows) ~= numel(values) || numel(item.columns) ~= numel(values)
            error('mdemat:io:InvalidRecord', 'Invalid sparse record representation.');
        end
        values = sparse(item.rows(:), item.columns(:), double(values), shape(1), shape(2));
        if strcmp(item.type, 'logical'), values = logical(values); end
    else
        if numel(values) ~= prod(shape)
            error('mdemat:io:InvalidRecord', 'Array byte count does not match its dimensions.');
        end
        values = reshape(values, shape);
    end
    output.(names{i}) = values;
end
end

function values = decodeBytes(text, type)
bytes = matlab.net.base64decode(char(text));
if strcmp(type, 'logical')
    if any(bytes ~= 0 & bytes ~= 1)
        error('mdemat:io:InvalidRecord','Logical bytes must be zero or one.');
    end
    values = logical(bytes(:));
    return;
end
values = typecast(bytes(:), type);
[~, ~, endian] = computer;
if endian == 'B', values = swapbytes(values); end
values = values(:);
end
