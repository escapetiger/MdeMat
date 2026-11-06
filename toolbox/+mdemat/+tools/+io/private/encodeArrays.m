function output = encodeArrays(arrays)
output = struct();
names = fieldnames(arrays);
for i = 1:numel(names)
    data = arrays.(names{i});
    values = data(:);
    item = struct('type', class(data), 'shape', size(data), ...
        'sparse', issparse(data), 'complex', ~isreal(data));
    if issparse(data)
        [rows, columns, values] = find(data);
        item.rows = rows(:).';
        item.columns = columns(:).';
        % find returns double values for sparse logical input.
        if islogical(data), values = logical(values); end
    end
    if isreal(data)
        % real(logical(...)) promotes to double; retain the primitive class.
        item.real = encodeBytes(values);
    else
        item.real = encodeBytes(real(values));
        item.imaginary = encodeBytes(imag(values));
    end
    output.(names{i}) = item;
end
end

function text = encodeBytes(values)
values = full(values);
if islogical(values), values = uint8(values); end
[~, ~, endian] = computer;
if endian == 'B', values = swapbytes(values); end
text = matlab.net.base64encode(typecast(values(:), 'uint8'));
end
