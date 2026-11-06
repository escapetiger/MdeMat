function descriptions = describeJsonArrays(encoded, limits)
% Validate the JSON primitive grammar and byte counts before base64 allocation.
if ~isstruct(encoded) || ~isscalar(encoded)
    error('mdemat:io:InvalidRecord','Invalid JSON array container.');
end
descriptions = struct();
names = fieldnames(encoded);
if numel(names) > limits.MaxNodes
    error('mdemat:io:ResourceLimit','Too many primitive arrays.');
end
for index = 1:numel(names)
    name = names{index};
    item = encoded.(name);
    if isempty(regexp(name,'^n\d+$','once')) || ~isstruct(item) || ~isscalar(item) || ...
            ~all(isfield(item,{'type','shape','sparse','complex','real'}))
        error('mdemat:io:InvalidRecord','Invalid JSON primitive description.');
    end
    spec = primitiveSpec(item.type);
    if ~flag(item.sparse) || ~flag(item.complex) || ...
            (item.sparse && ~spec.Sparse) || (item.complex && ~spec.Complex)
        error('mdemat:io:InvalidRecord','Invalid primitive sparsity/complex flags.');
    end
    expected = {'type','shape','sparse','complex','real'};
    if item.complex, expected{end+1} = 'imaginary'; end
    if item.sparse, expected = [expected,{'rows','columns'}]; end
    if ~isequal(sort(fieldnames(item)),sort(expected(:)))
        error('mdemat:io:InvalidRecord','Primitive fields do not match its flags.');
    end
    [shape,count] = recordShape(item.shape,limits,item.sparse);
    if item.sparse
        if numel(shape) ~= 2 || ~isnumeric(item.rows) || ~isnumeric(item.columns) || ...
                ~isreal(item.rows) || ~isreal(item.columns) || ...
                ~(isvector(item.rows) || isempty(item.rows)) || ...
                ~(isvector(item.columns) || isempty(item.columns)) || ...
                numel(item.rows) ~= numel(item.columns)
            error('mdemat:io:InvalidRecord','Invalid sparse index arrays.');
        end
        rows = item.rows(:);
        columns = item.columns(:);
        count = numel(rows);
        bytes = count*(spec.Bytes*(1+double(item.complex))+8)+(shape(2)+1)*8;
        if count > limits.MaxElements || bytes > limits.MaxArrayBytes
            error('mdemat:io:ResourceLimit','Sparse reconstruction exceeds the resource budget.');
        end
        if any(~isfinite(rows) | rows < 1 | ...
                rows ~= fix(rows) | rows > shape(1)) || ...
                any(~isfinite(columns) | columns < 1 | columns ~= fix(columns) | columns > shape(2)) || ...
                size(unique([rows,columns],'rows'),1) ~= count
            error('mdemat:io:InvalidRecord','Sparse indices must be unique and inside the matrix.');
        end
    else
        bytes = count*spec.Bytes*(1+double(item.complex));
    end
    if bytes > limits.MaxArrayBytes
        error('mdemat:io:ResourceLimit','Primitive array byte budget exceeded.');
    end
    byteCount(item.real,count*spec.Bytes);
    if item.complex, byteCount(item.imaginary,count*spec.Bytes); end
    descriptions.(name) = struct('Type',item.type,'Bytes',bytes);
end
end

function valid = flag(value)
valid = (islogical(value) || isnumeric(value)) && isscalar(value) && ...
    isreal(value) && isfinite(value) && (value == 0 || value == 1);
end

function byteCount(text, expected)
if ~ischar(text) || ~(isrow(text) || isempty(text)) || mod(numel(text),4) ~= 0 || ...
        (~isempty(text) && isempty(regexp(text,'^[A-Za-z0-9+/]*={0,2}$','once')))
    error('mdemat:io:InvalidRecord','Invalid base64 byte payload.');
end
padding = 0;
if ~isempty(text)
    padding = double(text(end) == '=');
    if numel(text) >= 2, padding = padding+double(text(end-1) == '='); end
end
if 3*(numel(text)/4)-padding ~= expected
    error('mdemat:io:InvalidRecord','Primitive byte count does not match its description.');
end
end
