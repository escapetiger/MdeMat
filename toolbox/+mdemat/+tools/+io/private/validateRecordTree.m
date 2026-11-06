function references = validateRecordTree(tree, descriptors, limits)
% Validate the complete description before reconstructing any container.
references = {};
nodes = 0;
elements = 0;
bytes = 0;
visit(tree,1);
    function visit(node,depth)
        nodes = nodes+1;
        if nodes > limits.MaxNodes || depth > limits.MaxDepth
            error('mdemat:io:ResourceLimit','Record node/depth budget exceeded.');
        end
        if ~isstruct(node) || ~isscalar(node) || ~isfield(node,'kind') || ...
                ~ischar(node.kind) || ~isrow(node.kind)
            error('mdemat:io:InvalidRecord','Invalid record node.');
        end
        switch node.kind
            case {'array','char'}
                exactFields(node,{'kind','variable'});
                name = node.variable;
                if ~ischar(name) || ~isrow(name) || isempty(regexp(name,'^n\d+$','once')) || ...
                        ~isfield(descriptors,name)
                    error('mdemat:io:InvalidRecord','Missing primitive array reference.');
                end
                description = descriptors.(name);
                if strcmp(node.kind,'char') && ~strcmp(description.Type,'uint16')
                    error('mdemat:io:InvalidRecord','Characters require uint16 UTF-16 units.');
                end
                if ~ismember(name,references)
                    references{end+1} = name;
                    bytes = bytes+description.Bytes;
                    if bytes > limits.MaxTotalBytes
                        error('mdemat:io:ResourceLimit','Decoded byte budget exceeded.');
                    end
                end
            case 'string'
                exactFields(node,{'kind','shape','text','missing'});
                [~,count] = recordShape(node.shape,limits);
                text = textItems(node.text);
                mask = node.missing;
                if ~(islogical(mask) || isnumeric(mask)) || ~isreal(mask) || ...
                        any(~isfinite(mask(:)) | (mask(:) ~= 0 & mask(:) ~= 1)) || ...
                        numel(mask) ~= count || numel(text) ~= count
                    error('mdemat:io:InvalidRecord','String mask/text does not match its shape.');
                end
                elements = elements+count;
                textBytes = 2*sum(cellfun(@numel,text))+8*count;
                bytes = bytes+textBytes;
                if textBytes > limits.MaxArrayBytes || bytes > limits.MaxTotalBytes
                    error('mdemat:io:ResourceLimit','String reconstruction byte budget exceeded.');
                end
            case {'cell','struct'}
                expected = {'kind','shape','items'};
                if strcmp(node.kind,'struct'), expected{end+1} = 'fields'; end
                exactFields(node,expected);
                [~,count] = recordShape(node.shape,limits);
                children = itemCells(node.items);
                if strcmp(node.kind,'struct')
                    fields = textItems(node.fields);
                    if numel(unique(fields)) ~= numel(fields) || ~all(cellfun(@isvarname,fields))
                        error('mdemat:io:InvalidRecord','Invalid struct fields.');
                    end
                    expectedCount = count*numel(fields);
                else
                    expectedCount = count;
                end
                if expectedCount > limits.MaxNodes
                    error('mdemat:io:ResourceLimit','Container child budget exceeded.');
                end
                if numel(children) ~= expectedCount
                    error('mdemat:io:InvalidRecord','Container children do not match its shape.');
                end
                elements = elements+count;
                if elements > limits.MaxElements
                    error('mdemat:io:ResourceLimit','Container element budget exceeded.');
                end
                for child = 1:numel(children), visit(children{child},depth+1); end
            otherwise
                error('mdemat:io:InvalidRecord','Unknown record node kind.');
        end
        if elements > limits.MaxElements
            error('mdemat:io:ResourceLimit','Container element budget exceeded.');
        end
    end
end

function exactFields(value, expected)
if ~isequal(sort(fieldnames(value)),sort(expected(:)))
    error('mdemat:io:InvalidRecord','Record node fields do not match its kind.');
end
end

function values = textItems(values)
if isempty(values), values = cell(0,1);
elseif ischar(values) && isrow(values), values = {values};
elseif iscell(values) && all(cellfun(@(value) ischar(value) && ...
        (isrow(value) || isempty(value)),values),'all'), values = values(:);
else, error('mdemat:io:InvalidRecord','Invalid record text list.');
end
end

function values = itemCells(values)
if isempty(values), values = cell(0,1);
elseif isstruct(values), values = num2cell(values(:));
elseif iscell(values), values = values(:);
else, error('mdemat:io:InvalidRecord','Invalid record children.');
end
end
