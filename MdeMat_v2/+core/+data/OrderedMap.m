classdef OrderedMap < handle
    % ORDEREDMAP Key-value map that preserves insertion order.
    %
    %   OrderedMap provides containers.Map-like interface but preserves
    %   insertion order of elements. Supports basic operations (get, set,
    %   remove, isKey) and adds methods for order manipulation.
    %
    % Examples:
    %   % Create and populate ordered map
    %   map = core.data.OrderedMap();
    %   map('first') = 1;
    %   map('second') = 2;
    %   map('third') = 3;
    %
    %   % Iterate through keys in insertion order
    %   keys = map.keys();
    %   % keys will be {'first', 'second', 'third'} in that order
    %
    %   % Check if key exists
    %   exists = map.isKey('third');  % exists = true

    properties
        keys % Ordered keys
    end

    properties (Dependent)
        values % Ordered values
        items % Ordered items
        keyType % Key type
        valueType % Value type
    end

    properties (Access = protected)
        map % Internal containers.Map for fast key-value lookups
    end

    methods
        function obj = OrderedMap(varargin)
            % ORDEREDMAP Constructor for OrderedMap class.
            %
            %   obj = OrderedMap() creates empty OrderedMap.
            %
            %   obj = OrderedMap(keys, values) creates OrderedMap from cell
            %   arrays.
            %
            %   obj = OrderedMap('key1', value1, 'key2', value2, ...)
            %   creates OrderedMap from key-value pairs.
            %
            % Inputs:
            %   varargin - Input arguments
            %
            % Outputs:
            %   obj - Constructed OrderedMap object

            if nargin == 0
                obj.constructEmptyMap();
                obj.keys = {};
                obj.map = containers.Map();
            elseif nargin == 2
                keys = varargin{1};
                values = varargin{2};
                if iscell(keys) && iscell(values)
                    obj = obj.constructFromArrays(keys, values);
                end
            elseif mod(nargin, 2) == 0
                nPairs = nargin / 2;
                keys = cell(1, nPairs);
                values = cell(1, nPairs);
                for iPair = 1:nPairs
                    keys{iPair} = varargin{2*iPair-1};
                    values{iPair} = varargin{2*iPair};
                end
                obj = obj.constructFromArrays(keys, values);

            else
                core.except.assert(false, 'InvalidInput', ...
                    ['Invalid number of arguments. Expected 0, 2, or an even ', ...
                    'number of arguments for key-value pairs. Got %d arguments.'], ...
                    nargin);
            end
        end

        function result = get.values(obj)
            % GET.VALUES Return cell array of values in insertion order.
            %
            %   result = values(obj) returns all values in OrderedMap in
            %   original insertion order.
            %
            % Inputs:
            %   obj - The OrderedMap object
            %
            % Outputs:
            %   result - Cell array of values in insertion order

            result = cell(size(obj.keys));
            for i = 1:length(obj.keys)
                result{i} = obj.map(obj.keys{i});
            end
        end

        function result = get.items(obj)
            % GET.ITEMS Return matrix of key-value pairs for iteration.
            %
            %   result = items(obj) returns cell matrix where first column
            %   contains keys and second column contains corresponding
            %   values.
            %
            % Inputs:
            %   obj - The OrderedMap object
            %
            % Outputs:
            %   result - n×2 cell matrix where:
            %            result{i,1} is the key of the i-th element
            %            result{i,2} is the value of the i-th element

            keyList = obj.keys();
            result = cell(length(keyList), 2);
            result(:, 1) = keyList(:);

            for i = 1:length(keyList)
                result{i, 2} = obj.map(keyList{i});
            end
        end

        function result = get.keyType(obj)
            % GET.KEYTYPE Return the key type of the map.
            %
            %   result = keyType(obj) returns key type string.
            %
            % Inputs:
            %   obj - The OrderedMap object
            %
            % Outputs:
            %   result - Key type string

            result = obj.map.KeyType;
        end

        function result = get.valueType(obj)
            % GET.VALUETYPE Return the value type of the map.
            %
            %   result = valueType(obj) returns value type string.
            %
            % Inputs:
            %   obj - The OrderedMap object
            %
            % Outputs:
            %   result - Value type string

            result = obj.map.ValueType;
        end

        function result = length(obj)
            % LENGTH Return the number of key-value pairs in the map.
            %
            %   result = length(obj) returns number of key-value pairs.
            %
            % Inputs:
            %   obj - The OrderedMap object
            %
            % Outputs:
            %   result - Number of key-value pairs in the map

            result = length(obj.map.keys);
        end

        function result = count(obj)
            % COUNT Return the number of key-value pairs in the map.
            %
            %   result = count(obj) returns number of key-value pairs.
            %
            % Inputs:
            %   obj - The OrderedMap object
            %
            % Outputs:
            %   result - Number of key-value pairs in the map

            result = length(obj);
        end

        function result = size(obj, dim)
            % SIZE Return the size of the map.
            %
            %   result = size(obj) returns size of the map.
            %
            %   result = size(obj, dim) returns size of specified dimension.
            %
            % Inputs:
            %   obj - The OrderedMap object
            %   dim - (Optional) Dimension to query
            %
            % Outputs:
            %   result - Size of the map ([1, n] where n is the length)

            s = [1, length(obj.keys)];
            if nargin < 2
                result = s;
            else
                if dim <= 2
                    result = s(dim);
                else
                    result = 1;
                end
            end
        end

        function result = isKey(obj, key)
            % ISKEY Check if a key exists in the map.
            %
            %   result = isKey(obj, key) determines whether specified @a
            %   key exists in the OrderedMap.
            %
            % Inputs:
            %   obj - The OrderedMap object
            %   key - Key to check for existence
            %
            % Outputs:
            %   result - True if the key exists, false otherwise

            result = obj.map.isKey(key);
        end

        function varargout = subsref(obj, s)
            % SUBSREF Overloaded method for subscript referencing.
            %
            %   Enables map access using index notation: map('key') or
            %   map.methodName().

            switch s(1).type
                case '()'
                    key = s(1).subs{1};
                    core.except.assert(obj.isKey(key), 'KeyNotFound', ...
                        'Key "%s" not found in OrderedMap.', ...
                        core.data.toString(key));
                    value = obj.map(key);

                    if length(s) > 1
                        value = subsref(value, s(2:end));
                    end

                    varargout{1} = value;

                case '.'
                    if length(s) == 1
                        [varargout{1:nargout}] = builtin('subsref', obj, s);
                    else
                        if nargout == 0
                            builtin('subsref', obj, s);
                        else
                            [varargout{1:nargout}] = builtin('subsref', obj, s);
                        end
                    end

                otherwise
                    core.except.assert(0, 'invalidInput', ...
                        'Unsupported subscript reference type.');
            end
        end

        function obj = subsasgn(obj, s, value)
            % SUBSASGN Overloaded method for subscript assignment.

            switch s(1).type
                case '()'
                    key = s(1).subs{1};
                    obj.setItem(key, value);
                case '.'
                    obj = builtin('subsasgn', obj, s, value);
                otherwise
                    core.except.assert(0, 'invalidInput', ...
                        'Unsupported subscript assignment type.');
            end
        end

        function obj = setItem(obj, key, value)
            % SETITEM Set a key-value pair in the map.
            %
            %   obj = setItem(obj, key, value) sets the @a value for
            %   specified @a key. If key exists, value is updated without
            %   changing insertion order.
            %
            % Inputs:
            %   obj - The OrderedMap object
            %   key - Key to set (any valid containers.Map key type)
            %   value - Value to associate with the key
            %
            % Outputs:
            %   obj - The OrderedMap object

            isNewKey = ~obj.isKey(key);
            obj.map(key) = value;
            if isNewKey, obj.keys{end+1} = key; end
        end

        function obj = remove(obj, key)
            % REMOVE Remove a key-value pair from the map.
            %
            %   obj = remove(obj, key) removes specified @a key and its
            %   corresponding value from OrderedMap.
            %
            % Inputs:
            %   obj - The OrderedMap object
            %   key - Key to remove
            %
            % Outputs:
            %   obj - The OrderedMap object

            core.except.assert(obj.isKey(key), 'KeyNotFound', ...
                'Cannot remove key "%s" because it does not exist in the map.', ...
                core.data.toString(key));

            remove(obj.map, key);

            for i = 1:length(obj.keys)
                if isequal(obj.keys{i}, key)
                    obj.keys(i) = [];
                    break;
                end
            end
        end

        function obj = clear(obj)
            % CLEAR Remove all key-value pairs from the map.
            %
            %   obj = clear(obj) removes all key-value pairs from the map.
            %
            % Inputs:
            %   obj - The OrderedMap object
            %
            % Outputs:
            %   obj - The OrderedMap object

            obj.map = containers.Map();
            obj.keys = {};
        end

        function obj = move(obj, key, newPosition)
            % MOVE Move a key to a new position in the order.
            %
            %   obj = move(obj, key, newPosition) changes position of existing
            %   @a key to @a newPosition in the ordering.
            %
            % Inputs:
            %   obj - The OrderedMap object
            %   key - Key to move
            %   newPosition - New position index (1-based)
            %
            % Outputs:
            %   obj - The OrderedMap object

            core.except.assert(obj.isKey(key), 'KeyNotFound', ...
                'Cannot move key "%s" because it does not exist in the map.', ...
                core.data.toString(key));

            n = length(obj.keys);
            core.except.assert(newPosition >= 1 && newPosition <= n, ...
                'InvalidPosition', 'Position must be between 1 and %d, got %d.', ...
                n, newPosition);

            currentPosition = 0;
            for i = 1:n
                if isequal(obj.keys{i}, key)
                    currentPosition = i;
                    break;
                end
            end

            if currentPosition == newPosition
                return;
            end

            keyToMove = obj.keys{currentPosition};
            obj.keys(currentPosition) = [];

            if newPosition == 1
                obj.keys = [{keyToMove}, obj.keys];
            elseif newPosition > length(obj.keys)
                obj.keys{end+1} = keyToMove;
            else
                obj.keys = [obj.keys(1:newPosition-1), ...
                    {keyToMove}, ...
                    obj.keys(newPosition:end)];
            end
        end

        function obj = moveToFront(obj, key)
            % MOVETOFRONT Move a key to the beginning of the order.
            %
            %   obj = moveToFront(obj, key) moves specified @a key to first
            %   position in the ordering.
            %
            % Inputs:
            %   obj - The OrderedMap object
            %   key - Key to move to the front
            %
            % Outputs:
            %   obj - The OrderedMap object

            if obj.isKey(key), obj.move(key, 1); end
        end

        function obj = moveToEnd(obj, key)
            % MOVETOEND Move a key to the end of the order.
            %
            %   obj = moveToEnd(obj, key) moves specified @a key to last
            %   position in the ordering.
            %
            % Inputs:
            %   obj - The OrderedMap object
            %   key - Key to move to the end
            %
            % Outputs:
            %   obj - The OrderedMap object

            if obj.isKey(key), obj.move(key, length(obj.keys)); end
        end

        function tbl = toTable(obj)
            % TOTABLE Convert the OrderedMap to a table for easy viewing.
            %
            %   tbl = toTable(obj) creates MATLAB table with two columns:
            %   'Key' and 'Value', ordered according to insertion order.
            %
            % Inputs:
            %   obj - The OrderedMap object
            %
            % Outputs:
            %   tbl - MATLAB table with keys and values

            if isempty(obj.keys)
                tbl = table('Size', [0, 2], ...
                    'VariableTypes', {'cell', 'cell'}, ...
                    'VariableNames', {'Key', 'Value'});
                return;
            end

            tbl = table(obj.keys(:), obj.values(:), ...
                'VariableNames', {'Key', 'Value'});
        end

        function obj = disp(obj)
            % DISP Display the OrderedMap contents.
            %
            %   obj = disp(obj) customizes display output for OrderedMap.
            %
            % Inputs:
            %   obj - The OrderedMap object
            %
            % Outputs:
            %   obj - The OrderedMap object

            n = length(obj);
            if n == 0
                disp('  OrderedMap with 0 key/value pairs');
                return;
            end

            fprintf('  OrderedMap with %d key/value pairs:\n\n', n);

            keyList = obj.keys();

            maxKeyLen = 0;
            for i = 1:length(keyList)
                key = keyList{i};
                if isnumeric(key) || islogical(key)
                    keyStr = num2str(key);
                elseif ischar(key)
                    keyStr = key;
                else
                    keyStr = class(key);
                end
                maxKeyLen = max(maxKeyLen, length(keyStr));
            end

            fmt = ['    %-', num2str(maxKeyLen), 's: %s\n'];

            for i = 1:length(keyList)
                key = keyList{i};
                value = obj.map(key);

                if isnumeric(key) || islogical(key)
                    keyStr = num2str(key);
                elseif ischar(key)
                    keyStr = key;
                else
                    keyStr = class(key);
                end

                if isnumeric(value) && isscalar(value)
                    valueStr = num2str(value);
                elseif islogical(value) && isscalar(value)
                    if value
                        valueStr = 'true';
                    else
                        valueStr = 'false';
                    end
                elseif ischar(value)
                    if length(value) > 50
                        valueStr = [value(1:47), '...'];
                    else
                        valueStr = value;
                    end
                else
                    if isnumeric(value) || islogical(value) || iscell(value)
                        sizeStr = sprintf('%dx', size(value));
                        sizeStr = sizeStr(1:end-1); % Remove trailing 'x'
                        valueStr = [class(value), ' [', sizeStr, ']'];
                    end
                end

                fprintf(fmt, keyStr, valueStr);
            end
        end
    end

    methods (Access = private)
        function obj = constructEmptyMap(obj)
            % CONSTRUCTEMPTYMAP Create an empty OrderedMap.

            obj.keys = {};
            obj.map = containers.Map();
        end

        function obj = constructFromArrays(obj, keys, values)
            % CONSTRUCTFROMARRAYS Create OrderedMap from cell arrays.

            core.except.assert(length(keys) == length(values), ...
                'LengthMismatch', ...
                'Keys and values must have the same length.');

            if isempty(keys)
                obj = obj.constructEmptyMap();
                return;
            end
            obj.keys = keys(:).';
            obj.map = containers.Map(keys, values);
        end
    end
end