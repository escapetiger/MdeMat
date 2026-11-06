classdef CompiledFunction < mdemat.data.function.Function
    % COMPILEDFUNCTION Base class for all compiled functions.
    %
    %   CompiledFunction provides a bridge between symbolic function
    %   definitions and efficient numerical evaluation by automatically
    %   loading precompiled function handles and metadata from the
    %   bundled numeric resource directory. Symbolic development metadata is
    %   never loaded by this runtime class.
    %
    %   The class supports automatic metadata discovery, caching, and
    %   efficient function evaluation using precompiled function handles
    %   stored in MAT files.
    %   The bundled spherical-harmonic resource contains degree 0 through 5
    %   (36 real modes); larger requests fail with CapacityExceeded.
    %
    % See Also:
    %   mdemat.data.function.Function

    properties (Constant, Access = protected)
        MetaPath = 'resources/basis' % Relative to the distributable toolbox root
        MetaExt = '.mat' % File extension for metadata files
    end

    properties
        Metadata % Loaded metadata structure containing function handles
    end

    methods
        function obj = CompiledFunction(nDims, nCodims, metaName)
            % COMPILEDFUNCTION Constructor for the CompiledFunction class.
            %
            %   obj = CompiledFunction(nDims, nCodims, metaName) creates a
            %   CompiledFunction with specified input dimension @a nDims,
            %   output dimension @a nCodims, and metadata file @a metaName.

            arguments
                nDims(1, 1) {mustBeInteger}
                nCodims(1, 1) {mustBeInteger}
                metaName{mustBeTextScalar}
            end

            obj@mdemat.data.function.Function(nDims=nDims, nCodims=nCodims);
            obj.load(metaName);
        end

        function obj = load(obj, fileName)
            % LOAD Load compiled metadata from specified file.
            %
            %   obj = load(obj, fileName) loads precompiled function
            %   metadata from the specified file and updates the object's
            %   function dimensions based on the loaded metadata.

            arguments
                obj mdemat.common.approximation.CompiledFunction
                fileName {mustBeTextScalar}
            end

            name = char(fileName);
            if endsWith(name, obj.MetaExt)
                name = name(1:end-length(obj.MetaExt));
            end
            catalog = {'chebyshev', 'monic_unit_legendre', 'normal_fourier', ...
                'normal_hermite', 'normal_legendre', 'normal_spherical_harmonic'};
            isInterpolation = ~isempty(regexp(name, ...
                '^unit_gauss_(legendre|lobatto)_lagrange_[1-6]$', 'once'));
            assert(ismember(name, catalog) || isInterpolation, ...
                'MdeMat:BasisResource:InvalidName', ...
                'Only declared bundled basis names are accepted.');

            resourceRoot = mdemat.tools.environment.toolboxRoot();
            fullPath = fullfile(resourceRoot, obj.MetaPath, [name, obj.MetaExt]);
            assert(isfile(fullPath), 'MdeMat:BasisResource:FileNotFound', ...
                'Numeric basis resource not found: %s.', fullPath);
            data = load(fullPath, 'metadata');
            assert(isfield(data, 'metadata'), 'MdeMat:BasisResource:InvalidSchema', ...
                'Numeric resource must contain metadata.');
            metadata = data.metadata;
            metadata.handles = obj.compileResource(metadata, obj.NDims, obj.NCodims);
            assert(strcmp(metadata.basisName, name), 'MdeMat:BasisResource:InvalidSchema', ...
                'The resource basis name does not match its catalog entry.');
            metadata.functionText = metadata.functionText(1:obj.NCodims, :);
            metadata.nCodims = obj.NCodims;
            metadata.resourceFile = fullPath;
            obj.Metadata = metadata;
        end
    end

    methods (Static, Hidden)
        function handles = compileResource(metadata, nDims, nCodims)
            % Internal schema boundary shared with the development exporter.
            % Text is accepted only in the restricted numeric grammar below.
            fields = {'schemaVersion'; 'basisName'; 'nDims'; 'nCodims'; ...
                'derivativeOrders'; 'functionText'};
            valid = isstruct(metadata) && isscalar(metadata) && ...
                isequal(sort(fieldnames(metadata)), sort(fields));
            assert(valid, 'MdeMat:BasisResource:InvalidSchema', ...
                'Resource fields do not match the numeric schema.');
            valid = isequal(metadata.schemaVersion, 1) && ...
                ischar(metadata.basisName) && isrow(metadata.basisName) && ...
                isa(metadata.nDims, 'double') && isscalar(metadata.nDims) && ...
                metadata.nDims == nDims && ismember(nDims, [1, 2]) && ...
                isa(metadata.nCodims, 'double') && isscalar(metadata.nCodims) && ...
                isfinite(metadata.nCodims) && metadata.nCodims >= 0 && ...
                metadata.nCodims == fix(metadata.nCodims) && ...
                iscell(metadata.functionText) && ismatrix(metadata.functionText) && ...
                size(metadata.functionText, 1) == metadata.nCodims && ...
                size(metadata.functionText, 2) >= 1 && ...
                isequal(metadata.derivativeOrders, 0:size(metadata.functionText, 2)-1) && ...
                (nDims == 1 || size(metadata.functionText, 2) == 1);
            assert(valid, 'MdeMat:BasisResource:InvalidSchema', ...
                'Numeric resource dimensions or derivative orders are invalid.');
            assert(nCodims <= metadata.nCodims, ...
                'MdeMat:BasisResource:CapacityExceeded', ...
                'Basis %s contains at most %d modes; %d were requested. Regenerate the resource.', ...
                metadata.basisName, metadata.nCodims, nCodims);
            handles = cell(nCodims, size(metadata.functionText, 2));
            for entry = 1:numel(metadata.functionText)
                text = metadata.functionText{entry};
                assert(ischar(text) && (isrow(text) || isempty(text)), ...
                    'MdeMat:BasisResource:InvalidExpression', ...
                    'Resource expressions must be character vectors.');
                if isempty(text)
                    continue;
                end
                compiled = mdemat.common.approximation.CompiledFunction.compileExpression(text, nDims);
                [mode, derivative] = ind2sub(size(metadata.functionText), entry);
                if mode <= nCodims
                    handles{mode, derivative} = compiled;
                end
            end
        end
    end

    methods (Static, Access = private)
        function handle = compileExpression(text, nDims)
            % Only scalar literals, declared inputs, pi, arithmetic and the
            % reviewed elementary functions are allowed. No indexing, text,
            % assignment, nested handles, field access or general calls.
            signature = regexp(text, '^@\(([^)]*)\)(.+)$', 'tokens', 'once');
            assert(~isempty(signature), 'MdeMat:BasisResource:InvalidExpression', ...
                'Expected a self-contained anonymous numeric function.');
            inputs = strtrim(strsplit(signature{1}, ','));
            functions = {'sin', 'cos', 'sqrt'};
            validNames = cellfun(@(name) ~isempty(regexp(name, ...
                '^[A-Za-z]\w*$', 'once')), inputs);
            valid = numel(inputs) == nDims && all(validNames) && ...
                numel(unique(inputs)) == nDims && ...
                ~any(ismember(inputs, [functions, {'pi'}]));
            assert(valid, 'MdeMat:BasisResource:InvalidExpression', ...
                'Anonymous input names do not match the resource dimension.');
            body = regexprep(signature{2}, '\s+', '');
            pattern = ['(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?', ...
                '|[A-Za-z]\w*|\.\*|\./|\.\^|[+\-*/^()]'];
            tokens = regexp(body, pattern, 'match');
            valid = strcmp([tokens{:}], body);
            for index = 1:numel(tokens)
                token = tokens{index};
                if ~isletter(token(1))
                    continue;
                end
                isFunction = ismember(token, functions);
                isVariable = ismember(token, inputs) || strcmp(token, 'pi');
                nextIsParenthesis = index < numel(tokens) && strcmp(tokens{index+1}, '(');
                valid = valid && (isFunction || isVariable) && ...
                    (isFunction == nextIsParenthesis);
            end
            assert(valid, 'MdeMat:BasisResource:InvalidExpression', ...
                'Expression contains an undeclared variable or unsupported operation.');
            try
                handle = str2func(text);
            catch exception
                error('MdeMat:BasisResource:InvalidExpression', ...
                    'Malformed numeric expression: %s', exception.message);
            end
        end
    end

    methods (Access = protected)
        function Y = evalImpl(obj, X)
            % EVALIMPL Implementation of function evaluation.
            %
            %   Y = evalImpl(obj, X) evaluates the compiled function
            %   at the specified points using the loaded metadata handles.

            handles = obj.Metadata.handles;
            Y = zeros(obj.NCodims, size(X, 2));

            Z = mat2cell(X, ones(1, size(X, 1)), size(X, 2));

            for j = 1:size(handles, 1)
                f = handles{j, 1};
                if ~isempty(f), Y(j, :) = f(Z{:}); end
            end
        end

        function dY = diffImpl(obj, X, order)
            % DIFFIMPL Implementation of function derivative
            % evaluation.
            %
            %   dY = diffImpl(obj, X, order) computes derivatives
            %   of the compiled function at the specified points using
            %   precompiled derivative handles from the loaded metadata.

            handles = obj.Metadata.handles;

            maxOrder = size(handles, 2) - 1;

            m = sum(order);

            assert(all(m <= maxOrder, 'all'), 'mdemat:common:approximation:CompiledFunction:InvalidInput', ...
                'Derivative order %d exceeds maximum available order %d.', ...
                m, maxOrder);

            k = m + 1;
            dY = zeros(obj.NCodims, size(X, 2));
            Z = mat2cell(X, ones(1, size(X, 1)), size(X, 2));
            for j = 1:size(handles, 1)
                f = handles{j, k};
                if ~isempty(f), dY(j, :) = f(Z{:}); end
            end
        end
    end
end
