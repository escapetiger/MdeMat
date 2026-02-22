classdef CompiledFunction < core.function.Function
    % COMPILEDFUNCTION Base class for all compiled functions.
    %
    %   CompiledFunction provides a bridge between symbolic function
    %   definitions and efficient numerical evaluation by automatically
    %   loading precompiled function handles and metadata from the
    %   configuration directory.
    %
    %   The class supports automatic metadata discovery, caching, and
    %   efficient function evaluation using precompiled function handles
    %   stored in MAT files.
    %
    % See Also:
    %   core.function.Function

    properties (Constant, Access = protected)
        MetaPath = '+core/+symbolic/metadata' % Directory containing metadata files
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

            obj@core.function.Function(nDims=nDims, nCodims=nCodims);
            obj.load(metaName);
        end

        function obj = load(obj, fileName)
            % LOAD Load compiled metadata from specified file.
            %
            %   obj = load(obj, fileName) loads precompiled function
            %   metadata from the specified file and updates the object's
            %   function dimensions based on the loaded metadata.

            arguments
                obj core.function.CompiledFunction
                fileName {mustBeTextScalar}
            end

            fileName = char(fileName);
            if ~endsWith(fileName, obj.MetaExt)
                fileName = [fileName, obj.MetaExt];
            end
            fullPath = fullfile(obj.MetaPath, fileName);
            
            core.except.assert(exist(fullPath, 'file') == 2, ...
                'FileNotFound', 'Metadata file not found: %s', fullPath);

            data = load(fullPath, 'metadata');

            core.except.assert(isfield(data, 'metadata'), ...
                'InvalidFile', 'Metadata file must contain a metadata field.');

            nd = obj.NDims;
            nc = obj.NCodims;
            metadata = data.metadata;
            metadata.nDims = nd;
            metadata.NCodims = nc;
            metadata.variables = metadata.variables(1:nd);
            metadata.expressions = metadata.expressions(1:nc);
            metadata.handles = metadata.handles(1:nc, 1:max(2, nc));
            obj.Metadata = metadata;
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

            core.except.assert(m <= maxOrder, 'InvalidInput', ...
                'Derivative order %d exceeds maximum available order %d.', ...
                m, maxOrder);

            k = m + 1;
            dY = zeros(obj.NCodims, size(X, 2));
            for j = 1:size(handles, 1)
                f = handles{j, k};
                if ~isempty(f), dY(j, :) = f(X); end
            end
        end
    end
end