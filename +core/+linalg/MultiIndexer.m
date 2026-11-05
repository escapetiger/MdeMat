classdef MultiIndexer < core.linalg.Indexer
    % MULTIINDEXER Indexers handling multi-dimensional indexing.
    %
    %   MultiIndexer generates multi-indices for tensor operations and
    %   provides methods to convert between multi-indices and linear
    %   indices.
    %
    % Examples:
    %   % Create 3D multi-indexer and generate indices
    %   indexer = core.linalg.MultiIndexer([2,3,4]);
    %   M = indexer.generate();
    %   L = indexer.multiToLinear(M);
    %
    % See Also:
    %   core.linalg.Indexer, core.linalg.CachedMultiIndexer

    properties
        shape % Tensor dimensions as row vector of positive integers
    end

    properties (Dependent)
        nDims % Number of dimensions, derived from shape
    end

    properties (Access = protected)
        strides % Tensor strides for linear index computation
    end

    methods
        function obj = MultiIndexer(shape, style)
            % MULTIINDEXER Constructor for multi-dimensional indexer.
            %
            %   obj = MultiIndexer() creates empty MultiIndexer.
            %
            %   obj = MultiIndexer(shape) creates with specified @a shape.
            %
            %   obj = MultiIndexer(shape, style) creates with @a shape and
            %   @a style.
            %
            % Inputs:
            %   shape - Tensor shape (optional)
            %   style - Storage ordering style (optional, default: 'F')
            %
            % Outputs:
            %   obj - Constructed MultiIndexer object

            if nargin < 2, style = 'F'; end
            if nargin < 1, shape = []; end

            obj@core.linalg.Indexer(style);

            if ~isempty(shape), obj.setShape(shape); end
        end

        function n = get.nDims(obj)
            % GET.NDIMS Get number of dimensions based on shape.

            n = length(obj.shape);
        end

        function obj = setShape(obj, shape)
            % SETSHAPE Set tensor shape and stride for indexer.
            %
            %   obj = setShape(obj, shape) validates and sets @a shape
            %   property, then updates strides based on storage style.
            %
            % Inputs:
            %   obj - The MultiIndexer object
            %   shape - Tensor shape
            %
            % Outputs:
            %   obj - The MultiIndexer object

            core.except.assert(all(shape > 0), ...
                'InvalidShape', ...
                'Shape must be a vector of positive integers.');

            obj.shape = shape(:).';
            obj.setStrides();
        end

        function M = generate(obj)
            % GENERATE Return all one-based multi-indices.
            %
            %   M = generate(obj) returns matrix where each row is
            %   multi-index satisfying 1 <= M(i,j) <= obj.shape(j).
            %
            % Inputs:
            %   obj - The MultiIndexer object
            %
            % Outputs:
            %   M - Matrix where each row represents a multi-index

            core.except.assert(~isempty(obj.shape), ...
                'MissingShape', ...
                'Shape must be set before using this method.');

            if obj.nDims == 1
                M = (1:obj.shape).';
                return;
            end

            N = arrayfun(@(x) 1:x, obj.shape, 'Un', 0);
            if strcmp(obj.style, 'F')
                [N{:}] = ndgrid(N{:});
            else
                [N{end:-1:1}] = ndgrid(N{end:-1:1});
            end
            N = cellfun(@(x) x(:), N, 'Un', 0);
            M = cat(2, N{:});
        end

        function L = multiToLinear(obj, M, bc)
            % MULTITOLINEAR Convert multi-indices to linear indices.
            %
            %   L = multiToLinear(obj, M) converts multi-indices @a M to
            %   linear indices based on tensor dimensions.
            %
            %   L = multiToLinear(obj, M, bc) uses boundary condition @a
            %   bc.
            %
            % Inputs:
            %   obj - The MultiIndexer object
            %   M - Multi-indices (matrix or cell array)
            %   bc - Boundary condition (optional, default: 1)
            %        * 0: periodic boundary (wrap-around)
            %        * 1: strict boundary (default, out-of-range set to 0)
            %
            % Outputs:
            %   L - Linear indices (vector)

            if nargin < 3, bc = 1; end

            core.except.assert(~isempty(obj.shape), ...
                'MissingShape', ...
                'Shape must be set before using this method.');

            if iscell(M)
                M = obj.factorToMulti(M);
            else
                core.except.assert( ...
                    ismatrix(M) && size(M, 2) == obj.nDims, ...
                    'InvalidInput', ...
                    'Multi-indices must has %d columns.', obj.nDims);
            end

            if bc == 0
                M = mod(M - 1, obj.shape) + 1;
            end

            L = (M - 1) * obj.strides(:) + 1;

            if bc == 1
                L(any(M < 1, 2) | any(M > obj.shape, 2)) = 0;
            end
        end

        function M = linearToMulti(obj, L)
            % LINEARTOMULTI Convert linear indices to multi-indices.
            %
            %   M = linearToMulti(obj, L) converts linear indices @a L to
            %   multi-indices for tensor with specified shape.
            %
            % Inputs:
            %   obj - The MultiIndexer object
            %   L - Linear indices (vector)
            %
            % Outputs:
            %   M - Multi-indices (matrix)

            core.except.assert(~isempty(obj.shape), ...
                'MissingShape', ...
                'Shape must be set before using this method.');

            upper = prod(obj.shape);
            core.except.assert(all(L <= upper) && all(L >= 1), ...
                'IndexOutOfBounds', ...
                'Linear index out of range (min: 1, max: %d).', upper);

            M = mod(floor((L(:) - 1)./obj.strides), obj.shape) + 1;
        end

        function M = factorToMulti(obj, F)
            % FACTORTOMULTI Convert factor indices to multi-indices.
            %
            %   M = factorToMulti(obj, F) converts factor indices @a F to
            %   multi-indices for tensor with specified shape.
            %
            % Inputs:
            %   obj - The MultiIndexer object
            %   F - Cell array where each entry is factor indices
            %
            % Outputs:
            %   M - Multi-indices (matrix)

            core.except.assert(length(F) == obj.nDims, ...
                'DimensionMismatch', ...
                'Factor number must be (%d)', obj.nDims);

            if strcmp(obj.style, 'F')
                [F{:}] = ndgrid(F{:});
            else
                [F{end:-1:1}] = ndgrid(F{end:-1:1});
            end
            M = cell2mat(cellfun(@(x) x(:), F, 'Un', 0));
        end
    end

    methods (Access = protected)
        function obj = setStrides(obj)
            % SETSTRIDES Set tensor strides for indexer.

            if isempty(obj.shape)
                obj.strides = [];
                return;
            end

            obj.strides = ones(1, obj.nDims);

            if strcmp(obj.style, 'F')
                obj.strides(2:end) = cumprod(obj.shape(1:end-1), 2);
            else
                obj.strides(1:end-1) = flip(cumprod(flip(obj.shape(2:end), 2), 2), 2);
            end
        end
    end
end
