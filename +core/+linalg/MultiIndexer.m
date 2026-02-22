classdef MultiIndexer < core.linalg.Indexer
    % MULTIINDEXER Indexers handling multi-dimensional indexing.
    %
    %   MultiIndexer generates multi-indices for tensor operations and
    %   provides methods to convert between multi-indices and linear
    %   indices. Supports both column-major ('F') and row-major ('C')
    %   storage orders with various boundary conditions for index
    %   conversion.
    %
    % See also:
    %   core.linalg.Indexer, core.linalg.CachedMultiIndexer

    properties
        ValidPads = {'empty', 'wrap', 'edge', 'reflect', 'symmetric'}
    end

    properties
        Shape % Tensor dimensions as row vector of positive integers
    end

    properties (Dependent)
        NDims % Number of dimensions, derived from Shape
    end

    properties (Access = protected)
        Strides % Tensor strides for linear index computation
    end

    methods
        function obj = MultiIndexer(options)
            % MULTIINDEXER Construct an instance of MultiIndexer.
            %
            %   obj = MultiIndexer() creates empty MultiIndexer.
            %
            %   obj = MultiIndexer(Name=Value) creates with specified
            %   options.

            arguments
                options.shape(1, :) {mustBeInteger, mustBePositive} = []
                options.style{mustBeTextScalar} = 'F'
            end

            obj@core.linalg.Indexer(style = options.style);
            obj.setShape(options.shape);
        end

        function n = get.NDims(obj)
            % GET.NDIMS Returns the value of the dependent property
            % 'NDims'.

            n = length(obj.Shape);
        end

        function obj = setShape(obj, shape)
            % SETSHAPE Set tensor shape and stride for indexer.

            arguments
                obj core.linalg.MultiIndexer
                shape(1, :) {mustBeInteger, mustBePositive}
            end

            obj.Shape = shape;
            obj.updateStrides();
        end

        function M = generate(obj)
            % GENERATE Return all one-based multi-indices.
            %
            %   M = generate(obj) returns matrix where each row is
            %   multi-index satisfying 1 <= M(i,j) <= obj.Shape(j).

            arguments
                obj core.linalg.MultiIndexer
            end

            core.except.assert(~isempty(obj.Shape), 'MissingShape', ...
                'Shape must be set before using this method.');

            if obj.NDims == 1
                M = (1:obj.Shape).';
                return;
            end

            N = arrayfun(@(x) 1:x, obj.Shape, 'Un', 0);
            if strcmp(obj.Style, 'F')
                [N{:}] = ndgrid(N{:});
            else
                [N{end:-1:1}] = ndgrid(N{end:-1:1});
            end
            N = cellfun(@(x) x(:), N, 'Un', 0);
            M = cat(2, N{:});
        end

        function L = multiToLinear(obj, M, options)
            % MULTITOLINEAR Convert multi-indices to linear indices.
            %
            %   L = multiToLinear(obj, M) converts multi-indices @a M to
            %   linear indices based on tensor dimensions.
            %
            %   L = multiToLinear(obj, M, Name=Value) uses boundary
            %   condition specified by options.
            %
            % Notes:
            %   Padding modes include 'empty' (out-of-range set to 0),
            %   'wrap' (wraps around), 'edge' (clamps to edges),
            %   'reflect' (reflects excluding edges), and 'symmetric'
            %   (reflects including edges).

            arguments
                obj core.linalg.MultiIndexer
                M {mustBeInteger}
                options.pad{mustBeTextScalar} = 'empty'
            end

            core.except.assert(ismember(options.pad, obj.ValidPads), ...
                'InvalidInput', 'Padding mode is invalid.');

            core.except.assert(~isempty(obj.Shape), 'MissingShape', ...
                'Shape must be set before using this method.');

            nd = obj.NDims;
            core.except.assert(ismatrix(M) && size(M, 2) == nd, ...
                'InvalidInput', 'Multi-indices must have %d columns.', nd);

            ns = obj.Shape;
            switch lower(options.pad)
                case 'wrap'
                    M = mod(M - 1, ns) + 1;
                case 'edge'
                    M = min(M, ns);
                    M = max(M, ones(1, nd));
                case 'reflect'
                    period = 2 * ns - 2;
                    P = mod(M - 1, period);
                    mask = P >= ns;
                    Q = repmat(period, size(P, 1), 1);
                    P(mask) = Q(mask) - P(mask);
                    M = P + 1;
                case 'symmetric'
                    period = 2 * ns;
                    P = mod(M - 1, period);
                    mask = P >= ns;
                    Q = repmat(period-1, size(P, 1), 1);
                    P(mask) = Q(mask) - P(mask);
                    M = P + 1;
                otherwise
            end
            L = (M - 1) * obj.Strides(:) + 1;
            L(any(M < 1, 2) | any(M > obj.Shape, 2)) = 0;
        end

        function M = linearToMulti(obj, L)
            % LINEARTOMULTI Convert linear indices to multi-indices.
            %
            %   M = linearToMulti(obj, L) converts linear indices @a L to
            %   multi-indices for tensor with specified shape.
            %
            % Notes:
            %   Linear indices must be in range [1, prod(Shape)]. Returns
            %   matrix where each row represents a multi-index.

            arguments
                obj core.linalg.MultiIndexer
                L {mustBePositive, mustBeInteger}
            end

            core.except.assert(~isempty(obj.Shape), 'MissingShape', ...
                'Shape must be set before using this method.');

            upper = prod(obj.Shape);
            core.except.assert(all(L <= upper), 'IndexOutOfBounds', ...
                'Linear index out of range (min: 1, max: %d).', upper);

            M = mod(floor((L(:) - 1)./obj.Strides), obj.Shape) + 1;
        end

        function M = factorToMulti(obj, F)
            % FACTORTOMULTI Convert factor indices to multi-indices.
            %
            %   M = factorToMulti(obj, F) converts factor indices @a F to
            %   multi-indices for tensor with specified shape.
            %
            % Notes:
            %   Factor indices @a F must be a cell array with length equal
            %   to NDims. Each cell contains indices for the corresponding
            %   dimension.

            arguments
                obj core.linalg.MultiIndexer
                F {mustBeA(F, 'cell')}
            end

            nd = obj.NDims;
            core.except.assert(length(F) == nd, 'DimensionMismatch', ...
                'Factor number must be (%d)', nd);

            if strcmp(obj.Style, 'F')
                [F{:}] = ndgrid(F{:});
            else
                [F{end:-1:1}] = ndgrid(F{end:-1:1});
            end
            M = cell2mat(cellfun(@(x) x(:), F, 'Un', 0));
        end
    end

    methods (Access = protected)
        function obj = updateStrides(obj)
            % UPDATESTRIDES Update tensor strides for indexer.

            if isempty(obj.Shape)
                obj.Strides = [];
                return;
            end

            obj.Strides = ones(1, obj.NDims);

            if strcmp(obj.Style, 'F')
                obj.Strides(2:end) = cumprod(obj.Shape(1:end-1), 2);
            else
                obj.Strides(1:end-1) = flip(cumprod(flip(obj.Shape(2:end), ...
                    2), 2), 2);
            end
        end
    end
end
