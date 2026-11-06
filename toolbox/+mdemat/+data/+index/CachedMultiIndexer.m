classdef CachedMultiIndexer < mdemat.data.index.Indexer
    % CACHEDMULTIINDEXER Base class for all indexers with cached indices.
    %
    %   CachedMultiIndexer provides foundation for index generation
    %   strategies that benefit from caching generated indices. Subclasses
    %   implement specific index generation while inheriting caching
    %   functionality.
    %
    % See Also:
    %   mdemat.data.index.Indexer, mdemat.data.index.MultiIndexer

    properties (Dependent)
        NDims % Number of dimensions, derived from cache
        NIndices % Number of cached indices, derived from cache
    end

    properties
        Cache = [] % Cache matrix where each row represents a multi-index
    end

    methods
        function obj = CachedMultiIndexer(options)
            % CACHEDMULTIINDEXER Construct an instance of
            % CachedMultiIndexer.
            %
            %   obj = CachedMultiIndexer() creates CachedMultiIndexer with
            %   default storage ordering style.
            %
            %   obj = CachedMultiIndexer(Name=Value) creates
            %   CachedMultiIndexer with specified options.

            arguments
                options.style {mustBeTextScalar} = 'F'
            end

            obj@mdemat.data.index.Indexer(Style = options.style);
        end

        function n = get.NDims(obj)
            % GET.NDIMS Returns the value of the dependent property
            % 'NDims'.

            n = size(obj.Cache, 2);
        end

        function n = get.NIndices(obj)
            % GET.NINDICES Returns the value of the dependent property
            % 'NIndices'.

            n = size(obj.Cache, 1);
        end

        function L = multiToLinear(obj, M)
            % MULTITOLINEAR Convert multi-indices to linear indices.
            %
            %   L = multiToLinear(obj, M) converts multi-indices @a M to
            %   linear indices within cached set of indices.

            arguments
                obj mdemat.data.index.CachedMultiIndexer
                M {mustBeInteger}
            end

            nd = obj.NDims;
            assert(all(size(M, 2) == nd, 'all'), 'mdemat:data:index:CachedMultiIndexer:InvalidInput', ...
                'Multi-indices must has %d columns.', nd);

            C = obj.Cache;
            assert(all(~isempty(C), 'all'), ...
                'mdemat:data:index:CachedMultiIndexer:MissingCache', '%s', 'Cache is empty.');

            L = zeros(size(M, 1), 1);
            for i = 1:size(M, 1)
                I = M(i, :);
                J = find(all(C == repmat(I, size(C, 1), 1), 2));
                if ~isempty(J), L(i) = J(1); end
            end
        end

        function M = linearToMulti(obj, L)
            % LINEARTOMULTI Convert linear indices to cached multi-indices.
            %
            %   M = linearToMulti(obj, L) converts linear indices @a L to
            %   multi-indices from cached set.

            arguments
                obj mdemat.data.index.CachedMultiIndexer
                L {mustBePositive, mustBeInteger}
            end

            C = obj.Cache;
            assert(all(~isempty(C), 'all'), 'mdemat:data:index:CachedMultiIndexer:MissingCache', ...
                '%s', 'Cache is empty.');

            I = (L > 0) & (L <= size(C, 1));

            assert(all(all(I), 'all'), 'mdemat:data:index:CachedMultiIndexer:InvalidInput', ...
                'Linear index %d is out of range (valid range: 1-%d).', ...
                L(find(~I, 1)), size(C, 1));

            M = C(L, :);
        end
    end
end