classdef CachedMultiIndexer < core.linalg.Indexer
    % CACHEDMULTIINDEXER Base class for all indexers with cached indices.
    %
    %   CachedMultiIndexer provides foundation for index generation
    %   strategies that benefit from caching generated indices. Subclasses
    %   implement specific index generation while inheriting caching
    %   functionality.
    %
    % See Also:
    %   core.linalg.Indexer, core.linalg.MultiIndexer

    properties (Dependent)
        nDims % Number of dimensions, derived from cache
        nIndices % Number of cached indices, derived from cache
    end

    properties
        cache = [] % Cache matrix where each row represents a multi-index
    end

    methods
        function obj = CachedMultiIndexer(style)
            % CACHEDMULTIINDEXER Constructor for cached multi-indexer.
            %
            %   obj = CachedMultiIndexer() creates CachedMultiIndexer with
            %   default storage ordering style. 
            % 
            %   obj = CachedMultiIndexer(style) creates CachedMultiIndexer
            %   with specified @a style.
            %
            % Inputs:
            %   style - Storage ordering style (optional, default: 'F')
            %
            % Outputs:
            %   obj - Constructed CachedMultiIndexer object

            if nargin < 1, style = 'F'; end
            obj@core.linalg.Indexer(style);
        end

        function n = get.nDims(obj)
            % GET.NDIMS Get number of dimensions based on cache.

            n = size(obj.cache, 2);
        end

        function n = get.nIndices(obj)
            % GET.NINDICES Get number of indices based on cache.

            n = size(obj.cache, 1);
        end

        function L = multiToLinear(obj, M)
            % MULTITOLINEAR Convert multi-indices to linear indices.
            %
            %   L = multiToLinear(obj, M) converts multi-indices @a M to
            %   linear indices within cached set of indices.
            %
            % Inputs:
            %   obj - The CachedMultiIndexer object
            %   M - Matrix where each row is a multi-index
            %
            % Outputs:
            %   L - Column vector of linear indices

            core.except.assert( ...
                ismatrix(M) && size(M, 2) == obj.nDims, ...
                'InvalidInput', ...
                'Multi-indices must has %d columns.', obj.nDims);

            C = obj.cache;
            core.except.assert(~isempty(C), ...
                'MissingCache', 'Cache is empty.');

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
            %
            % Inputs:
            %   obj - The CachedMultiIndexer object
            %   L - Vector of linear indices
            %
            % Outputs:
            %   M - Matrix where each row is the corresponding multi-index

            C = obj.cache;
            core.except.assert(~isempty(C), 'MissingCache', ...
                'Cache is empty.');

            I = (L > 0) & (L <= size(C, 1));

            core.except.assert(all(I), 'InvalidInput', ...
                'Linear index %d is out of range (valid range: 1-%d).', ...
                L(find(~I, 1)), size(C, 1));

            M = C(L, :);
        end
    end

    methods (Abstract)
        % SETCACHE Generate and cache indices.
        obj = setCache(obj, varargin)
    end
end