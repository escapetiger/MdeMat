classdef LxMultiIndexer < core.linalg.CachedMultiIndexer
    % LXMULTIINDEXER Generate L-infinity constrained multi-indices.
    %
    %   LxMultiIndexer generates one-based multi-indices such that maximum
    %   value in each index is within specified bounds. Useful for tensor
    %   approximation with L-infinity norm constraints.
    %
    % See Also:
    %   core.linalg.CachedMultiIndexer, core.linalg.L1MultiIndexer

    methods
        function obj = setCache(obj, nDims, maxIdx, options)
            % SETCACHE Generate and cache L-infinity constrained
            % multi-indices.
            %
            %   obj = setCache(obj, nDims, maxIdx) generates all
            %   multi-indices with @a nDims dimensions and max value <= @a
            %   maxIdx.
            %
            %   obj = setCache(obj, nDims, maxIdx, minIdx=minIdx) generates
            %   indices with values between @a minIdx and @a maxIdx.

            arguments
                obj core.linalg.LxMultiIndexer
                nDims{mustBeInteger, mustBePositive}
                maxIdx(1, :) {mustBeInteger, mustBePositive}
                options.minIdx(1, :) {mustBeInteger, mustBePositive} = []
            end

            nd = nDims;
            b = maxIdx;

            if isempty(options.minIdx)
                a = 1;
            else
                a = options.minIdx;
            end

            if isscalar(a), a = repmat(a, 1, nd); end
            if isscalar(b), b = repmat(b, 1, nd); end

            core.except.assert(all(b >= a), 'InvalidInput', ...
                'maxIdx must be greater than or equal to minIdx.');

            if nd == 1
                obj.Cache = (a(1):b(1)).';
                return;
            end

            N = arrayfun(@(i) a(i):b(i), 1:nd, 'Un', 0);
            if strcmp(obj.Style, 'F')
                [N{:}] = ndgrid(N{:});
            else
                [N{end:-1:1}] = ndgrid(N{end:-1:1});
            end

            N = cellfun(@(x) x(:), N, 'Un', 0);
            obj.Cache = cat(2, N{:});
        end
    end
end