classdef L1MultiIndexer < core.linalg.CachedMultiIndexer
    % L1MULTIINDEXER Generate L1-constrained multi-indices.
    %
    %   L1MultiIndexer generates one-based multi-indices such that sum of
    %   indices is within specified bounds. Useful for sparse tensor
    %   approximation with L1-norm constraints.
    %
    % See Also:
    %   core.linalg.CachedMultiIndexer, core.linalg.LxMultiIndexer

    methods
        function obj = setCache(obj, nDims, maxSum, options)
            % SETCACHE Generate and cache L1-constrained multi-indices.
            %
            %   obj = setCache(obj, nDims=nDims, maxSum=maxSum)
            %   generates all multi-indices with @a nDims dimensions and
            %   sum <= @a sMax.
            %
            %   obj = setCache(obj, nDims=nDims, maxSum=maxSum,
            %   minSum=minSum) generates  multi-indices with @a nDims
            %   dimensions and sum between @a minSum and @a maxSum.

            arguments
                obj core.linalg.L1MultiIndexer
                nDims{mustBeInteger, mustBePositive}
                maxSum{mustBeInteger, mustBePositive}
                options.minSum{mustBeScalarOrEmpty, mustBeInteger} = []
            end

            nd = nDims;
            b = maxSum;

            if isempty(options.minSum)
                a = nd;
            else
                a = options.minSum;
                core.except.assert(a >= nd, 'InvalidInput', ...
                    'minSum must be greater than or equal to nDims.');
            end

            core.except.assert(b >= a, 'InvalidInput', ...
                'maxSum must be greater than or equal to minSum.');

            if nd == 1
                obj.Cache = (a:b).';
                return;
            end

            nm = 0;
            for i = a:b
                nm = nm + nchoosek(i-1, nd-1);
            end

            M = zeros(nm, nd);
            i = 0;
            for s = a:b
                D = nchoosek(1:(s - 1), nd-1);
                G = [D(:, 1), diff(D, 1, 2), s - D(:, end)];
                if strcmp(obj.Style, 'C')
                    [G, ~] = sortrows(G, 1:nd);
                else
                    [G, ~] = sortrows(G, nd:-1:1);
                end
                k = size(G, 1);
                M(i + (1:k), :) = G;
                i = i + k;
            end
            obj.Cache = M;
        end
    end
end