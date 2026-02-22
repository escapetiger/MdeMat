classdef L1MultiIndexer < core.linalg.CachedMultiIndexer
    % L1MULTIINDEXER Generate L1-constrained multi-indices.
    %
    %   L1MultiIndexer generates one-based multi-indices such that sum of
    %   indices is within specified bounds. Useful for sparse tensor
    %   approximation with L1-norm constraints.
    %
    % Examples:
    %   % Create L1-constrained indexer
    %   indexer = core.linalg.L1MultiIndexer();
    %   
    %   % Generate 2D indices with sum <= 4
    %   indexer.setCache(2, 4);
    %   
    %   % Generate 3D indices with sum between 3 and 5
    %   indexer.setCache(3, 3, 5);
    %
    % See Also:
    %   core.linalg.CachedMultiIndexer, core.linalg.LxMultiIndexer
    
    methods
        function obj = L1MultiIndexer(style)
            % L1MULTIINDEXER Constructor for L1MultiIndexer.
            %
            %   obj = L1MultiIndexer() creates L1MultiIndexer with default
            %   storage ordering style.
            %
            %   obj = L1MultiIndexer(style) creates with specified @a style.
            %
            % Inputs:
            %   style - Storage ordering style (optional, default: 'F')
            %
            % Outputs:
            %   obj - Constructed L1MultiIndexer object

            if nargin < 1, style = 'F'; end
            obj@core.linalg.CachedMultiIndexer(style);
        end

        function obj = setCache(obj, varargin)
            % SETCACHE Generate and cache L1-constrained multi-indices.
            %
            %   obj = setCache(obj, nDims, sMax) generates all
            %   multi-indices with @a nDims dimensions and sum <= @a sMax.
            %
            %   obj = setCache(obj, nDims, sMin, sMax) generates indices
            %   with sum between @a sMin and @a sMax.
            %
            % Inputs:
            %   obj - The L1MultiIndexer object
            %   varargin - Input arguments
            %<   nDims - Number of dimensions for indices
            %<   sMax - Maximum sum of indices (must be >= nDims)
            %<   sMin - Minimum sum of indices (optional, default: nDims)
            %
            % Outputs:
            %   obj - The L1MultiIndexer object
            
            core.except.assert(length(varargin) >= 2, ...
                'InvalidInput', ...
                'Expected at least nDims and sMax arguments.');
            
            d = varargin{1};
            
            core.except.assert(d >= 1, 'InvalidInput', ...
                'Number of dimensions must be a positive integer.');
            
            if length(varargin) == 2
                b = varargin{2};
                a = d;
            else
                a = varargin{2};
                b = varargin{3};
            end
            
            core.except.assert(all(b >= 1), ...
                'InvalidInput', 'Upper bound must be a positive integer.');
                
            core.except.assert(all(a >= 1), ...
                'InvalidInput', 'Lower bound must be a positive integer.');
            
            core.except.assert(a >= d, 'InvalidInput', ...
                'sMin must be greater than or equal to nDims.');
            
            core.except.assert(b >= a, 'InvalidInput', ...
                'sMax must be greater than or equal to sMin.');
            
            if d == 1
                obj.cache = (1:b).';
                return;
            end
            
            m = 0;
            for n = a:b
                m = m + nchoosek(n - 1, d - 1);
            end
            
            M = zeros(m, d);
            i = 0;
            for s = a:b
                D = nchoosek(1:(s-1), d-1);
                G = [D(:, 1), diff(D, 1, 2), s - D(:, end)];
                if strcmp(obj.style, 'C')
                    [G, ~] = sortrows(G, 1:d);
                else
                    [G, ~] = sortrows(G, d:-1:1);
                end
                k = size(G, 1);
                M(i + (1:k), :) = G;
                i = i + k;
            end
            obj.cache = M;
        end
    end
end