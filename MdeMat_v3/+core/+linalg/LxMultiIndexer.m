classdef LxMultiIndexer < core.linalg.CachedMultiIndexer
    % LXMULTIINDEXER Generate L-infinity constrained multi-indices.
    %
    %   LxMultiIndexer generates one-based multi-indices such that maximum
    %   value in each index is within specified bounds. Useful for tensor
    %   approximation with L-infinity norm constraints.
    %
    % Examples:
    %   % Create L-infinity constrained indexer
    %   indexer = core.linalg.LxMultiIndexer();
    %   
    %   % Generate 2D indices with max <= 2
    %   indexer.setCache(2, 2);
    %   
    %   % Generate 3D indices with each component between 2 and 3
    %   indexer.setCache(3, 2, 3);
    %
    % See Also:
    %   core.linalg.CachedMultiIndexer, core.linalg.L1MultiIndexer
    
    methods
        function obj = LxMultiIndexer(style)
            % LXMULTIINDEXER Constructor for LxMultiIndexer.
            %
            %   obj = LxMultiIndexer() creates LxMultiIndexer with default
            %   storage ordering style.
            %
            %   obj = LxMultiIndexer(style) creates with specified @a
            %   style.
            %
            % Inputs:
            %   style - Storage ordering style (optional, default: 'F')
            %
            % Outputs:
            %   obj - Constructed LxMultiIndexer object

            if nargin < 1, style = 'F'; end
            obj@core.linalg.CachedMultiIndexer(style);
        end
        
        function obj = setCache(obj, varargin)
            % SETCACHE Generate and cache L-infinity constrained
            % multi-indices.
            %
            %   obj = setCache(obj, nDims, iMax) generates all multi-indices
            %   with @a nDims dimensions and max value <= @a iMax.
            % 
            %   obj = setCache(obj, nDims, iMin, iMax) generates indices with
            %   values between @a iMin and @a iMax.
            %
            % Inputs:
            %   obj - The LxMultiIndexer object
            %   varargin - Input arguments
            %<   nDims - Number of dimensions for indices
            %<   iMax - Maximum index value (positive integer)
            %<   iMin - Minimum index value (optional, default: 1)
            %
            % Outputs:
            %   obj - The LxMultiIndexer object
            
            core.except.assert(length(varargin) >= 2, ...
                'InvalidInput', ...
                'Expected at least nDims and iMax arguments.');
            
            d = varargin{1};
            
            core.except.assert(d >= 1, ...
                'InvalidInput', ...
                'Number of dimensions must be a positive integer.');
            
            if length(varargin) == 2
                b = varargin{2};
                a = ones(1, d);
            else
                a = varargin{2};
                b = varargin{3};
            end

            if isscalar(a), a = repmat(a, 1, d); end
            if isscalar(b), b = repmat(b, 1, d); end

            core.except.assert(all(a > 0), ...
                'InvalidInput', ...
                'Lower bounds must be positive integers.');
                
            core.except.assert(all(b > 0), ...
                'InvalidInput', ...
                'Upper bounds must be positive integers.');

            if d == 1
                obj.cache = (a(1):b(1)).';
                return;
            end
            
            N = arrayfun(@(i) a(i):b(i), 1:d, 'Un', 0);
            if strcmp(obj.style, 'F')
                [N{:}] = ndgrid(N{:});
            else
                [N{end:-1:1}] = ndgrid(N{end:-1:1});
            end

            N = cellfun(@(x) x(:), N, 'Un', 0);
            obj.cache = cat(2, N{:});
        end
    end
end