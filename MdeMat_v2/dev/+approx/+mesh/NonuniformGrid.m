classdef NonuniformGrid < approx.mesh.Grid
    % NONUNIFORMGRID Nonuniform multidimensional structured grid.
    %
    %   NonuniformGrid represents a structured grid with non-uniform
    %   spacing defined by coordinate vectors along each dimension. This
    %   class implements memory-efficient tensor-based storage while
    %   allowing arbitrary element sizes and distributions.
    %
    % Examples:
    %   % 1D nonuniform grid with refined regions
    %   x = [0, 0.1, 0.15, 0.2, 0.5, 0.8, 0.9, 1.0];
    %   grid1D = NonuniformGrid({x});
    %   
    %   % 2D nonuniform grid with different spacing in each direction
    %   x = 0:0.1:1;
    %   y = [0, 0.05, 0.2, 0.5, 1.0];
    %   grid2D = NonuniformGrid({x, y});
    %
    % Notes:
    %   The grid is defined by cell arrays of coordinate vectors, one per
    %   dimension. Element centroids and spacings are computed
    %   automatically from the coordinate vectors.
    %
    % See also:
    %   approx.mesh.Grid, approx.mesh.UniformGrid,
    %   approx.mesh.SeparableGraph

    methods
        function obj = NonuniformGrid(X)
            % NONUNIFORMGRID Constructor for NonuniformGrid.
            %
            %   obj = NonuniformGrid(X) creates a nonuniform grid from
            %   the specified coordinate vectors.
            %
            % Inputs:
            %   X - Cell array of coordinate vectors, one per dimension.
            %
            % Outputs:
            %   obj - Constructed NonuniformGrid object

            core.except.assert(iscell(X), 'InvalidInput', ...
                'Input X must be a cell array.');

            obj@approx.mesh.Grid(length(X));
            obj.initialize(X);
        end

        function Y = collocate(obj, X, I)
            % COLLOCATE Map reference points to physical coordinates.
            %
            %   Y = collocate(obj, X) maps reference points X from the
            %   reference element [-1,1]^d to physical grid coordinates.
            %
            %   Y = collocate(obj, X, I) maps points to specific elements
            %   identified by multi-indices I.
            %
            % Inputs:
            %   obj - The NonuniformGrid object
            %   X - Reference points (cell array or matrix)
            %   I - Element multi-indices (optional, matrix)
            %
            % Outputs:
            %   Y - Physical coordinates (cell array or matrix)

            if nargin < 3, I = []; end

            d = obj.nDims;
            if ~isempty(I)
                core.except.assert(size(I, 2) == d, ...
                    'InvalidInput', 'Indices dimension mismatch.');
            end

            if ismatrix(X) && ~iscell(X)
                X = num2cell(reshape(X, d, []), 2);
            end

            Y = cell(d, 1);
            for i = 1:d
                if isempty(I)
                    a = obj.centroids{i};
                    h = obj.spacings{i};
                    b = bsxfun(@times, h, X{i}(:));
                else
                    a = obj.centroids{i}(I(:, i));
                    h = obj.spacings{i};
                    b = bsxfun(@times, h(I(:, i)), X{i}(:));
                end
                Y{i} = reshape(bsxfun(@plus, a, b), 1, []);
            end

            if ~isempty(I), Y = cell2mat(Y); end
        end
    
        function newObj = refine(obj, k)
            % REFINE Refine the nonuniform grid.
            %
            %   newObj = refine(obj, k) refines the grid by subdividing each
            %   element according to the subdivision pattern specified in
            %   @a k.
            %
            % Inputs:
            %   obj - The NonuniformGrid object
            %   k - Number of subdivisions per element in each dimension. 
            %
            % Outputs:
            %   newObj - The refined NonuniformGrid object

            d = obj.nDims;
            x = cell(1, d);
            for i = 1:d
                y = obj.centroids{i};
                r = obj.resolution(i);
                s = obj.spacings{i};
                a = y - s/2;
                b = y + s/2;
                h = (b-a) ./ k{i};
                x{i} = zeros(1, sum(k{i})+1);
                x{i}(1:(k{i}(1)+1)) = a(1) : h(1): b(1);
                l = k{i}(1) + 1;
                for j = 2:r
                    x{i}(l+(1:k{i}(j))) = a(j) + h(j):h(j):b(j);
                    l = l + k{i}(j);
                end
            end
            newObj = approx.mesh.NonuniformGrid(x);
        end
    end

    methods (Access = protected)
        function obj = initialize(obj, X)
            % INITIALIZE Initialize the nonuniform grid from coordinates.
            %
            %   obj = initialize(obj, X) sets up the grid structure from
            %   coordinate vectors, computing element centroids, spacings,
            %   and other grid properties.
            %
            % Inputs:
            %   obj - The NonuniformGrid object
            %   X - Cell array of coordinate vectors
            %
            % Outputs:
            %   obj - The NonuniformGrid object
            
            d = obj.nDims;
            obj.resolution = zeros(1, d);
            obj.centroids = cell(1, d);
            obj.spacings = cell(1, d);

            for i = 1:d
                x = X{i};
                core.except.assert( ...
                    isvector(x) && ~isempty(x), ...
                    'InvalidInput', ...
                    'Each element of X must be a non-empty vector.');
                obj.resolution(i) = numel(x) - 1;
                obj.spacings{i} = diff(x);
                obj.centroids{i} = (x(1:end-1) + x(2:end)) / 2;
            end

            obj.bbox = cellfun(@(x) [min(x), max(x)], X, 'Un', 0);
            obj.bbox = horzcat(obj.bbox{:});
            obj.setElements();
            obj.setBoundary();
            obj.setIndexer();
        end

        function h = getMinSpacing(obj)
            % GETMINSPACING Get the minimum element spacing.
            %
            %   h = getMinSpacing(obj) returns the minimum spacing among
            %   all elements in the grid across all dimensions.
            %
            % Inputs:
            %   obj - The NonuniformGrid object
            %
            % Outputs:
            %   h - Minimum element spacing (positive scalar)

            h = min(cellfun(@min, obj.spacings));
        end

        function h = getAllElementMagnitudes(obj)
            % GETALLELEMENTMAGNITUDES Get element magnitudes.
            %
            %   h = getAllElementMagnitudes(obj) computes the magnitude
            %   (volume/area/length) of each element in the grid.
            %
            % Inputs:
            %   obj - The NonuniformGrid object
            %
            % Outputs:
            %   h - Vector of element magnitudes

            h = obj.spacings;
            switch obj.indexer.style
                case 'C'
                    h = core.linalg.kronecker(h{:});
                case 'F'
                    h = core.linalg.kronecker(h{end:-1:1});
            end
            h = h(:);
        end

        function G = getSeparableGraph(obj)
            % GETSEPARABLEGRAPH Get the underlying separable graph.
            %
            %   G = getSeparableGraph(obj) creates a separable graph
            %   representation of the grid connectivity structure.
            %
            % Inputs:
            %   obj - The NonuniformGrid object
            %
            % Outputs:
            %   G - SeparableGraph representing grid connectivity

            d = obj.nDims;
            h = obj.spacings;
            c = obj.centroids;
            V = arrayfun(@(i) [c{i} - h{i} / 2, c{i}(end) + h{i}(end) / 2].', 1:d, 'Un', 0);
            E = arrayfun(@(n) [(1:n).', (2:(n + 1)).'], obj.resolution, 'Un', 0);
            G = approx.mesh.SeparableGraph(V, E);
        end
    end
end