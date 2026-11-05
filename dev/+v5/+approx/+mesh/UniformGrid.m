classdef UniformGrid < approx.mesh.Grid
    % UNIFORMGRID Uniform multidimensional structured grid.
    %
    %   UniformGrid represents a structured grid with uniform spacing
    %   defined by the number of elements per dimension and domain bounds.
    %   This class implements memory-efficient tensor-based storage with
    %   constant element sizes throughout the domain.
    %
    % Examples:
    %   % 2D uniform grid: 10×10 elements on unit square
    %   grid2D = UniformGrid([10, 10], [0, 0, 1, 1]);
    %   
    %   % 3D uniform grid: 5×8×12 elements on custom domain
    %   grid3D = UniformGrid([5, 8, 12], [-1, 0, -2, 1, 2, 3]);
    %   
    %   % Query grid properties
    %   h = grid2D.measure;  % Minimum element size
    %   nElems = grid2D.nTotalElements;
    %
    % Notes:
    %   Uniform grids have constant element spacing within each dimension.
    %   The bounding box format is [a1, b1, a2, b2, ..., ad, bd] where
    %   [ai, bi] defines the interval for dimension i.
    %
    % See also:
    %   approx.mesh.Grid, approx.mesh.NonuniformGrid,
    %   approx.mesh.SeparableGraph

    methods
        function obj = UniformGrid(n, bbox)
            % UNIFORMGRID Constructor for UniformGrid.
            %
            %   obj = UniformGrid(n, bbox) creates a uniform grid with
            %   the specified number of elements and domain bounds.
            %
            % Inputs:
            %   n - Vector specifying number of elements per dimension
            %   bbox - Bounding box coordinates [a1, b1, a2, b2, ...]
            %
            % Outputs:
            %   obj - Constructed UniformGrid obj.

            obj@approx.mesh.Grid(length(n));
            obj.bbox = bbox;

            core.except.assert(numel(bbox) == 2*obj.nDims, ...
                'DimensionMismatch', ...
                'bbox must have length 2*length(n).');

            core.except.assert(all(n > 0) && all(mod(n, 1) == 0), ...
                'InvalidInput', ...
                'n must contain positive integers.');

            n = n(:)';
            bbox = bbox(:)';

            a = bbox(1:2:end);
            b = bbox(2:2:end);
            a = a(:).';
            b = b(:).';

            core.except.assert(all(b > a), ...
                'InvalidInput', ...
                'Upper bounds must be greater than lower bounds.');

            obj.initialize(n, a, b);
        end

        function Y = collocate(obj, X, I)
            % COLLOCATE Map reference points to physical coordinates.
            %
            %   Y = collocate(obj, X) maps reference points X from the
            %   reference element to physical grid coordinates across
            %   all elements.
            %
            %   Y = collocate(obj, X, I) maps points to specific elements
            %   identified by multi-indices I.
            %
            % Inputs:
            %   obj - The UniformGrid object
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
                X = num2cell(reshape(X, d, []), 2).';
            end

            Y = cell(d, 1);
            for i = 1:d
                h = obj.spacings{i};
                if isempty(I)
                    a = obj.centroids{i};
                    b = bsxfun(@times, h, X{i}(:));
                else
                    a = obj.centroids{i}(I(:, i));
                    b = bsxfun(@times, h, X{i}(:));
                end
                Y{i} = reshape(bsxfun(@plus, a, b), 1, []);
            end

            if ~isempty(I), Y = cell2mat(Y); end
        end

        function newObj = refine(obj, k)
            % REFINE Refine the uniform grid.
            %
            %   newObj = refine(obj, k) refines the grid by subdividing
            %   each element into k subelements per dimension, maintaining
            %   uniformity.
            %
            % Inputs:
            %   obj - The UniformGrid object
            %   k - Refinement factor per dimension (vector or scalar)
            %
            % Outputs:
            %   newObj - The refined UniformGrid object

            d = obj.nDims;
            if isscalar(k)
                k = repmat(k, 1, d);
            end
            
            core.except.assert(length(k) == d, ...
                'DimensionMismatch', ...
                'Refinement factor must match grid dimensions.');
            
            n = obj.resolution;
            newObj = approx.mesh.UniformGrid(n .* k, obj.bbox);
        end
    end

    methods (Access = protected)
        function obj = initialize(obj, n, a, b)
            % INITIALIZE Initialize the uniform grid structure.
            %
            %   obj = initialize(obj, n, a, b) sets up the grid with
            %   specified resolution and bounds, computing centroids and
            %   spacings.
            %
            % Inputs:
            %   obj - The UniformGrid object
            %   n - Number of elements per dimension
            %   a - Lower bounds per dimension
            %   b - Upper bounds per dimension
            %
            % Outputs:
            %   obj - The UniformGrid object

            d = obj.nDims;
            V = arrayfun(@(i) linspace(a(i), b(i), n(i)+1), 1:d, 'Un', 0);
            c = arrayfun(@(i) (V{i}(1:end - 1) + V{i}(2:end))/2, 1:d, 'Un', 0);
            h = num2cell((b - a)./n);
            obj.resolution = n;
            obj.centroids = c;
            obj.spacings = h;
            obj.setElements();
            obj.setBoundary();
            obj.setIndexer();
        end

        function h = getMinSpacing(obj)
            % GETMINSPACING Get the minimum element spacing.
            %
            %   h = getMinSpacing(obj) returns the minimum spacing among
            %   all dimensions in the uniform grid.
            %
            % Inputs:
            %   obj - The UniformGrid object
            %
            % Outputs:
            %   h - Minimum element spacing (positive scalar)

            h = min([obj.spacings{:}]);
        end

        function h = getAllElementMagnitudes(obj)
            % GETALLELEMENTMAGNITUDES Get element magnitudes.
            %
            %   h = getAllElementMagnitudes(obj) computes the magnitude
            %   (volume/area/length) for all elements. For uniform grids,
            %   all elements have the same magnitude.
            %
            % Inputs:
            %   obj - The UniformGrid object
            %
            % Outputs:
            %   h - Scalar element magnitude (all elements are identical)

            h = prod(cell2mat(obj.spacings));
        end

        function G = getSeparableGraph(obj)
            % GETSEPARABLEGRAPH Get the underlying separable graph.
            %
            %   G = getSeparableGraph(obj) creates a separable graph
            %   representation of the uniform grid connectivity.
            %
            % Inputs:
            %   obj - The UniformGrid object
            %
            % Outputs:
            %   G - SeparableGraph representing grid connectivity

            h = obj.spacings;
            c = obj.centroids;
            V = arrayfun(@(i) (c{i}(1) - h{i} / 2):h{i}:(c{i}(end) + h{i} / 2), 1:obj.nDims, 'Un', 0);
            E = arrayfun(@(n) [(1:n).', (2:(n + 1)).'], obj.resolution, 'Un', 0);
            G = approx.mesh.SeparableGraph(V, E);
        end
    end
end