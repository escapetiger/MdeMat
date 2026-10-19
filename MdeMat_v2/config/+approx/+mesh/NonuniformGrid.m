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
        function obj = NonuniformGrid(nodes)
            % NONUNIFORMGRID Constructor for NonuniformGrid.
            %
            %   obj = NonuniformGrid(nodes) creates a nonuniform grid from
            %   the specified coordinate vectors.
            %
            % Inputs:
            %   nodes - Cell array of coordinate vectors, one per dimension.
            %
            % Outputs:
            %   obj - Constructed NonuniformGrid object

            core.except.assert(iscell(nodes), 'InvalidInput', ...
                'Nodes must be a cell array.');

            resolution = cellfun(@(x) length(x)-1, nodes);
            spacings = cellfun(@(x) diff(x), nodes, 'Un', 0);
            centroids = cellfun(@(x) (x(1:end-1) + x(2:end))/2, nodes, 'Un', 0);
            bbox = cellfun(@(x) [min(x), max(x)], nodes, 'Un', 0);
            bbox = horzcat(bbox{:});
            obj@approx.mesh.Grid(bbox, resolution, centroids, spacings);
        end

        function graphObj = graphify(obj)
            % GRAPHIFY Convert to a separable graph.
            %
            %   graphObj = graphify(obj) creates a separable graph
            %   representation of the nonuniform grid connectivity.
            %
            % Inputs:
            %   obj - The UniformGrid object
            %
            % Outputs:
            %   graphObj - SeparableGraph object

            d = obj.nDims;
            n = obj.resolution;
            h = obj.spacings;
            c = obj.centroids;
            V = arrayfun(@(i) [c{i} - h{i} / 2, c{i}(end) + h{i}(end) / 2].', 1:d, 'Un', 0);
            E = arrayfun(@(n) [(1:n).', (2:(n + 1)).'], n, 'Un', 0);
            graphObj = approx.mesh.SeparableGraph(V, E);
        end

        function Y = collocate(obj, X, L)
            % COLLOCATE Map reference points to physical coordinates.
            %
            %   Y = collocate(obj, X) maps reference points @a X from the
            %   reference element \f$[-1/2,1/2]^d\f$ to physical grid
            %   coordinates.
            %
            %   Y = collocate(obj, X, L) maps points to specific elements
            %   identified by linear indices @a L.
            %
            % Inputs:
            %   obj - The NonuniformGrid object
            %   X - Reference points (cell array or matrix)
            %   L - Element linear indices (optional, matrix)
            %
            % Outputs:
            %   Y - Physical coordinates (cell array or matrix)

            if nargin < 3, L = []; end

            d = obj.nDims;

            if ismatrix(X) && ~iscell(X)
                X = num2cell(reshape(X, d, []), 2);
            end

            Y = cell(d, 1);
            if isempty(L)
                for i = 1:d
                    a = obj.centroids{i};
                    h = obj.spacings{i};
                    b = bsxfun(@times, h, X{i}(:));
                    Y{i} = reshape(bsxfun(@plus, a, b), 1, []);
                end
            else
                M = obj.indexer.linearToMulti(L);
                for i = 1:d
                    a = obj.centroids{i}(M(:, i));
                    h = obj.spacings{i};
                    b = bsxfun(@times, h(M(:, i)), X{i}(:));
                    Y{i} = reshape(bsxfun(@plus, a, b), 1, []);
                end
            end

            if ~isempty(L), Y = cell2mat(Y); end
        end

        function newObj = refine(obj, nLevels)
            % REFINE Refine the nonuniform grid.
            %
            %   newObj = refine(obj, nLevels) creates a refined nonuniform
            %   grid by subdividing each element according to the specified
            %   refinement level.
            %
            % Inputs:
            %   obj - The NonuniformGrid object
            %   nLevels - Number of refinement levels (scalar)
            %
            % Outputs:
            %   newObj - The refined NonuniformGrid object

            if nLevels == 0
                newObj = obj;
                return;
            end

            d = obj.nDims;
            p = 2^nLevels;
            x = cell(1, d);
            for i = 1:d
                y = obj.centroids{i};
                n = obj.resolution(i);
                s = obj.spacings{i};
                a = y - s / 2;
                b = y + s / 2;
                h = (b - a) ./ p;
                x{i} = zeros(1, 1+n*p);
                x{i}(1) = a(1);
                for j = 1:n
                    x{i}(1 + (j - 1) * p + (1:p)) = (a(j) + h(j)):h(j):b(j);
                end
            end
            newObj = approx.mesh.NonuniformGrid(x);
        end

        function h = computeMeasure(obj)
            % COMPUTEMEASURE Get the minimum element spacing.
            %
            %   h = computeMeasure(obj) returns the minimum spacing among
            %   all dimensions in the uniform grid.
            %
            % Inputs:
            %   obj - The UniformGrid object
            %
            % Outputs:
            %   h - Minimum element spacing (positive scalar)

            h = min(cellfun(@min, obj.spacings));
        end

        function detJ = computeElementJacobianDeterminants(obj)
            % COMPUTEELEMENTJACOBIANDETERMINANTS Compute element Jacobian determinants.
            %
            %   detJ = computeElementJacobianDeterminants(obj) computes the
            %   Jacobian determinant for coordinate transformation from
            %   reference element to physical elements. For nonuniform grids,
            %   elements have different Jacobian determinants.
            %
            % Inputs:
            %   obj - The NonuniformGrid object
            %
            % Outputs:
            %   detJ - Jacobian determinants for all elements (nElements × 1)
            %
            % Notes:
            %   For element (i1, i2, ..., id) with spacings h = [h1(i1), h2(i2), ..., hd(id)]:
            %   detJ = h1(i1) * h2(i2) * ... * hd(id)

            h = cell(1, obj.nDims);
            [h{:}] = ndgrid(obj.spacings{:});
            h = cellfun(@(x) x(:).', h, 'Un', 0);
            detJ = prod(cat(1, h{:}), 1);
            detJ = detJ(:);
        end

        function J = computeElementJacobians(obj)
            % COMPUTEELEMENTJACOBIANS Compute element Jacobian matrices.
            %
            %   J = computeElementJacobians(obj) computes the Jacobian
            %   matrix for coordinate transformation. For nonuniform grids,
            %   elements have different diagonal Jacobian matrices.
            %
            % Inputs:
            %   obj - The NonuniformGrid object
            %
            % Outputs:
            %   J - Jacobian matrices (nDims × nDims × nElements)
            %
            % Notes:
            %   For element (i1, i2, ..., id) with spacings h = [h1(i1), h2(i2), ..., hd(id)]:
            %   J(:,:,element) = diag([h1(i1), h2(i2), ..., hd(id)])

            d = obj.nDims;
            h = cell(1, d);
            [h{:}] = ndgrid(obj.spacings{:});
            h = cellfun(@(x) x(:).', h, 'Un', 0);
            h = cat(1, h{:});
            J = reshape(h, d, 1, []) .* eye(d, d);
        end

        function invJ = computeElementInverseJacobians(obj)
            % COMPUTEELEMENTINVERSEJACOBIANS Compute element inverse Jacobian matrices.
            %
            %   invJ = computeElementInverseJacobians(obj) computes the
            %   inverse Jacobian matrix for coordinate transformation.
            %   For nonuniform grids, elements have different diagonal
            %   inverse Jacobian matrices.
            %
            % Inputs:
            %   obj - The NonuniformGrid object
            %
            % Outputs:
            %   invJ - Inverse Jacobian matrices (nDims × nDims × nElements)
            %
            % Notes:
            %   For element (i1, i2, ..., id) with spacings h = [h1(i1), h2(i2), ..., hd(id)]:
            %   invJ(:,:,element) = diag([1/h1(i1), 1/h2(i2), ..., 1/hd(id)])

            d = obj.nDims;
            h = cell(1, d);
            [h{:}] = ndgrid(obj.spacings{:});
            h = cellfun(@(x) x(:).', h, 'Un', 0);
            h = 1 ./ cat(1, h{:});
            invJ = reshape(h, d, 1, []) .* eye(d, d);
        end

        function detJFace = computeFaceJacobianDeterminants(obj, faceIndex)
            % COMPUTEFACEJACOBIANDETERMINANTS Compute face Jacobian determinants.
            %
            %   detJFace = computeFaceJacobianDeterminants(obj, faceIndex)
            %   computes the Jacobian determinant for coordinate
            %   transformation on the specified face for all elements.
            %   For nonuniform grids, elements have different face
            %   determinants.
            %
            % Inputs:
            %   obj - The NonuniformGrid object
            %   faceIndex - Face index (1 to 2*nDims)
            %
            % Outputs:
            %   detJFace - Face Jacobian determinants for all elements (1 × nElements)
            %
            % Notes:
            %   For face perpendicular to dimension d:
            %   detJFace = product of spacings in other dimensions for each element

            core.except.assert(faceIndex >= 1 && faceIndex <= 2*obj.nDims, ...
                'InvalidInput', 'Face index must be between 1 and 2*nDims.');

            %< Determine which dimension is perpendicular to this face
            dim = ceil(faceIndex/2);

            %< Other dimensions (not perpendicular to face)
            otherDims = 1:obj.nDims;
            otherDims(dim) = [];

            if isempty(otherDims)
                %< 1D case: face is a point, determinant is 1 for all elements
                detJFace = ones(1, obj.nElements);
            else
                %< Product of spacings in other dimensions for each element
                h = cell(1, obj.nDims);
                [h{:}] = ndgrid(obj.spacings{:});
                h = cellfun(@(x) x(:).', h(otherDims), 'Un', 0);
                detJFace = prod(cat(1, h{:}), 1);
            end
        end

        function normals = computeOutwardNormals(obj, faceIndex)
            % COMPUTEOUTWARDNORMALS Compute outward normal vectors.
            %
            %   normals = computeOutwardNormals(obj, faceIndex) computes
            %   the outward unit normal vector for the specified face.
            %   For nonuniform grids, the normal is constant and axis-aligned
            %   (same for all boundary elements of the same face).
            %
            % Inputs:
            %   obj - The NonuniformGrid object
            %   faceIndex - Face index (1 to 2*nDims)
            %
            % Outputs:
            %   normals - Outward unit normal vector (nDims × 1)
            %
            % Notes:
            %   For axis-aligned faces:
            %   - Face 2*d-1: normal = -e_d (negative d-direction)
            %   - Face 2*d:   normal = +e_d (positive d-direction)
            %   where e_d is the unit vector in dimension d
            %
            %   Note: All boundary elements on the same face share the same normal

            core.except.assert(faceIndex >= 1 && faceIndex <= 2*obj.nDims, ...
                'InvalidInput', 'Face index must be between 1 and 2*nDims.');

            %< Determine which dimension and direction
            dim = ceil(faceIndex/2);
            isPositive = (mod(faceIndex, 2) == 0);

            %< Create unit normal vector
            normals = zeros(obj.nDims, 1);
            if isPositive
                normals(dim) = 1; %< Positive direction
            else
                normals(dim) = -1; %< Negative direction
            end
        end
    end
end