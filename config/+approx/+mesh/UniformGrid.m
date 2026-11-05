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
    %   grid2D = UniformGrid([0, 0, 1, 1], [10, 10]);
    %
    %   % 3D uniform grid: 5×8×12 elements on custom domain
    %   grid3D = UniformGrid([-1, 0, -2, 1, 2, 3], [5, 8, 12]);
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
        function obj = UniformGrid(bbox, resolution)
            % UNIFORMGRID Constructor for UniformGrid.
            %
            %   obj = UniformGrid(bbox, resolution) creates a uniform grid
            %   with the specified domain bounds @a bbox and grid
            %   resolution @a resolution.
            %
            % Inputs:
            %   bbox - Bounding box coordinates [a1, b1, a2, b2, ...]
            %   resolution - Grid resolution (d x 1 vector)
            %
            % Outputs:
            %   obj - Constructed UniformGrid object.

            bbox = bbox(:).';
            d = length(bbox) / 2;
            a = bbox(1:2:end);
            b = bbox(2:2:end);
            core.except.assert(all(b > a), ...
                'InvalidInput', ...
                'Upper bounds must be greater than lower bounds.');

            nodes = arrayfun(@(i) linspace(a(i), b(i), resolution(i)+1), 1:d, 'Un', 0);
            centroids = cellfun(@(x) (x(1:end-1) + x(2:end))/2, nodes, 'Un', 0);
            spacings = num2cell((b - a)./resolution);
            obj@approx.mesh.Grid(bbox, resolution, centroids, spacings);
        end

        function graphObj = graphify(obj)
            % GRAPHIFY Convert to a separable graph.
            %
            %   graphObj = graphify(obj) creates a separable graph
            %   representation of the uniform grid connectivity.
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
            V = arrayfun(@(i) (c{i}(1) - h{i} / 2):h{i}:(c{i}(end) + h{i} / 2), 1:d, 'Un', 0);
            E = arrayfun(@(k) [(1:k).', (2:(k + 1)).'], n, 'Un', 0);
            graphObj = approx.mesh.SeparableGraph(V, E);
        end

        function Y = collocate(obj, X, L)
            % COLLOCATE Map reference points to physical coordinates.
            %
            %   Y = collocate(obj, X) maps reference points X from the
            %   reference element to physical grid coordinates across
            %   all elements.
            %
            %   Y = collocate(obj, X, L) maps points to specific elements
            %   identified by linear indices L.
            %
            % Inputs:
            %   obj - The UniformGrid object
            %   X - Reference points (cell array or matrix)
            %   L - Element linear indices (optional, vector)
            %
            % Outputs:
            %   Y - Physical coordinates (cell array or matrix)

            if nargin < 3, L = []; end

            d = obj.nDims;

            if ismatrix(X) && ~iscell(X)
                X = num2cell(reshape(X, d, []), 2).';
            end

            Y = cell(d, 1);
            if isempty(L)
                for i = 1:d
                    h = obj.spacings{i};
                    a = obj.centroids{i};
                    b = bsxfun(@times, h, X{i}(:));
                    Y{i} = reshape(bsxfun(@plus, a, b), 1, []);
                end
            else
                M = obj.indexer.linearToMulti(L);
                for i = 1:d
                    h = obj.spacings{i};
                    a = obj.centroids{i}(M(:, i));
                    b = bsxfun(@times, h, X{i}(:));
                    Y{i} = reshape(bsxfun(@plus, a, b), 1, []);
                end
            end

            if ~isempty(L), Y = cell2mat(Y); end
        end

        function newObj = refine(obj, nLevels)
            % REFINE Refine the uniform grid.
            %
            %   newObj = refine(obj, nLevels) creates a refined uniform
            %   grid by subdividing each element according to the specified
            %   refinement level.
            %
            % Inputs:
            %   obj - The UniformGrid object
            %   nLevels - Number of refinement levels (scalar)
            %
            % Outputs:
            %   newObj - The refined UniformGrid object

            if nLevels == 0
                newObj = obj;
                return;
            end

            n = obj.resolution;
            p = 2^nLevels;
            newObj = approx.mesh.UniformGrid(obj.bbox, n*p);
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

            h = min([obj.spacings{:}]);
        end
        
        function detJ = computeElementJacobianDeterminants(obj)
            % COMPUTEELEMENTJACOBIANDETERMINANTS Compute element Jacobian
            % determinants.
            %
            %   detJ = computeElementJacobianDeterminants(obj) computes the
            %   Jacobian determinant for coordinate transformation from
            %   reference element to physical elements. For uniform grids,
            %   all elements share the same Jacobian determinant.
            %
            % Inputs:
            %   obj - The UniformGrid object
            %
            % Outputs:
            %   detJ - Shared Jacobian determinant (scalar)
            %
            % Notes:
            %   For uniform grids with spacing h = [h1, h2, ..., hd]:
            %   detJ = h1 * h2 * ... * hd
            
            detJ = prod(cell2mat(obj.spacings));
        end
        
        function J = computeElementJacobians(obj)
            % COMPUTEELEMENTJACOBIANS Compute element Jacobian matrices.
            %
            %   J = computeElementIJacobians(obj) computes the Jacobian
            %   matrix for coordinate transformation. For uniform grids,
            %   all elements share the same diagonal Jacobian matrix.
            %
            % Inputs:
            %   obj - The UniformGrid object
            %
            % Outputs:
            %   J - Shared Jacobian matrix (nDims × nDims)
            %
            % Notes:
            %   For uniform grids with spacing h = [h1, h2, ..., hd]:
            %   J = diag([h1, h2, ..., hd])
            
            J = diag(cell2mat(obj.spacings));
        end

        function invJ = computeElementInverseJacobians(obj)
            % COMPUTEELEMENTINVERSEJACOBIANS Compute element inverse
            % Jacobian matrices.
            %
            %   invJ = computeElementInverseJacobians(obj) computes the
            %   inverse Jacobian matrix for coordinate transformation.
            %   For uniform grids, all elements share the same diagonal
            %   inverse Jacobian matrix.
            %
            % Inputs:
            %   obj - The UniformGrid object
            %
            % Outputs:
            %   invJ - Shared inverse Jacobian matrix (nDims × nDims)
            %
            % Notes:
            %   For uniform grids with spacing h = [h1, h2, ..., hd]:
            %   invJ = diag([1/h1, 1/h2, ..., 1/hd])

            invJ = diag(1 ./ cell2mat(obj.spacings));
        end
        
        function detJFace = computeFaceJacobianDeterminants(obj, i)
            % COMPUTEFACEJACOBIANDETERMINANTS Compute face Jacobian
            % determinants.
            %
            %   detJFace = computeFaceJacobianDeterminants(obj, i) computes
            %   the Jacobian determinant for coordinate transformation on
            %   the specified face. For uniform grids, all faces of the
            %   same type share the same determinant.
            %
            % Inputs:
            %   obj - The UniformGrid object
            %   i - Face index (1 to 2*nDims)
            %
            % Outputs:
            %   detJFace - Shared face Jacobian determinant (scalar)
            %
            % Notes:
            %   For face perpendicular to dimension d with spacing hd:
            %   detJFace = product of spacings in other dimensions
            
            core.except.assert(i >= 1 && i <= 2*obj.nDims, ...
                'InvalidInput', 'Face index must be between 1 and 2*nDims.');
            
            %< Determine which dimension is perpendicular to this face
            dim = ceil(i / 2);
            
            %< Face Jacobian determinant is product of spacings in other dimensions
            spacings = cell2mat(obj.spacings);
            otherDims = 1:obj.nDims;
            otherDims(dim) = [];
            
            if isempty(otherDims)
                %< 1D case: face is a point, determinant is 1
                detJFace = 1;
            else
                detJFace = prod(spacings(otherDims));
            end
        end
        
        function normals = computeOutwardNormals(obj, i)
            % COMPUTEOUTWARDNORMALS Compute outward normal vectors.
            %
            %   normals = computeOutwardNormals(obj, i) computes
            %   the outward unit normal vector for the specified face.
            %   For uniform grids, the normal is constant and axis-aligned.
            %
            % Inputs:
            %   obj - The UniformGrid object
            %   i - Face index (1 to 2*nDims)
            %
            % Outputs:
            %   normals - Outward unit normal vector (nDims × 1)
            %
            % Notes:
            %   For axis-aligned faces:
            %   - Face 2*d-1: normal = -e_d (negative d-direction)
            %   - Face 2*d:   normal = +e_d (positive d-direction)
            %   where e_d is the unit vector in dimension d
            
            core.except.assert(i >= 1 && i <= 2*obj.nDims, ...
                'InvalidInput', 'Face index must be between 1 and 2*nDims.');
            
            %< Determine which dimension and direction
            dim = ceil(i / 2);
            isPositive = (mod(i, 2) == 0);
            
            %< Create unit normal vector
            normals = zeros(obj.nDims, 1);
            if isPositive
                normals(dim) = 1;   %< Positive direction
            else
                normals(dim) = -1;  %< Negative direction
            end
        end
    end
end