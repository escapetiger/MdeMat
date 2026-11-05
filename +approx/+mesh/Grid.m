classdef Grid < approx.mesh.Mesh
    % GRID Base class for multidimensional structured grids.
    %
    %   Grid provides a foundational framework for representing structured
    %   computational grids using memory-efficient tensor-based storage.
    %   Unlike unstructured meshes that store explicit connectivity,
    %   structured grids leverage their regular topology to achieve O(d×n)
    %   memory usage instead of O(n^d) for d-dimensional grids with n
    %   elements per dimension.
    %
    % See also:
    %   approx.mesh.Mesh, approx.mesh.SeparableGraph

    properties (Access = public)
        resolution % Element count per dimension (1 × nDims vector)
        centroids % Element center coordinates per dimension (1 × nDims cell)
        spacings % Element spacing per dimension (1 × nDims cell)
        boundary % Boundary element specifications (1 × 2nDims cell)
        indexer % Multi-index converter (core.linalg.MultiIndexer)
    end

    properties (Dependent)
        bbox % Bounding box
    end

    methods
        function obj = Grid(vertices)
            % GRID Constructor for Grid base class.
            %
            %   obj = Grid(vertices) creates a non-uniform grid from vertex
            %   coordinate vectors that define element boundaries in each
            %   dimension, automatically computing element properties.
            %
            % Inputs:
            %   vertices - Vertex coordinate cell array (1 × nDims cell)
            %
            % Outputs:
            %   obj - Constructed Grid object

            elements = cellfun(@(x) 1:(length(x) - 1), vertices, 'Un', 0);
            obj@approx.mesh.Mesh(vertices, elements);
            obj.resolution = cellfun(@(x) length(x)-1, vertices);
            obj.centroids = cellfun(@(x) (x(1:end-1) + x(2:end))/2, vertices, 'Un', 0);
            obj.spacings = cellfun(@(x) diff(x), vertices, 'Un', 0);
            obj.indexer = core.linalg.MultiIndexer(obj.resolution);
            obj.boundary = cell(1, 2*obj.nDims);
            for i = 1:obj.nDims
                f = arrayfun(@(n) 1:n, obj.resolution, 'Un', 0);
                f{i} = 1;
                obj.boundary{2*i-1} = f;
                f{i} = obj.resolution(i);
                obj.boundary{2*i} = f;
            end
        end

        function IE = getInternalElements(obj)
            % GETINTERIORELEMENTS Get linear indices of internal elements.
            %
            %   IE = getInternalElements(obj) computes linear indices of
            %   internal elements that do not touch any boundary.
            %
            % Inputs:
            %   obj - The Grid object
            %
            % Outputs:
            %   IE - Linear indices of internal elements (vector)

            ME = cellfun(@(I) I(2:end-1), obj.elements, 'Un', 0);
            ME = obj.indexer.factorToMulti(ME);
            IE = obj.indexer.multiToLinear(ME);
        end

        function IE = getBoundaryElements(obj, IFE)
            % GETBOUNDARYELEMENTS Get linear indices of boundary elements.
            %
            %   IE = getBoundaryElements(obj) computes the linear indices
            %   of all elements that have at least one face on the
            %   boundary.
            %
            %   IE = getBoundaryElements(obj, IFE) finds elements whose
            %   IFE-th face is on the boundary.
            %
            % Inputs:
            %   obj - The Grid object
            %   IFE - Boundary face index (1 to 2×nDims)
            %
            % Outputs:
            %   IE - Linear indices of boundary elements (vector)

            if nargin < 2 || isempty(IFE)
                IE = cellfun(@(IF) obj.indexer.multiToLinear(IF), obj.boundary, 'Un', 0);
                IE = vertcat(IE{:});
                IE = unique(IE, 'stable');
                return;
            end

            ME = obj.indexer.factorToMulti(obj.boundary{IFE});
            IE = obj.indexer.multiToLinear(ME);
        end

        function JE = findInternalElements(obj, IE, IFE)
            % FINDINTERIORELEMENTS Find internal elements across faces.
            %
            %   JE = findInternalElements(obj, IE, IFE) finds elements
            %   whose IFE-th face is not on the boundary within the
            %   specified elements @a IE.
            %
            % Inputs:
            %   obj - The Grid object
            %   IE - Linear indices of primary elements (vector)
            %   IFE - Local face index (1 to 2×nDims)
            %
            % Outputs:
            %   JE - Linear indices of internal elements (vector)

            ME = obj.indexer.linearToMulti(IE);
            dim = ceil(IFE/2);
            if mod(IFE, 2) == 1
                mask = ME(:, dim) > 1;
            else
                mask = ME(:, dim) < obj.resolution(dim);
            end
            JE = obj.indexer.multiToLinear(ME(mask, :));
        end

        function JE = findBoundaryElements(obj, IE, IFE)
            % FINDBOUNDARYELEMENTS Find boundary elements across faces.
            %
            %   JE = findBoundaryElements(obj, IE, IFE) finds elements
            %   whose IFE-th face is on the boundary within the specified
            %   elements @a IE.
            %
            % Inputs:
            %   obj - The Grid object
            %   IE - Linear indices of primary elements (vector)
            %   IFE - Local face index (1 to 2×nDims)
            %
            % Outputs:
            %   JE - Linear indices of boundary elements (vector)

            ME = obj.indexer.linearToMulti(IE);
            dim = ceil(IFE/2);
            if mod(IFE, 2) == 1
                mask = ME(:, dim) == 1;
            else
                mask = ME(:, dim) == obj.resolution(dim);
            end
            JE = obj.indexer.multiToLinear(ME(mask, :));
        end

        function [JE, JFE] = findNeighborElements(obj, IE, IFE, BC)
            % FINDNEIGHBORELEMENTS Find neighboring elements across faces.
            %
            %   [JE, JFE] = findNeighborElements(obj, IE, IFE, BC) computes
            %   neighbor elements sharing the specified local face of given
            %   elements, with appropriate boundary condition treatment.
            %
            % Inputs:
            %   obj - The Grid object
            %   IE - Linear indices of primary elements (vector)
            %   IFE - Local face index (1 to 2×nDims)
            %   BC - Boundary condition ('periodic' or 'dirichlet' or [])
            %
            % Outputs:
            %   JE - Linear indices of neighbor elements (vector)
            %   JFE - Local face indices of neighbors (vector)

            ME = obj.indexer.linearToMulti(IE);

            dim = ceil(IFE/2);

            if mod(IFE, 2) == 1
                ME(:, dim) = ME(:, dim) - 1;
                JFE = 2 * dim;
            else
                ME(:, dim) = ME(:, dim) + 1;
                JFE = 2 * dim - 1;
            end

            switch lower(BC)
                case 'periodic'
                    JE = obj.indexer.multiToLinear(ME, 'wrap');
                case 'dirichlet'
                    JE = obj.indexer.multiToLinear(ME, 'empty');
                otherwise
                    JE = obj.indexer.multiToLinear(ME, 'empty');
            end
        end

        function graphObj = graphify(obj)
            % GRAPHIFY Convert grid to separable graph representation.
            %
            %   graphObj = graphify(obj) creates a SeparableGraph that
            %   represents the grid's connectivity structure using 1D
            %   factors with variable vertex spacing.
            %
            % Inputs:
            %   obj - The Grid object
            %
            % Outputs:
            %   graphObj - SeparableGraph with tensor-product structure

            d = obj.nDims;
            n = obj.resolution;
            h = obj.spacings;
            c = obj.centroids;
            V = arrayfun(@(i) [c{i} - h{i} / 2, c{i}(end) + h{i}(end) / 2].', 1:d, 'Un', 0);
            E = arrayfun(@(k) [(1:k).', (2:(k + 1)).'], n, 'Un', 0);
            graphObj = approx.mesh.SeparableGraph(V, E);
        end

        function Y = collocate(obj, X, IE)
            % COLLOCATE Map reference coordinates to physical coordinates.
            %
            %   Y = collocate(obj, X) maps reference coordinates from the
            %   standard hypercube $[-1/2,1/2]^d$ to physical coordinates
            %   across all grid elements using element-specific
            %   transformations.
            %
            %   Y = collocate(obj, X, IE) maps coordinates to specific
            %   elements identified by linear indices IE.
            %
            % Inputs:
            %   obj - The Grid object
            %   X - Reference coordinates (d × nPoints or cell array)
            %   IE - Element linear indices (optional, vector)
            %
            % Outputs:
            %   Y - Physical coordinates (d × nPoints or cell array)

            if nargin < 3, IE = []; end

            d = obj.nDims;
            if ~iscell(X)
                C = cell(1, d);
                [C{:}] = ndgrid(obj.centroids{:});
                C = reshape(cat(d+1, C{:}), [], d).';
                if ~isempty(IE), C = C(:, IE); end
                C = reshape(C, d, 1, []);

                h = cell(1, d);
                [h{:}] = ndgrid(obj.spacings{:});
                h = reshape(cat(d+1, h{:}), [], d).';
                if ~isempty(IE), h = h(:, IE); end
                h = reshape(h, d, 1, []);
                
                X = reshape(X, d, [], 1);
                Y = C + h .* X;
                Y = reshape(Y, d, []);
            else
                Y = cell(d, 1);
                if isempty(IE)
                    for i = 1:d
                        C = obj.centroids{i};
                        h = obj.spacings{i};
                        Y{i} = C(:).' + h(:).' .* X{i}(:).';
                    end
                else
                    M = obj.indexer.linearToMulti(IE);
                    for i = 1:d
                        C = obj.centroids{i}(M(:, i));
                        h = obj.spacings{i}(M(:, i));
                        Y{i} = C(:).' + h(:).' .* X{i}(:).';
                    end
                end
            end
        end

        function newObj = refine(obj, nLevels)
            % REFINE Create refined grid preserving spacing patterns.
            %
            %   newObj = refine(obj, nLevels) creates a refined grid by
            %   uniformly subdividing each element while preserving the
            %   relative spacing characteristics of the original grid.
            %
            % Inputs:
            %   obj - The Grid object
            %   nLevels - Number of refinement levels (non-negative integer)
            %
            % Outputs:
            %   newObj - Refined Grid object

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
            newObj = approx.mesh.Grid(x);
        end

        function h = computeMeasure(obj)
            % COMPUTEMEASURE Compute minimum element spacing across the grid.
            %
            %   h = computeMeasure(obj) returns the minimum element spacing
            %   among all elements in all dimensions, providing a
            %   characteristic length scale for the grid.
            %
            % Inputs:
            %   obj - The Grid object
            %
            % Outputs:
            %   h - Minimum element spacing (positive scalar)

            h = min(cellfun(@min, obj.spacings));
        end

        function detJ = computeAllElementJacobianDeterminants(obj)
            % COMPUTEALLELEMENTJACOBIANDETERMINANTS Compute
            % element-specific Jacobian determinants.
            %
            %   detJ = computeAllElementJacobianDeterminants(obj) computes
            %   the Jacobian determinant for coordinate transformation for
            %   each element. For structured grids, determinants vary by
            %   element based on local spacing.
            %
            % Inputs:
            %   obj - The Grid object
            %
            % Outputs:
            %   detJ - Jacobian determinants for all elements (nElements × 1)

            h = cell(1, obj.nDims);
            [h{:}] = ndgrid(obj.spacings{:});
            h = cellfun(@(x) x(:).', h, 'Un', 0);
            detJ = prod(cat(1, h{:}), 1);
            detJ = detJ(:);
        end

        function J = computeAllElementJacobians(obj)
            % COMPUTEALLELEMENTJACOBIANS Compute element-specific Jacobian
            % matrices.
            %
            %   J = computeAllElementJacobians(obj) computes the Jacobian
            %   matrix for coordinate transformation for each element. For
            %   structured grids, Jacobians vary by element but remain
            %   diagonal.
            %
            % Inputs:
            %   obj - The Grid object
            %
            % Outputs:
            %   J - Jacobian matrices (d × d × nElements)

            d = obj.nDims;
            h = cell(1, d);
            [h{:}] = ndgrid(obj.spacings{:});
            h = cellfun(@(x) x(:).', h, 'Un', 0);
            h = cat(1, h{:});
            J = reshape(h, d, 1, []) .* eye(d, d);
        end

        function invJ = computeAllElementInverseJacobians(obj)
            % COMPUTEALLELEMENTINVERSEJACOBIANS Compute element-specific
            % inverse Jacobian matrices.
            %
            %   invJ = computeAllElementInverseJacobians(obj) computes the
            %   inverse Jacobian matrix for coordinate transformation for
            %   each element. For structured grids, inverse Jacobians vary
            %   by element.
            %
            % Inputs:
            %   obj - The Grid object
            %
            % Outputs:
            %   invJ - Inverse Jacobian matrices (d × d × nElements)

            d = obj.nDims;
            h = cell(1, d);
            [h{:}] = ndgrid(obj.spacings{:});
            h = cellfun(@(x) x(:).', h, 'Un', 0);
            h = 1 ./ cat(1, h{:});
            invJ = reshape(h, d, 1, []) .* eye(d, d);
        end

        function detFJ = computeAllFaceJacobianDeterminants(obj)
            % COMPUTEALLFACEJACOBIANDETERMINANTS Compute element-specific
            % face Jacobian determinants.
            %
            %   detFJ = computeAllFaceJacobianDeterminants(obj) computes
            %   Jacobian determinants for coordinate transformation on
            %   element faces. For structured grids, face determinants vary
            %   by element.
            %
            % Inputs:
            %   obj - The Grid object
            %
            % Outputs:
            %   detFJ - Face Jacobian determinants (cell array by face type)

            d = obj.nDims;
            if d == 1
                detFJ = repmat({ones(obj.nElements, 1)}, 1, 2);
                return;
            end

            detFJ = cell(1, 2*d);
            for i = 1:2 * d
                dim = ceil(i/2);
                h = cell(1, d);
                [h{:}] = ndgrid(obj.spacings{:});
                h = cellfun(@(x) x(:).', h(setdiff(1:d, dim)), 'Un', 0);
                detFJ{i} = prod(cat(1, h{:}), 1).';
            end
        end

        function outN = computeAllOutwardNormals(obj)
            % COMPUTEALLOUTWARDNORMALS Compute outward normal vectors for
            % all faces.
            %
            %   outN = computeAllOutwardNormals(obj) computes unit outward
            %   normal vectors for each face type. For structured grids,
            %   normals remain axis-aligned and constant for faces of the
            %   same type.
            %
            % Inputs:
            %   obj - The Grid object
            %
            % Outputs:
            %   outN - Outward normal vectors (cell array by face type)

            d = obj.nDims;
            outN = cell(1, 2*d);
            for i = 1:(2 * d)
                dim = ceil(i/2);
                outN{i} = zeros(d, 1);
                if mod(i, 2) == 0
                    outN{i}(dim) = +1;
                else
                    outN{i}(dim) = -1;
                end
            end
        end
    
        function v = computeUniqueDirection(obj)
            % COMPUTEUNIQUEDIRECTION Compute a unique direction vector.
            %
            %   v = computeUniqueDirection(obj) computes a unique unit
            %   direction vector to define the "plus" and "minus" side of a
            %   face.
            %
            % Inputs:
            %   obj - The Mesh object
            %
            % Outputs:
            %   v - Direction vector (nDims × 1)

            d = obj.nDims;
            if d == 1
                v = 1;
            else
                v = ones(d, 1) / sqrt(d);
            end
        end
    end

    methods
        function bbox = get.bbox(obj)
            % GET.BBOX Get the bounding box.

            xmin = cellfun(@(x) min(x), obj.vertices);
            xmax = cellfun(@(x) max(x), obj.vertices);
            bbox = reshape([xmin; xmax], 1, []);
        end
    end
end