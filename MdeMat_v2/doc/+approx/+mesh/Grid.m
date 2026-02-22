classdef Grid < handle
    % GRID Abstract base class for multidimensional structured grids.
    %
    %   Grid provides a foundation for representing multidimensional
    %   structured grids using memory-efficient tensor-based storage.
    %   This abstract class defines common properties and methods for
    %   querying grid characteristics, connectivity, and performing
    %   grid-based operations.
    %
    % See also:
    %   approx.mesh.UniformGrid, approx.mesh.NonuniformGrid,
    %   approx.mesh.SeparableGraph

    properties (Access = public)
        bbox % Bounding box coordinates [a1, b1, a2, b2, ...]
        resolution % Vector of grid dimensions per axis
        centroids % Cell array of element centroids per dimension
        spacings % Cell array of element spacings per dimension
        elements % Cell array of element indices per dimension
        boundary % Cell array of boundary element linear indices
        indexer % Multi-indexer for linear/multi-index conversion
    end

    properties (Dependent)
        nDims % Number of spatial dimensions (positive integer)
        nElements % Total number of elements in the grid
    end

    methods
        function obj = Grid(bbox, resolution, centroids, spacings)
            % GRID Constructor for Grid.
            %
            %   obj = Grid(bbox, resolution, centroids, spacings) creates a
            %   grid with the specified domain bounds @a bbox, grid
            %   resolution @a resolution, center points @a centroids and
            %   spacings @a spacings.
            %
            % Inputs:
            %   bbox - Bounding box coordinates [a1, b1, a2, b2, ...]
            %   resolution - Grid resolution (d x 1 vector)
            %   centroids - Center points (d x 1 cell)
            %   spacings - Spacings (d x 1 cell)
            %
            % Outputs:
            %   obj - Constructed UniformGrid object.

            obj.bbox = bbox;
            obj.resolution = resolution;
            obj.centroids = centroids;
            obj.spacings = spacings;
            obj.indexer = core.linalg.MultiIndexer(obj.resolution);
            obj.elements = arrayfun(@(n) 1:n, obj.resolution, 'Un', 0);
            obj.boundary = cell(1, 2*obj.nDims);
            for i = 1:obj.nDims
                f = arrayfun(@(k) 1:k, obj.resolution, 'Un', 0);
                f{i} = 1;
                obj.boundary{2*i-1} = f;
                f{i} = obj.resolution(i);
                obj.boundary{2*i} = f;
            end
        end

        function n = get.nDims(obj)
            % GET.NDIMS Get the number of dimensions.

            n = length(obj.resolution);
        end

        function n = get.nElements(obj)
            % GET.NELEMENTS Get the total number of elements.

            n = prod(obj.resolution);
        end

        function M = linearToMulti(obj, L)
            % LINEARTOMULTI Convert linear indices to multi-indices.
            %
            %   M = linearToMulti(obj, L) converts linear indices @a L to
            %   multi-indices for tensor with specified shape.
            %
            % Inputs:
            %   obj - The Grid object
            %   L - Linear indices (vector)
            %
            % Outputs:
            %   M - Multi-indices (matrix)

            M = obj.indexer.linearToMulti(L);
        end

        function L = multiToLinear(obj, M, bc)
            % MULTITOLINEAR Convert multi-indices to linear indices.
            %
            %   L = multiToLinear(obj, M) converts multi-indices @a M to
            %   linear indices based on tensor dimensions.
            %
            %   L = multiToLinear(obj, M, bc) uses boundary condition @a
            %   bc.
            %
            % Inputs:
            %   obj - The Grid object
            %   M - Multi-indices (matrix or cell array)
            %   bc - Boundary condition (optional, default: 1)
            %        * 0: periodic boundary (wrap-around)
            %        * 1: strict boundary (default, out-of-range set to 0)
            %
            % Outputs:
            %   L - Linear indices (vector)

            if nargin < 3 || isempty(bc)
                bc = 1;
            end

            L = obj.indexer.multiToLinear(M, bc);
        end
    
        function L = findInteriorElements(obj)
            % FINDINTERIORELEMENTS Find linear indices of interior
            % elements.
            %
            %   L = findInteriorElements(obj) computes the linear indices
            %   of interior elements.
            %
            % Inputs:
            %   obj - The Grid object
            %
            % Outputs:
            %   L - Linear indices of interior elements (vector)

            F = cellfun(@(I) I(2:end-1), obj.elements, 'Un', 0);
            M = obj.indexer.factorToMulti(F);
            L = obj.indexer.multiToLinear(M);
        end

        function L = findBoundaryElements(obj, i)
            % FINDBOUNDARYELEMENTS Find linear indices of boundary
            % elements.
            %
            %   L = findBoundaryElements(obj, i) computes the linear
            %   indices of elements sharing the i-th face with boundary.
            %
            % Inputs:
            %   obj - The Grid object
            %   i - Face index (scalar)
            %
            % Outputs:
            %   L - Linear indices (vector)

            F = obj.boundary{i};
            M = obj.indexer.factorToMulti(F);
            L = obj.indexer.multiToLinear(M);
        end
    
        function L = findNeighborElements(obj, i, K, bc)
            % FINDNEIGHBORELEMENTS Find linear indices of neighbor
            % elements.
            %
            %   L = findNeighborElements(obj, i) computes the linear
            %   indices of neighbor elements sharing the i-th face of prime
            %   elements.
            %
            % Inputs:
            %   obj - The Grid object
            %   i - Face index (scalar)
            %   K - Linear indices of prime elements (vector)
            %   bc - Boundary condition ('periodic' or 'strict')
            %
            % Outputs:
            %   L - Linear indices (vector)

            if nargin < 4 || isempty(bc)
                bc = 'strict';
            end

            M = obj.indexer.linearToMulti(K);
            d = ceil(i / 2);
            if mod(i, 2) == 1
                M(:, d) = M(:, d) - 1;
            else
                M(:, d) = M(:, d) + 1;
            end
            
            switch lower(bc)
                case 'periodic'
                    L = obj.indexer.multiToLinear(M, 0);
                case 'strict'
                    L = obj.indexer.multiToLinear(M, 1);
            end
        end
    end

    methods (Abstract)
        % GRAPHIFY Build the separable graph.
        graphObj = graphify(obj)

        % COLLOCATE Map reference nodes to physical nodes.
        Y = collocate(obj, X, L)

        % REFINE Refine the grid.
        newObj = refine(obj, nLevels)

        % COMPUTEMEASURE Compute the mesh measure.
        h = computeMeasure(obj)
        
        % COMPUTEELEMENTJACOBIANDETERMINANTS Compute element Jacobian
        % determinants.
        detJ = computeElementJacobianDeterminants(obj)
        
        % COMPUTEELEMENTJACOBIANS Compute element Jacobian matrices.
        J = computeElementJacobians(obj)

        % COMPUTEELEMENTINVERSEJACOBIANS Compute element inverse Jacobian
        % matrices.
        invJ = computeElementInverseJacobians(obj)
        
        % COMPUTEFACEJACOBIANDETERMINANTS Compute face Jacobian
        % determinants.
        detJFace = computeFaceJacobianDeterminants(obj, faceIndex)
        
        % COMPUTEOUTWARDNORMALS Compute outward normal vectors.
        normals = computeOutwardNormals(obj, faceIndex)
    end
end