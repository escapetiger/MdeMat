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

        function IE = findInteriorElements(obj)
            % FINDINTERIORELEMENTS Find linear indices of interior
            % elements.
            %
            %   IE = findInteriorElements(obj) computes the linear indices
            %   of interior elements.
            %
            % Inputs:
            %   obj - The Grid object
            %
            % Outputs:
            %   IE - Linear indices of interior elements (vector)

            ME = cellfun(@(I) I(2:end-1), obj.elements, 'Un', 0);
            ME = obj.indexer.factorToMulti(ME);
            IE = obj.indexer.multiToLinear(ME);
        end

        function IE = findBoundaryElements(obj, IFE)
            % FINDBOUNDARYELEMENTS Find linear indices of boundary
            % elements.
            %
            %   IE = findBoundaryElements(obj, IFE) computes the linear
            %   indices of elements sharing the @a IFE-th face with
            %   boundary.
            %
            % Inputs:
            %   obj - The Grid object
            %   IFE - Face index (scalar)
            %
            % Outputs:
            %   IE - Linear indices of boundary elememts (vector)

            IF = obj.boundary{IFE};
            ME = obj.indexer.factorToMulti(IF);
            IE = obj.indexer.multiToLinear(ME);
        end
    
        function [JE, JFE] = findNeighborElements(obj, IE, IFE, BC)
            % FINDNEIGHBORELEMENTS Find linear indices of neighbor
            % elements.
            %
            %   [JE, JFE] = findNeighborElements(obj, IE, IFE, BC) computes
            %   the linear indices of neighbor elements sharing the @a
            %   IFE-th face of prime elements @a IE with the specified
            %   boundary condition @a BC.
            %
            % Inputs:
            %   obj - The Grid object
            %   IE - Linear indices of prime elements (vector)
            %   IFE - Local face index of prime elements (scalar)
            %   BC - Boundary condition ('periodic' or 'dirichlet')
            %
            % Outputs:
            %   JE - Linear indices of neighbor elements (vector)
            %   JFE - Local face index of neighbor elements (vector) 

            ME = obj.indexer.linearToMulti(IE);
            dim = ceil(IFE / 2);
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
        detJ = computeElementJacobianDeterminants(obj, IE)
        
        % COMPUTEELEMENTJACOBIANS Compute element Jacobian matrices.
        J = computeElementJacobians(obj, IE)

        % COMPUTEELEMENTINVERSEJACOBIANS Compute element inverse Jacobian
        % matrices.
        invJ = computeElementInverseJacobians(obj, IE)
        
        % COMPUTEFACEJACOBIANDETERMINANTS Compute face Jacobian
        % determinants.
        detJFace = computeFaceJacobianDeterminants(obj, IE, IFE)
        
        % COMPUTEOUTWARDNORMALS Compute outward normal vectors.
        normals = computeOutwardNormals(obj, IE, IFE)
    end
end