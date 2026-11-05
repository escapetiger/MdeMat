classdef Grid < approx.mesh.Mesh
    % GRID Base class for multidimensional structured grids.
    %
    %   Grid provides a foundational framework for representing structured
    %   computational grids using memory-efficient tensor-based storage.
    %   Unlike unstructured meshes that store explicit connectivity,
    %   structured grids leverage their regular topology to achieve O(d*n)
    %   memory usage instead of O(n^d) for d-dimensional grids with n
    %   elements per dimension.
    %
    % See also:
    %   approx.mesh.Mesh, approx.mesh.SeparableGraph

    properties (Access = public)
        resolution  % Element count per dimension (1 × nDims vector)
        centroids   % Element center coordinates per dimension (1 × nDims cell)
        spacings    % Element spacing per dimension (1 × nDims cell)
        elements    % Element index ranges per dimension (1 × nDims cell)
        boundary    % Boundary element specifications (1 × 2nDims cell)
        indexer     % Multi-index converter (core.linalg.MultiIndexer)
    end

    properties (Dependent)
        nDims       % Number of spatial dimensions (positive integer)
        nElements   % Total number of elements in the grid (positive integer)
    end

    methods
        function obj = Grid(domain, reference, resolution, centroids, spacings)
            % GRID Constructor for Grid base class.
            %
            %   obj = Grid(domain, reference, resolution, centroids,
            %   spacings) creates a structured grid with the specified
            %   geometric domain, reference element, resolution, and
            %   spacing information.
            %
            % Inputs:
            %   domain - Physical domain geometry (core.geometry.Geometry)
            %   reference - Reference element geometry (core.geometry.Geometry)
            %   resolution - Elements per dimension (1 × nDims vector)
            %   centroids - Element centers per dimension (1 × nDims cell)
            %   spacings - Element spacings per dimension (1 × nDims cell)
            %
            % Outputs:
            %   obj - Constructed Grid object

            obj@approx.mesh.Mesh(domain, reference);
            
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
            % GET.NDIMS Get the number of spatial dimensions.

            n = length(obj.resolution);
        end

        function n = get.nElements(obj)
            % GET.NELEMENTS Get the total number of elements.

            n = prod(obj.resolution);
        end

        function IE = findInteriorElements(obj)
            % FINDINTERIORELEMENTS Find linear indices of interior elements.
            %
            %   IE = findInteriorElements(obj) computes linear indices of
            %   interior elements.
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
            % FINDBOUNDARYELEMENTS Find linear indices of boundary elements.
            %
            %   IE = findBoundaryElements(obj, IFE) computes linear indices
            %   of elements that share the specified boundary face,
            %   essential for applying boundary conditions in numerical
            %   methods.
            %
            % Inputs:
            %   obj - The Grid object
            %   IFE - Boundary face index (1 to 2×nDims)
            %
            % Outputs:
            %   IE - Linear indices of boundary elements (vector)

            IF = obj.boundary{IFE};
            ME = obj.indexer.factorToMulti(IF);
            IE = obj.indexer.multiToLinear(ME);
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
            %   IFE - Local face index (1 to 2*nDims)
            %   BC - Boundary condition ('periodic' or 'dirichlet')
            %
            % Outputs:
            %   JE - Linear indices of neighbor elements (vector)
            %   JFE - Local face indices of neighbors (vector)

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
end