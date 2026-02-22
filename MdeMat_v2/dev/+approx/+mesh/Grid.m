classdef Grid < handle
    % GRID Abstract base class for multidimensional structured grids.
    %
    %   Grid provides a foundation for representing multidimensional
    %   structured grids using memory-efficient tensor-based storage.
    %   This abstract class defines common properties and methods for
    %   querying grid characteristics, connectivity, and performing
    %   grid-based operations.
    %
    % Examples:
    %   % Cannot instantiate abstract class directly
    %   % Use concrete subclasses instead:
    %   uniformGrid = approx.mesh.UniformGrid([10, 10], [0, 0, 1, 1]);
    %   nonuniformGrid = approx.mesh.NonuniformGrid({0:0.1:1, [0, 0.2, 0.5, 1]});
    %
    % Notes:
    %   This is an abstract class that cannot be instantiated directly.
    %   Concrete subclasses must implement the collocate and refine methods,
    %   as well as the protected abstract methods for grid-specific operations.
    %
    % See also:
    %   approx.mesh.UniformGrid, approx.mesh.NonuniformGrid,
    %   approx.mesh.SeparableGraph

    properties (Access = public)
        nDims % Number of spatial dimensions (positive integer)
        bbox % Bounding box coordinates [a1, b1, a2, b2, ...]
        resolution % Vector of grid dimensions per axis
        centroids % Cell array of element centroids per dimension
        spacings % Cell array of element spacings per dimension
        elements % Cell array of element indices per dimension
        boundary % Cell array of boundary element indices
        indexer % Multi-indexer for linear/multi-index conversion
    end

    properties (Dependent)
        nTotalElements % Total number of elements in the grid
        allElementLinearIndices % Linear indices of all elements
        allElementMultiIndices % Multi-indices of all elements
        boundaryElementLinearIndices % Linear indices of boundary elements
        boundaryElementMultiIndices % Multi-indices of boundary elements
        measure % Minimum element spacing (grid measure)
        magnitudes % Vector of element magnitudes
        graph % Underlying separable graph representation
    end

    methods
        function obj = Grid(nDims)
            % GRID Constructor for Grid.
            %
            %   obj = Grid(nDims) creates a grid object for the specified
            %   number of dimensions. This constructor is called by
            %   concrete subclasses during their initialization.
            %
            % Inputs:
            %   nDims - Number of spatial dimensions (positive integer)
            %
            % Outputs:
            %   obj - Constructed Grid object

            obj.nDims = nDims;
        end

        function X = collocateOnBoundary(obj, i, xRef)
            d = ceil(i / 2);
            m = obj.allElementMultiIndices;
            if mod(d, 2) == 1
                e = m(:, d) == 1;
            else
                e = m(:, d) == obj.resolution(d);
            end
            m = m(e, :);
            X = obj.context.mesh.collocate(xRef, m);
        end
    end

    methods (Abstract)
        % COLLOCATE Map points to grid coordinates.
        Y = collocate(obj, nodes, indices)

        % REFINE Refine the grid resolution.
        obj = refine(obj, n)
    end

    methods (Abstract, Access = protected)
        % GETMINSPACING Get the minimum element spacing.
        h = getMinSpacing(obj)

        % GETALLELEMENTMAGNITUDES Get element magnitudes.
        h = getAllElementMagnitudes(obj)

        % GETSEPARABLEGRAPH Get the underlying graph representation.
        G = getSeparableGraph(obj)
    end

    methods
        function n = get.nTotalElements(obj)
            % GET.NTOTALELEMENTS Get the total number of elements.

            n = prod(obj.resolution);
        end

        function L = get.allElementLinearIndices(obj)
            % GET.ALLELEMENTLINEARINDICES Get linear indices of all elements.

            L = (1:obj.nTotalElements).';
        end

        function M = get.allElementMultiIndices(obj)
            % GET.ALLELEMENTMULTIINDICES Get multi-indices of all elements.

            M = obj.indexer.factorToMulti(obj.elements);
        end

        function L = get.boundaryElementLinearIndices(obj)
            % GET.BOUNDARYELEMENTLINEARINDICES Get boundary element indices.

            L = cellfun(@(M) obj.indexer.multiToLinear(M), ...
                obj.boundaryElementMultiIndices, 'Un', 0);
        end

        function M = get.boundaryElementMultiIndices(obj)
            % GET.BOUNDARYELEMENTMULTIINDICES Get boundary multi-indices.

            M = cellfun(@(F) obj.indexer.factorToMulti(F), obj.boundary, 'Un', 0);
        end

        function h = get.measure(obj)
            % GET.MEASURE Get the grid measure (minimum spacing).

            h = obj.getMinSpacing();
        end

        function h = get.magnitudes(obj)
            % GET.MAGNITUDES Get element magnitudes.

            h = obj.getAllElementMagnitudes();
        end

        function G = get.graph(obj)
            % GET.GRAPH Get the underlying graph representation.

            G = obj.getSeparableGraph();
        end
    end

    methods (Access = protected)
        function obj = setIndexer(obj)
            % SETINDEXER Initialize multi-indexer from resolution.
            %
            %   obj = setIndexer(obj) creates a MultiIndexer object for
            %   converting between linear and multi-dimensional indices
            %   based on the current grid resolution.
            %
            % Inputs:
            %   obj - The Grid object
            %
            % Outputs:
            %   obj - The Grid object

            obj.indexer = core.linalg.MultiIndexer(obj.resolution);
        end

        function obj = setElements(obj)
            % SETELEMENTS Initialize element factor indices.
            %
            %   obj = setElements(obj) creates factor index arrays for all
            %   elements based on the current grid resolution.
            %
            % Inputs:
            %   obj - The Grid object
            %
            % Outputs:
            %   obj - The Grid object

            obj.elements = arrayfun(@(n) 1:n, obj.resolution, 'Un', 0);
        end

        function obj = setBoundary(obj)
            % SETBOUNDARY Initialize boundary element indices.
            %
            %   obj = setBoundary(obj) creates factor index arrays for
            %   boundary elements on each face of the grid domain.
            %
            % Inputs:
            %   obj - The Grid object
            %
            % Outputs:
            %   obj - The Grid object

            obj.boundary = cell(1, 2*obj.nDims);
            for i = 1:obj.nDims
                b = arrayfun(@(r) 1:r, obj.resolution, 'Un', 0);
                b{i} = 1;
                obj.boundary{2*i-1} = b;
                b{i} = obj.resolution(i);
                obj.boundary{2*i} = b;
            end
        end
    end
end