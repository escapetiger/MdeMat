classdef SeparableGraph < core.linalg.Separable
    % SEPARABLEGRAPH Memory-efficient structured graph representation.
    %
    %   SeparableGraph provides a memory-efficient representation of
    %   structured graphs using tensor decomposition concepts. Instead of
    %   storing all vertices and edges explicitly, it stores projections
    %   along each dimension as separate Graph objects, enabling efficient
    %   operations on large structured grids.
    %
    % Examples:
    %   % Create 2D grid graph
    %   V = {[0, 1, 2], [0, 1]};  % 3×2 grid vertices
    %   E = {[1,2; 2,3], [1,2]};  % Edge connectivity per dimension
    %   sepGraph = SeparableGraph(V, E);
    %   
    %   % Query properties
    %   nVerts = sepGraph.nTotalVertices;
    %   nEdges = sepGraph.nTotalEdges;
    %   
    %   % Convert to full graph representation
    %   fullGraph = sepGraph.full();
    %
    % Notes:
    %   This representation is most efficient for structured grids where
    %   connectivity follows a tensor product pattern. Memory usage scales
    %   as O(sum of factor sizes) rather than O(product of factor sizes).
    %
    % See also:
    %   core.linalg.SeparableObject, approx.mesh.Graph,
    %   approx.mesh.Grid

    properties (Dependent)
        nDims           % Number of dimensions
        nVertices       % Vector of vertex counts per dimension
        nTotalVertices  % Total number of vertices in full graph
        nEdges          % Vector of edge counts per dimension
        nTotalEdges     % Total number of edges in full graph
        vertices        % Cell array of factor vertex coordinates
        edges           % Cell array of factor edge connections
    end
    
    methods
        function obj = SeparableGraph(V, E)
            % SEPARABLEGRAPH Constructor for SeparableGraph.
            %
            %   obj = SeparableGraph(V, E) creates a separable graph from
            %   factor vertex coordinates and edge connectivity arrays.
            %
            % Inputs:
            %   V - Cell array of vertex coordinate vectors per dimension
            %   E - Cell array of edge connectivity matrices per dimension
            %
            % Outputs:
            %   obj - Constructed SeparableGraph object
            
            core.except.assert(iscell(V), ...
                'InvalidVertices', 'Vertices must be a cell array.');
            
            core.except.assert(iscell(E), ...
                'InvalidEdges', 'Edges must be a cell array.');
            
            core.except.assert(length(V) == length(E), ...
                'DimensionMismatch', ...
                'Vertices and edges must have the same dimension');
            
            factors = cell(1, length(V));
            for i = 1:length(V)
                v = V{i};
                e = E{i};
                core.except.assert( ...
                    isvector(v) && isnumeric(v), ...
                    'InvalidVertices', ...
                    'Each entry in vertices must be a numeric vector');

                core.except.assert( ...
                    ismatrix(e) && isnumeric(e) && size(e, 2) == 2, ...
                    'InvalidEdges', ...
                    'Each entry in edges must be a numeric matrix with 2 columns');
                
                n = length(v);
                core.except.assert(all(e(:) <= n & e(:) >= 1), ...
                    'InvalidEdgeIndices', ...
                    'Edge indices must reference valid vertices.');
                
                factors{i} = approx.mesh.Graph(v(:), e);
            end
            
            obj@core.linalg.Separable(factors);
        end
        
        function n = get.nDims(obj)
            % GET.NDIMS Get the number of dimensions.

            n = length(obj.factors);
        end

        function n = get.nVertices(obj)
            % GET.NVERTICES Get vertex counts per dimension.

            n = arrayfun(@(x) x.nVertices, obj.factors);
        end

        function n = get.nTotalVertices(obj)
            % GET.NTOTALVERTICES Get total number of vertices.

            n = prod(obj.nVertices);
        end

        function n = get.nEdges(obj)
            % GET.NEDGES Get edge counts per dimension.

            n = arrayfun(@(x) x.nEdges, obj.factors);
        end

        function n = get.nTotalEdges(obj)
            % GET.NTOTALEDGES Get total number of edges.

            m = obj.nVertices;
            n = sum(obj.nEdges .* prod(m) ./ m);
        end
        
        function V = get.vertices(obj)
            % GET.VERTICES Get factor vertex coordinates.

            V = arrayfun(@(g) g.vertices, obj.factors, 'Un', 0);
        end
        
        function E = get.edges(obj)
            % GET.EDGES Get factor edge connectivity.

            E = arrayfun(@(g) g.edges, obj.factors, 'Un', 0);
        end

        function G = full(obj)
            % FULL Expand to full Graph representation.
            %
            %   G = full(obj) expands the separable representation to
            %   a complete Graph object containing all vertices and edges
            %   in the tensor product structure.
            %
            % Inputs:
            %   obj - The SeparableGraph object
            %
            % Outputs:
            %   G - Graph object with full vertex and edge arrays

            V = obj.getFullVertices();
            E = obj.getFullEdges();
            G = approx.mesh.Graph(V, E);
        end
    end
    
    methods (Access = protected)
        function V = getFullVertices(obj)
            % GETFULLVERTICES Generate full vertex array.
            %
            %   V = getFullVertices(obj) creates the complete vertex
            %   coordinate matrix by taking tensor products of factor
            %   vertex coordinates.
            %
            % Inputs:
            %   obj - The SeparableGraph object
            %
            % Outputs:
            %   V - Matrix of all vertex coordinates (nVertices × nDims)

            V = obj.vertices;
            [V{:}] = ndgrid(V{:});
            V = cellfun(@(x) x(:), V, 'Un', 0);
            V = cat(2, V{:});
        end

        function E = getFullEdges(obj)
            % GETFULLEDGES Generate full edge connectivity matrix.
            %
            %   E = getFullEdges(obj) creates the complete edge
            %   connectivity matrix by expanding factor edges to the full
            %   tensor product graph structure.
            %
            % Inputs:
            %   obj - The SeparableGraph object
            %
            % Outputs:
            %   E - Matrix of all edge connections (nEdges × 2)

            D = obj.nDims;

            if D == 1
                E = obj.factors(1).edges;
                return;
            end

            N = obj.nVertices;
            indexer = core.linalg.MultiIndexer(N);

            E = zeros(obj.nTotalEdges, 2);
            iEdge = 0;
            for d = 1:D
                factorEdges = obj.factors(d).edges;
                nFactorEdges = size(factorEdges, 1);
                
                if nFactorEdges == 0
                    continue;
                end
                
                I = 1:D;
                I(d) = [];
                S = N(I);
                P = prod(S);
                batchSize = min(10000, P);
                
                for batchStart = 1:batchSize:P
                    batchEnd = min(batchStart + batchSize - 1, P);
                    batchSize = batchEnd - batchStart + 1;
                    
                    indexer2 = core.linalg.MultiIndexer(S);
                    L = (batchStart:batchEnd)';
                    M = indexer2.linearToMulti(L);

                    for e = 1:nFactorEdges
                        v1 = factorEdges(e, 1);
                        v2 = factorEdges(e, 2);
                        
                        newEdges = zeros(batchSize, 2);
                        
                        for j = 1:batchSize
                            idx1 = ones(1, D);
                            idx2 = ones(1, D);
                            idx1(I) = M(j, :);
                            idx2(I) = M(j, :);
                            idx1(d) = v1;
                            idx2(d) = v2;
                            L1 = indexer.multiToLinear(idx1);
                            L2 = indexer.multiToLinear(idx2);
                            newEdges(j, :) = [L1, L2];
                        end
                        
                        E(iEdge+1:iEdge+batchSize, :) = newEdges;
                        iEdge = iEdge + batchSize;
                    end
                end
            end
        end
    end
end