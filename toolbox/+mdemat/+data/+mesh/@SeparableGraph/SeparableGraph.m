classdef SeparableGraph < mdemat.data.variable.Separable
    % SEPARABLEGRAPH Memory-efficient structured graph representation.
    %
    %   SeparableGraph provides a memory-efficient representation of
    %   structured graphs using tensor decomposition concepts. Instead of
    %   storing all vertices and edges explicitly, it stores projections
    %   along each dimension as separate Graph objects, enabling efficient
    %   operations on large structured grids.
    %
    %   This representation is most efficient for structured grids where
    %   connectivity follows a tensor product pattern. Memory usage scales
    %   as O(sum of factor sizes) rather than O(product of factor sizes).
    %
    % See also:
    %   mdemat.data.variable.Separable, mdemat.data.mesh.Graph

    properties (Dependent)
        NDims % Number of dimensions
        NVertices % Vector of vertex counts per dimension
        NTotalVertices % Total number of vertices in full graph
        NEdges % Vector of edge counts per dimension
        NTotalEdges % Total number of edges in full graph
        Vertices % Cell array of factor vertex coordinates
        Edges % Cell array of factor edge connections
    end

    methods
        function obj = SeparableGraph(V, E)
            % SEPARABLEGRAPH Constructor for SeparableGraph.
            %
            %   obj = SeparableGraph(V, E) creates a separable graph from
            %   factor vertex coordinates @a V and edge connectivity @a E arrays.

            arguments
                V {mustBeNonempty, mustBeA(V, 'cell')}
                E {mustBeNonempty, mustBeA(E, 'cell')}
            end

            assert(all(length(V) == length(E), 'all'), ...
                'mdemat:data:mesh:SeparableGraph:DimensionMismatch', ...
                '%s', 'Vertices and edges must have the same dimension');

            factors = cell(1, length(V));
            for i = 1:length(V)
                v = V{i};
                e = E{i};
                assert( ...
                    all(isvector(v) && isnumeric(v), 'all'), ...
                    'mdemat:data:mesh:SeparableGraph:InvalidVertices', ...
                    '%s', 'Each entry in vertices must be a numeric vector');

                assert( ...
                    all(ismatrix(e) && isnumeric(e) && size(e, 2) == 2, 'all'), ...
                    'mdemat:data:mesh:SeparableGraph:InvalidEdges', ...
                    '%s', 'Each entry in edges must be a numeric matrix with 2 columns');

                n = length(v);
                assert(all(isreal(e) && all(isfinite(e(:))) && ...
                    all(e(:) == fix(e(:))) && all(e(:) <= n & e(:) >= 1), 'all'), ...
                    'mdemat:data:mesh:SeparableGraph:InvalidEdgeIndices', ...
                    '%s', 'Edge indices must reference valid vertices.');

                factors{i} = mdemat.data.mesh.Graph(v(:), e);
            end

            obj@mdemat.data.variable.Separable(factors=factors);
        end

        function n = get.NDims(obj)
            % GET.NDIMS Get the number of dimensions.

            n = length(obj.Factors);
        end

        function n = get.NVertices(obj)
            % GET.NVERTICES Get vertex counts per dimension.

            n = arrayfun(@(x) x.NVertices, obj.Factors);
        end

        function n = get.NTotalVertices(obj)
            % GET.NTOTALVERTICES Get total number of vertices.

            n = prod(obj.NVertices);
        end

        function n = get.NEdges(obj)
            % GET.NEDGES Get edge counts per dimension.

            n = arrayfun(@(x) x.NEdges, obj.Factors);
        end

        function n = get.NTotalEdges(obj)
            % GET.NTOTALEDGES Get total number of edges.

            m = obj.NVertices;
            n = sum(obj.NEdges.*prod(m)./m);
        end

        function V = get.Vertices(obj)
            % GET.VERTICES Get factor vertex coordinates.

            V = arrayfun(@(g) g.Vertices, obj.Factors, 'Un', 0);
        end

        function E = get.Edges(obj)
            % GET.EDGES Get factor edge connectivity.

            E = arrayfun(@(g) g.Edges, obj.Factors, 'Un', 0);
        end

        function G = full(obj)
            % FULL Expand to full Graph representation.
            %
            %   G = full(obj) expands the separable representation to
            %   a complete Graph object containing all vertices and edges
            %   in the tensor product structure.

            arguments
                obj mdemat.data.mesh.SeparableGraph
            end

            V = obj.getFullVertices();
            E = obj.getFullEdges();
            G = mdemat.data.mesh.Graph(V, E);
        end
    end

    methods (Access = protected)
        function V = getFullVertices(obj)
            % GETFULLVERTICES Generate full vertex array.

            V = obj.Vertices;
            [V{:}] = ndgrid(V{:});
            V = cellfun(@(x) x(:), V, 'Un', 0);
            V = cat(2, V{:});
        end

        function E = getFullEdges(obj)
            % GETFULLEDGES Generate full edge connectivity matrix.

            D = obj.NDims;

            if D == 1
                E = obj.Factors(1).Edges;
                return;
            end

            N = obj.NVertices;
            indexer = mdemat.data.index.MultiIndexer(shape=N);

            E = zeros(obj.NTotalEdges, 2);
            iEdge = 0;
            for d = 1:D
                factorEdges = obj.Factors(d).Edges;
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
                    batchEnd = min(batchStart+batchSize-1, P);
                    batchSize = batchEnd - batchStart + 1;

                    indexer2 = mdemat.data.index.MultiIndexer(shape=S);
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
