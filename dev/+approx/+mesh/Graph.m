classdef Graph < handle
    % GRAPH Fundamental graph structure with vertices and edges.
    %
    %   Graph provides a basic representation of a graph data structure
    %   consisting of vertices (nodes) and edges (connections). This class
    %   serves as a foundation for more complex graph-based algorithms and
    %   mesh representations in computational applications.
    %
    %   Vertices are stored as an nVertices × nDims matrix where each row
    %   represents a vertex. Edges are stored as an nEdges × 2 matrix where
    %   each row contains indices of connected vertices.
    %
    % See also:
    %   approx.mesh.SeparableGraph, approx.mesh.Grid
    
    properties (Access = public)
        Vertices  % Matrix of vertex coordinates (nVertices × nDims)
        Edges     % Matrix of edge connections (nEdges × 2)
    end

    properties (Dependent)
        NDims     % Number of spatial dimensions
        NVertices % Number of vertices in the graph
        NEdges    % Number of edges in the graph
    end
    
    methods
        function obj = Graph(vertices, edges)
            % GRAPH Constructor for Graph.
            %
            %   obj = Graph() creates an empty graph with no vertices
            %   or edges.
            %
            %   obj = Graph(vertices, edges) creates a graph with specified
            %   @a vertices and @a edges.

            arguments
                vertices {mustBeNumeric} = []
                edges {mustBeNumeric} = []
            end
            
            obj.setVertices(vertices);
            obj.setEdges(edges);
        end

        function obj = setVertices(obj, vertices)
            % SETVERTICES Set the vertex coordinates.
            %
            %   obj = setVertices(obj, vertices) updates the graph vertices
            %   with the specified coordinate matrix @a vertices.

            arguments
                obj approx.mesh.Graph
                vertices {mustBeNumeric}
            end

            obj.Vertices = vertices;
        end

        function obj = setEdges(obj, edges)
            % SETEDGES Set the edge connectivity.
            %
            %   obj = setEdges(obj, edges) updates the graph edges with the
            %   specified connectivity matrix @a edges. Validates that edge 
            %   indices reference existing vertices.

            arguments
                obj approx.mesh.Graph
                edges {mustBeNumeric}
            end

            if ~isempty(edges)
                core.except.assert(ismatrix(edges) && size(edges, 2) == 2, ...
                    'InvalidEdges', ...
                    'Edges must be a numeric matrix with 2 columns');
                
                n = size(obj.Vertices, 1);
                core.except.assert(all(edges(:) <= n & edges(:) >= 1), ...
                    'InvalidEdgeIndices', ...
                    'Edge indices must reference valid vertices');
            end

            obj.Edges = edges;
        end

        function n = get.NDims(obj)
            % GET.NDIMS Get the number of spatial dimensions.

            n = size(obj.Vertices, 2);
        end

        function n = get.NVertices(obj)
            % GET.NVERTICES Get the number of vertices.

            n = size(obj.Vertices, 1);
        end

        function n = get.NEdges(obj)
            % GET.NEDGES Get the number of edges.

            n = size(obj.Edges, 1);
        end
    end
end