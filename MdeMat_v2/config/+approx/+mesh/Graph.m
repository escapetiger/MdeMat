classdef Graph < handle
    % GRAPH Fundamental graph structure with vertices and edges.
    %
    %   Graph provides a basic representation of a graph data structure
    %   consisting of vertices (nodes) and edges (connections). This class
    %   serves as a foundation for more complex graph-based algorithms and
    %   mesh representations in computational applications.
    %
    % Examples:
    %   % Create empty graph
    %   graph = Graph();
    %   
    %   % Create graph with vertices and edges
    %   vertices = [0, 0; 1, 0; 1, 1; 0, 1];  % 2D square
    %   edges = [1, 2; 2, 3; 3, 4; 4, 1];     % Square connectivity
    %   graph = Graph(vertices, edges);
    %   
    %   % Query graph properties
    %   nVerts = graph.nVertices;
    %   nEdges = graph.nEdges;
    %   dims = graph.nDims;
    %
    % Notes:
    %   Vertices are stored as an nVertices × nDims matrix where each row
    %   represents a vertex. Edges are stored as an nEdges × 2 matrix where
    %   each row contains indices of connected vertices.
    %
    % See also:
    %   approx.mesh.SeparableGraph, approx.mesh.Grid
    
    properties (Access = public)
        vertices  % Matrix of vertex coordinates (nVertices × nDims)
        edges     % Matrix of edge connections (nEdges × 2)
    end

    properties (Dependent)
        nDims     % Number of spatial dimensions
        nVertices % Number of vertices in the graph
        nEdges    % Number of edges in the graph
    end
    
    methods
        function obj = Graph(vertices, edges)
            % GRAPH Constructor for Graph.
            %
            %   obj = Graph() creates an empty graph with no vertices
            %   or edges.
            %
            %   obj = Graph(vertices, edges) creates a graph with specified
            %   vertices and edges.
            %
            % Inputs:
            %   vertices - Matrix where each row represents a vertex (optional)
            %   edges - Matrix where each row represents an edge (optional)
            %
            % Outputs:
            %   obj - Constructed Graph object
            
            if nargin < 1, vertices = []; end
            if nargin < 2, edges = []; end
            
            obj.setVertices(vertices);
            obj.setEdges(edges);
        end

        function obj = setVertices(obj, vertices)
            % SETVERTICES Set the vertex coordinates.
            %
            %   obj = setVertices(obj, vertices) updates the graph vertices
            %   with the specified coordinate matrix.
            %
            % Inputs:
            %   obj - The Graph object
            %   vertices - Matrix where each row represents a vertex coordinate
            %
            % Outputs:
            %   obj - The Graph object

            core.except.assert(~isempty(vertices) && ismatrix(vertices), ...
                'InvalidVertices', 'Vertices must be a numeric matrix.');

            obj.vertices = vertices;
        end

        function obj = setEdges(obj, edges)
            % SETEDGES Set the edge connectivity.
            %
            %   obj = setEdges(obj, edges) updates the graph edges with the
            %   specified connectivity matrix. Validates that edge indices
            %   reference existing vertices.
            %
            % Inputs:
            %   obj - The Graph object
            %   edges- Matrix where each row contains two vertex indices
            %
            % Outputs:
            %   obj - The Graph object

            if ~isempty(edges)
                core.except.assert(ismatrix(edges) && size(edges, 2) == 2, ...
                    'InvalidEdges', ...
                    'Edges must be a numeric matrix with 2 columns');
                
                n = size(obj.vertices, 1);
                core.except.assert(all(edges(:) <= n & edges(:) >= 1), ...
                    'InvalidEdgeIndices', ...
                    'Edge indices must reference valid vertices');
            end

            obj.edges = edges;
        end

        function n = get.nDims(obj)
            % GET.NDIMS Get the number of spatial dimensions.

            n = size(obj.vertices, 2);
        end

        function n = get.nVertices(obj)
            % GET.NVERTICES Get the number of vertices.

            n = size(obj.vertices, 1);
        end

        function n = get.nEdges(obj)
            % GET.NEDGES Get the number of edges.

            n = size(obj.edges, 1);
        end
    end
end