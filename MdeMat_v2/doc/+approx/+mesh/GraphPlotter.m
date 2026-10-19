classdef GraphPlotter < handle
    % GRAPHPLOTTER Visualization class for Graph objects.
    %
    %   GraphPlotter provides comprehensive visualization functionality for
    %   approx.mesh.Graph objects, enabling customizable plots of vertices
    %   and edges in 1D, 2D, and 3D spaces. The class supports various
    %   styling options for different visualization needs.
    %
    % Examples:
    %   % Create and plot a 2D graph
    %   vertices = [0, 0; 1, 0; 1, 1; 0, 1];
    %   edges = [1, 2; 2, 3; 3, 4; 4, 1];
    %   graph = approx.mesh.Graph(vertices, edges);
    %   plotter = approx.mesh.GraphPlotter();
    %   plotter.plot(graph);
    %
    %   % Plot with custom styling
    %   plotter.plot(graph, 'VertexColor', 'red', 'EdgeColor', 'blue', ...
    %                'VertexSize', 50, 'EdgeWidth', 2);
    %
    %   % Plot vertices only
    %   plotter.plot(graph, 'EdgeColor', 'none');
    %
    % See also:
    %   approx.mesh.Graph, approx.mesh.MeshPlotter

    methods
        function plot(obj, graph, varargin)
            % PLOT Visualize a Graph object.
            %
            %   plot(obj, graph) creates a visualization of the graph with
            %   default styling options.
            %
            %   plot(obj, graph, Name, Value, ...) allows customization of
            %   the visualization through name-value pair arguments.
            %
            % Inputs:
            %   obj - The GraphPlotter object
            %   graph - approx.mesh.Graph object to visualize
            %   varargin - Name-value pairs for customization
            %
            % Name-Value Pairs:
            %   'VertexColor' - Color for vertices (default: 'b')
            %   'EdgeColor' - Color for edges (default: 'k')
            %   'VertexSize' - Size of vertex markers (default: 36)
            %   'EdgeWidth' - Width of edge lines (default: 1)
            %   'ShowVertices' - Flag to show vertices (default: true)
            %   'ShowEdges' - Flag to show edges (default: true)

            % Validate input
            core.except.assert(isa(graph, 'approx.mesh.Graph'), ...
                'InvalidInput', 'Input must be an approx.mesh.Graph object.');

            % Parse input parameters
            p = inputParser;
            addParameter(p, 'VertexColor', 'b');
            addParameter(p, 'EdgeColor', 'k');
            addParameter(p, 'VertexSize', 36);
            addParameter(p, 'EdgeWidth', 1);
            addParameter(p, 'ShowVertices', true);
            addParameter(p, 'ShowEdges', true);
            parse(p, varargin{:});

            options = p.Results;

            % Prepare the plot
            holdState = ishold;
            if ~holdState
                cla;
            end
            hold on;

            % Get the dimension of the graph
            nDims = graph.nDims;

            % Plot based on dimension
            switch nDims
                case 1
                    plot1D(graph, options);
                case 2
                    plot2D(graph, options);
                case 3
                    plot3D(graph, options);
                otherwise
                    core.except.assert(false, 'UnsupportedDimension', ...
                        'Plotting is only supported for 1D, 2D and 3D graphs.');
            end

            % Set plot properties
            grid on;
            axis equal;

            % Add title
            title(sprintf('%dD Graph: %d vertices, %d edges', ...
                nDims, graph.nVertices, graph.nEdges));

            % Restore hold state
            if ~holdState
                hold off;
            end
        end
    end
end

function plot1D(graph, options)
% PLOT1D Visualize a 1D graph.
%
%   plot1D(graph, options) creates a 1D visualization by plotting vertices
%   along the x-axis at y=0 and connecting them with edges.
%
% Inputs:
%   graph - approx.mesh.Graph object with nDims == 1
%   options - Structure containing plotting options

% Plot edges if enabled
if options.ShowEdges && ~isequal(options.EdgeColor, 'none') && ...
        graph.nEdges > 0
    for iEdge = 1:graph.nEdges
        v1 = graph.edges(iEdge, 1);
        v2 = graph.edges(iEdge, 2);
        x = [graph.vertices(v1), graph.vertices(v2)];
        y = [0, 0];
        plot(x, y, '-', 'Color', options.EdgeColor, ...
            'LineWidth', options.EdgeWidth);
    end
end

% Plot vertices if enabled
if options.ShowVertices && ~isequal(options.VertexColor, 'none') && ...
        graph.nVertices > 0
    scatter(graph.vertices, zeros(graph.nVertices, 1), ...
        options.VertexSize, options.VertexColor, 'filled');
end

% Set y-axis limits for better visualization
ylim([-0.5, 0.5]);
xlabel('X');
ylabel('');

% Add reference line at y=0 if no edges are shown
if ~options.ShowEdges || isequal(options.EdgeColor, 'none')
    xLimits = xlim;
    plot(xLimits, [0, 0], 'k-', 'LineWidth', 0.5, 'Color', [0.7, 0.7, 0.7]);
end
end

function plot2D(graph, options)
% PLOT2D Visualize a 2D graph.
%
%   plot2D(graph, options) creates a 2D visualization by plotting vertices
%   as points and edges as line segments in the xy-plane.
%
% Inputs:
%   graph - approx.mesh.Graph object with nDims == 2
%   options - Structure containing plotting options

% Plot edges if enabled
if options.ShowEdges && ~isequal(options.EdgeColor, 'none') && ...
        graph.nEdges > 0
    for iEdge = 1:graph.nEdges
        v1 = graph.edges(iEdge, 1);
        v2 = graph.edges(iEdge, 2);
        x = [graph.vertices(v1, 1), graph.vertices(v2, 1)];
        y = [graph.vertices(v1, 2), graph.vertices(v2, 2)];
        plot(x, y, '-', 'Color', options.EdgeColor, ...
            'LineWidth', options.EdgeWidth);
    end
end

% Plot vertices if enabled
if options.ShowVertices && ~isequal(options.VertexColor, 'none') && ...
        graph.nVertices > 0
    scatter(graph.vertices(:, 1), graph.vertices(:, 2), ...
        options.VertexSize, options.VertexColor, 'filled');
end

xlabel('X');
ylabel('Y');
end

function plot3D(graph, options)
% PLOT3D Visualize a 3D graph.
%
%   plot3D(graph, options) creates a 3D visualization by plotting vertices
%   as points and edges as line segments in 3D space with appropriate
%   viewing angle.
%
% Inputs:
%   graph - approx.mesh.Graph object with nDims == 3
%   options - Structure containing plotting options

% Plot edges if enabled
if options.ShowEdges && ~isequal(options.EdgeColor, 'none') && ...
        graph.nEdges > 0
    for iEdge = 1:graph.nEdges
        v1 = graph.edges(iEdge, 1);
        v2 = graph.edges(iEdge, 2);
        x = [graph.vertices(v1, 1), graph.vertices(v2, 1)];
        y = [graph.vertices(v1, 2), graph.vertices(v2, 2)];
        z = [graph.vertices(v1, 3), graph.vertices(v2, 3)];
        plot3(x, y, z, '-', 'Color', options.EdgeColor, ...
            'LineWidth', options.EdgeWidth);
    end
end

% Plot vertices if enabled
if options.ShowVertices && ~isequal(options.VertexColor, 'none') && ...
        graph.nVertices > 0
    scatter3(graph.vertices(:, 1), graph.vertices(:, 2), ...
        graph.vertices(:, 3), options.VertexSize, ...
        options.VertexColor, 'filled');
end

xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);
end
