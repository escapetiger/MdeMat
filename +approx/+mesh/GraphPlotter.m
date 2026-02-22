classdef GraphPlotter < handle
    % GRAPHPLOTTER Visualization class for Graph objects.
    %
    %   GraphPlotter provides comprehensive visualization functionality for
    %   approx.mesh.Graph objects, enabling customizable plots of vertices
    %   and edges in 1D, 2D, and 3D spaces. The class supports various
    %   styling options for different visualization needs.
    %
    % See also:
    %   approx.mesh.Graph, approx.mesh.MeshPlotter

    methods
        function plot(obj, graph, options)
            % PLOT Visualize a Graph object.
            %
            %   plot(obj, graph) creates a visualization of the graph with
            %   default styling options.
            %
            %   plot(obj, graph, options) allows customization of the
            %   visualization through options structure.

            arguments
                obj approx.mesh.GraphPlotter
                graph approx.mesh.Graph
                options.VertexColor = 'b'
                options.EdgeColor = 'k'
                options.VertexSize{mustBePositive} = 36
                options.EdgeWidth{mustBePositive} = 1
                options.ShowVertices{mustBeNumericOrLogical} = true
                options.ShowEdges{mustBeNumericOrLogical} = true
            end

            %< Prepare the plot
            holdState = ishold;
            if ~holdState
                cla;
            end
            hold on;

            %< Get the dimension of the graph
            nDims = graph.NDims;

            %< Plot based on dimension
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

            %< Set plot properties
            grid on;
            axis equal;

            %< Add title
            title(sprintf('%dD Graph: %d vertices, %d edges', ...
                nDims, graph.NVertices, graph.NEdges));

            %< Restore hold state
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

%< Plot edges if enabled
if options.ShowEdges && ~isequal(options.EdgeColor, 'none') && ...
        graph.NEdges > 0
    for iEdge = 1:graph.NEdges
        v1 = graph.Edges(iEdge, 1);
        v2 = graph.Edges(iEdge, 2);
        x = [graph.Vertices(v1), graph.Vertices(v2)];
        y = [0, 0];
        plot(x, y, '-', 'Color', options.EdgeColor, ...
            'LineWidth', options.EdgeWidth);
    end
end

%< Plot vertices if enabled
if options.ShowVertices && ~isequal(options.VertexColor, 'none') && ...
        graph.NVertices > 0
    scatter(graph.Vertices, zeros(graph.NVertices, 1), ...
        options.VertexSize, options.VertexColor, 'filled');
end

%< Set y-axis limits for better visualization
ylim([-0.5, 0.5]);
xlabel('X');
ylabel('');

%< Add reference line at y=0 if no edges are shown
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

%< Plot edges if enabled
if options.ShowEdges && ~isequal(options.EdgeColor, 'none') && ...
        graph.NEdges > 0
    for iEdge = 1:graph.NEdges
        v1 = graph.Edges(iEdge, 1);
        v2 = graph.Edges(iEdge, 2);
        x = [graph.Vertices(v1, 1), graph.Vertices(v2, 1)];
        y = [graph.Vertices(v1, 2), graph.Vertices(v2, 2)];
        plot(x, y, '-', 'Color', options.EdgeColor, ...
            'LineWidth', options.EdgeWidth);
    end
end

%< Plot vertices if enabled
if options.ShowVertices && ~isequal(options.VertexColor, 'none') && ...
        graph.NVertices > 0
    scatter(graph.Vertices(:, 1), graph.Vertices(:, 2), ...
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

%< Plot edges if enabled
if options.ShowEdges && ~isequal(options.EdgeColor, 'none') && ...
        graph.NEdges > 0
    for iEdge = 1:graph.NEdges
        v1 = graph.Edges(iEdge, 1);
        v2 = graph.Edges(iEdge, 2);
        x = [graph.Vertices(v1, 1), graph.Vertices(v2, 1)];
        y = [graph.Vertices(v1, 2), graph.Vertices(v2, 2)];
        z = [graph.Vertices(v1, 3), graph.Vertices(v2, 3)];
        plot3(x, y, z, '-', 'Color', options.EdgeColor, ...
            'LineWidth', options.EdgeWidth);
    end
end

%< Plot vertices if enabled
if options.ShowVertices && ~isequal(options.VertexColor, 'none') && ...
        graph.NVertices > 0
    scatter3(graph.Vertices(:, 1), graph.Vertices(:, 2), ...
        graph.Vertices(:, 3), options.VertexSize, ...
        options.VertexColor, 'filled');
end

xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);
end
