classdef SeparableGraphPlotter < approx.mesh.GraphPlotter
    % SEPARABLEGRAPHPLOTTER Visualization class for SeparableGraph objects.
    %
    %   SeparableGraphPlotter provides memory-efficient visualization
    %   functionality for approx.mesh.SeparableGraph objects without
    %   requiring conversion to full Graph representations. The class
    %   supports direct tensor-product visualization, subsampling for
    %   large grids, and projection for high-dimensional cases.
    %
    % See also:
    %   approx.mesh.GraphPlotter, approx.mesh.SeparableGraph

    methods
        function plot(obj, graph, options)
            % PLOT Visualize a SeparableGraph object.
            %
            %   plot(obj, graph) creates a visualization of the separable
            %   graph using memory-efficient tensor-product structure.
            %
            %   plot(obj, graph, options) allows customization through
            %   options structure.

            arguments
                obj
                graph
                options.VertexColor = 'b'
                options.EdgeColor = 'k'
                options.VertexSize{mustBePositive} = 36
                options.EdgeWidth{mustBePositive} = 1
                options.ShowVertices{mustBeNumericOrLogical} = true
                options.ShowEdges{mustBeNumericOrLogical} = true
                options.MaxVertices{mustBePositive} = 10000
                options.ProjectionDims{mustBePositive, mustBeInteger} = [1, 2, 3]
            end

            %< Check if this is a SeparableGraph
            if isa(graph, 'approx.mesh.SeparableGraph')
                %< Choose visualization strategy based on graph properties
                if graph.NTotalVertices > options.MaxVertices
                    %< Large grids: use subsampling
                    obj.plotWithSubsampling(graph, options);
                elseif graph.NDims <= 3
                    %< Small 1D/2D/3D grids: direct tensor visualization
                    obj.plotTensorDirect(graph, options);
                else
                    %< High-dimensional grids: use projection
                    obj.plotWithProjection(graph, options);
                end
            else
                %< Fall back to parent GraphPlotter for regular Graph objects
                plot@approx.mesh.GraphPlotter(obj, graph, options);
            end
        end
    end

    methods (Access = private)
        function plotTensorDirect(~, graph, options)
            % PLOTTENSORDIRECT Direct visualization using tensor structure.
            %
            %   plotTensorDirect(obj, graph, options) visualizes the
            %   separable graph directly using its tensor-product structure
            %   without expanding to full representation.

            %< Prepare the plot
            holdState = ishold;
            if ~holdState
                cla;
            end
            hold on;

            %< Plot based on dimensionality
            switch graph.NDims
                case 1
                    plotTensor1D(graph, options);
                case 2
                    plotTensor2D(graph, options);
                case 3
                    plotTensor3D(graph, options);
            end

            %< Set common plot properties
            axis equal;
            grid on;

            %< Add informative title
            title(sprintf('%dD SeparableGraph: %d total vertices, %d total edges', ...
                graph.NDims, graph.NTotalVertices, graph.NTotalEdges));

            %< Restore hold state
            if ~holdState
                hold off;
            end
        end

        function plotWithSubsampling(obj, graph, options)
            % PLOTWITHSUBSAMPLING Visualize large graphs using subsampling.
            %
            %   plotWithSubsampling(obj, graph, options) creates a
            %   subsampled version of the graph for visualization when
            %   the total number of vertices exceeds the threshold.

            %< Calculate subsampling factor
            subsampleFactor = ceil(sqrt(graph.NTotalVertices/options.MaxVertices));

            %< Create subsampled factors
            subsampledVertices = cell(1, graph.NDims);
            subsampledEdges = cell(1, graph.NDims);

            for d = 1:graph.NDims
                vertices = graph.Vertices{d};
                edges = graph.Edges{d};

                %< Subsample vertices uniformly
                n = length(vertices);
                indices = 1:subsampleFactor:n;
                subsampledVertices{d} = vertices(indices);

                %< Adjust edges for subsampled vertices
                if ~isempty(edges)
                    %< Map original edge indices to subsampled indices
                    [~, newV1] = ismember(edges(:, 1), indices);
                    [~, newV2] = ismember(edges(:, 2), indices);

                    %< Keep only edges where both vertices are retained
                    validEdges = (newV1 > 0) & (newV2 > 0) & (newV1 ~= newV2);
                    subsampledEdges{d} = [newV1(validEdges), newV2(validEdges)];
                else
                    subsampledEdges{d} = [];
                end
            end

            %< Create subsampled graph and plot
            subsampledGraph = approx.mesh.SeparableGraph(subsampledVertices, subsampledEdges);
            obj.plotTensorDirect(subsampledGraph, options);

            %< Update title to indicate subsampling
            currentTitle = get(get(gca, 'Title'), 'String');
            newTitle = sprintf('%s (subsampled by factor %d)', currentTitle, subsampleFactor);
            title(newTitle);
        end

        function plotWithProjection(obj, graph, options)
            % PLOTWITHPROJECTION Visualize high-dimensional graphs using
            % projection.
            %
            %   plotWithProjection(obj, graph, options) creates a
            %   lower-dimensional projection of high-dimensional separable
            %   graphs for visualization.


            %< Validate and adjust projection dimensions
            projDims = options.ProjectionDims;
            projDims = projDims(projDims <= graph.NDims);
            projDims = projDims(1:min(3, length(projDims)));

            if isempty(projDims)
                projDims = 1:min(3, graph.NDims);
            end

            %< Create projected graph factors
            projVertices = graph.Vertices(projDims);
            projEdges = graph.Edges(projDims);

            %< Create and plot projected graph
            projGraph = approx.mesh.SeparableGraph(projVertices, projEdges);
            obj.plotTensorDirect(projGraph, options);

            %< Update title to indicate projection
            currentTitle = get(get(gca, 'Title'), 'String');
            newTitle = sprintf('%s (projected onto dims %s)', ...
                currentTitle, mat2str(projDims));
            title(newTitle);
        end
    end
end


function plotTensor1D(graph, options)
% PLOTTENSOR1D Plot 1D separable graph.
%
%   plotTensor1D(obj, graph, options) creates a 1D
%   visualization by plotting vertices along the x-axis.

vertices = graph.Vertices{1};
edges = graph.Edges{1};

%< Plot vertices
if options.ShowVertices && ~isequal(options.VertexColor, 'none')
    scatter(vertices, zeros(size(vertices)), ...
        options.VertexSize, options.VertexColor, 'filled');
end

%< Plot edges
if options.ShowEdges && ~isequal(options.EdgeColor, 'none') && ~isempty(edges)
    for iEdge = 1:size(edges, 1)
        v1 = edges(iEdge, 1);
        v2 = edges(iEdge, 2);
        plot([vertices(v1), vertices(v2)], [0, 0], '-', ...
            'Color', options.EdgeColor, 'LineWidth', options.EdgeWidth);
    end
end

%< Set axis properties
ylim([-0.5, 0.5]);
xlabel('X');
end

function plotTensor2D(graph, options)
% PLOTTENSOR2D Plot 2D separable graph using tensor structure.
%
%   plotTensor2D(obj, graph, options) creates a 2D grid
%   visualization using the tensor-product structure without
%   expanding to full vertex coordinates.

vertices1 = graph.Vertices{1};
vertices2 = graph.Vertices{2};
edges1 = graph.Edges{1};
edges2 = graph.Edges{2};

%< Plot vertices using tensor product
if options.ShowVertices && ~isequal(options.VertexColor, 'none')
    [X, Y] = meshgrid(vertices1, vertices2);
    scatter(X(:), Y(:), options.VertexSize, ...
        options.VertexColor, 'filled');
end

%< Plot edges in first dimension (horizontal lines)
if options.ShowEdges && ~isequal(options.EdgeColor, 'none') && ~isempty(edges1)
    for j = 1:length(vertices2)
        for iEdge = 1:size(edges1, 1)
            v1 = edges1(iEdge, 1);
            v2 = edges1(iEdge, 2);
            plot([vertices1(v1), vertices1(v2)], ...
                [vertices2(j), vertices2(j)], '-', ...
                'Color', options.EdgeColor, 'LineWidth', options.EdgeWidth);
        end
    end
end

%< Plot edges in second dimension (vertical lines)
if options.ShowEdges && ~isequal(options.EdgeColor, 'none') && ~isempty(edges2)
    for i = 1:length(vertices1)
        for iEdge = 1:size(edges2, 1)
            v1 = edges2(iEdge, 1);
            v2 = edges2(iEdge, 2);
            plot([vertices1(i), vertices1(i)], ...
                [vertices2(v1), vertices2(v2)], '-', ...
                'Color', options.EdgeColor, 'LineWidth', options.EdgeWidth);
        end
    end
end

%< Set axis properties
xlabel('X');
ylabel('Y');
end

function plotTensor3D(graph, options)
% PLOTTENSOR3D Plot 3D separable graph using tensor structure.
%
%   plotTensor3D(graph, options) creates a 3D grid
%   visualization using the tensor-product structure.

vertices1 = graph.Vertices{1};
vertices2 = graph.Vertices{2};
vertices3 = graph.Vertices{3};
edges1 = graph.Edges{1};
edges2 = graph.Edges{2};
edges3 = graph.Edges{3};

%< Plot vertices using tensor product
if options.ShowVertices && ~isequal(options.VertexColor, 'none')
    [X, Y, Z] = meshgrid(vertices1, vertices2, vertices3);
    scatter3(X(:), Y(:), Z(:), options.VertexSize, ...
        options.VertexColor, 'filled');
end

%< Plot edges in first dimension (X direction)
if options.ShowEdges && ~isequal(options.EdgeColor, 'none') && ~isempty(edges1)
    for j = 1:length(vertices2)
        for k = 1:length(vertices3)
            for iEdge = 1:size(edges1, 1)
                v1 = edges1(iEdge, 1);
                v2 = edges1(iEdge, 2);
                plot3([vertices1(v1), vertices1(v2)], ...
                    [vertices2(j), vertices2(j)], ...
                    [vertices3(k), vertices3(k)], '-', ...
                    'Color', options.EdgeColor, 'LineWidth', options.EdgeWidth);
            end
        end
    end
end

%< Plot edges in second dimension (Y direction)
if options.ShowEdges && ~isequal(options.EdgeColor, 'none') && ~isempty(edges2)
    for i = 1:length(vertices1)
        for k = 1:length(vertices3)
            for iEdge = 1:size(edges2, 1)
                v1 = edges2(iEdge, 1);
                v2 = edges2(iEdge, 2);
                plot3([vertices1(i), vertices1(i)], ...
                    [vertices2(v1), vertices2(v2)], ...
                    [vertices3(k), vertices3(k)], '-', ...
                    'Color', options.EdgeColor, 'LineWidth', options.EdgeWidth);
            end
        end
    end
end

%< Plot edges in third dimension (Z direction)
if options.ShowEdges && ~isequal(options.EdgeColor, 'none') && ~isempty(edges3)
    for i = 1:length(vertices1)
        for j = 1:length(vertices2)
            for iEdge = 1:size(edges3, 1)
                v1 = edges3(iEdge, 1);
                v2 = edges3(iEdge, 2);
                plot3([vertices1(i), vertices1(i)], ...
                    [vertices2(j), vertices2(j)], ...
                    [vertices3(v1), vertices3(v2)], '-', ...
                    'Color', options.EdgeColor, 'LineWidth', options.EdgeWidth);
            end
        end
    end
end

%< Set axis properties
xlabel('X');
ylabel('Y');
zlabel('Z');
view(3);
end