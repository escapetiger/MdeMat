classdef MeshPlotter < handle
    % MESHPLOTTER Optimized visualization class for Mesh objects.
    %
    %   MeshPlotter provides fast visualization functionality for
    %   approx.mesh.Mesh objects using MATLAB's built-in plotting functions
    %   such as triplot, trimesh, trisurf, and vectorized plotting for
    %   maximum performance.
    %
    % Examples:
    %   % Create and plot a 2D simplex mesh
    %   vertices = [0, 0; 1, 0; 0.5, 1; 1, 1];
    %   elements = [1, 2, 3; 2, 4, 3];
    %   mesh = approx.mesh.SimplexMesh(vertices, elements);
    %   plotter = approx.mesh.MeshPlotter();
    %   plotter.plot(mesh);
    %
    %   % Plot with custom styling
    %   plotter.plot(mesh, 'ElementColor', 'cyan', 'BoundaryColor', 'red', ...
    %                'ShowFaces', true, 'Alpha', 0.5);
    %
    %   % Plot only boundary and vertices
    %   plotter.plot(mesh, 'ShowEdges', false, 'ShowElements', false);
    %
    % See also:
    %   approx.mesh.Mesh, approx.mesh.SimplexMesh, triplot, trimesh, trisurf

    properties (Access = protected)
        graphPlotter % GraphPlotter for visualizing the mesh's graph (GraphPlotter)
    end

    methods
        function obj = MeshPlotter()
            % MESHPLOTTER Constructor for MeshPlotter.
            %
            %   obj = MeshPlotter() creates a mesh plotter instance with
            %   an associated graph plotter for rendering vertices and edges.
            %
            % Outputs:
            %   obj - Constructed MeshPlotter object

            obj.graphPlotter = approx.mesh.GraphPlotter();
        end

        function plot(obj, mesh, varargin)
            % PLOT Visualize a Mesh object using optimized MATLAB functions.
            %
            %   plot(obj, mesh) creates a visualization of the mesh with
            %   default styling options using fast built-in functions.
            %
            %   plot(obj, mesh, Name, Value, ...) allows customization of
            %   the visualization through name-value pair arguments.
            %
            % Inputs:
            %   obj - The MeshPlotter object
            %   mesh - approx.mesh.Mesh object to visualize
            %   varargin - Name-value pairs for customization
            %
            % Name-Value Pairs:
            %   'VertexColor' - Color for vertices (default: 'b')
            %   'EdgeColor' - Color for edges (default: 'k')
            %   'FaceColor' - Color for faces (default: 'none')
            %   'ElementColor' - Color for elements (default: 'none')
            %   'BoundaryColor' - Color for boundary faces (default: 'r')
            %   'VertexSize' - Size of vertex markers (default: 36)
            %   'EdgeWidth' - Width of edge lines (default: 1)
            %   'Alpha' - Transparency for faces and elements (default: 0.3)
            %   'ShowVertices' - Flag to show vertices (default: true)
            %   'ShowEdges' - Flag to show edges (default: true)
            %   'ShowFaces' - Flag to show faces (default: false)
            %   'ShowElements' - Flag to show elements (default: false)
            %   'ShowBoundary' - Flag to show boundary (default: true)
            %   'UseFastRendering' - Use optimized built-ins (default: true)

            % Validate input
            core.except.assert(isa(mesh, 'approx.mesh.Mesh'), ...
                'InvalidInput', 'Input must be an approx.mesh.Mesh object.');

            % Parse input parameters
            p = inputParser;
            addParameter(p, 'VertexColor', 'b');
            addParameter(p, 'EdgeColor', 'k');
            addParameter(p, 'FaceColor', 'none');
            addParameter(p, 'ElementColor', 'none');
            addParameter(p, 'BoundaryColor', 'r');
            addParameter(p, 'VertexSize', 36);
            addParameter(p, 'EdgeWidth', 1);
            addParameter(p, 'Alpha', 0.3);
            addParameter(p, 'ShowVertices', true);
            addParameter(p, 'ShowEdges', true);
            addParameter(p, 'ShowFaces', false);
            addParameter(p, 'ShowElements', false);
            addParameter(p, 'ShowBoundary', true);
            addParameter(p, 'UseFastRendering', true);
            parse(p, varargin{:});

            options = p.Results;

            % Prepare the plot
            holdState = ishold;
            if ~holdState
                cla;
            end
            hold on;

            % Get the dimension of the mesh
            nDims = mesh.nDims;

            % Construct connectivity graph of mesh
            if p.Results.ShowEdges
                graph = mesh.graphify();
            else
                graph = [];
            end

            % Choose rendering method
            if options.UseFastRendering
                % Use optimized MATLAB built-ins
                switch nDims
                    case {1, 2, 3}
                        obj.plotMeshFast(mesh, graph, options);
                    otherwise
                        obj.plotHighDimensional(mesh, graph, options);
                end
            else
                % Use original method for compatibility
                switch nDims
                    case {1, 2, 3}
                        obj.plotMeshOriginal(mesh, graph, options);
                    otherwise
                        obj.plotHighDimensional(mesh, graph, options);
                end
            end

            % Set plot properties
            axis equal;
            grid on;

            % Add title with mesh information
            title(sprintf('%dD Mesh: %d vertices, %d elements, %d faces (%d boundary)', ...
                nDims, mesh.nVertices, mesh.nElements, mesh.nFaces, mesh.nBoundaryFaces));

            % Add legend
            addLegend(options);

            % Set view based on dimension
            if nDims == 3
                view(3);
            end

            % Restore hold state
            if ~holdState
                hold off;
            end
        end
    end

    methods (Access = protected)
        function plotMeshFast(~, mesh, graph, options)
            % PLOTMESHFAST Visualize mesh using optimized MATLAB functions.
            %
            %   plotMeshFast(obj, mesh, graph, options) renders mesh components
            %   using vectorized MATLAB built-in functions for maximum speed.
            %
            % Inputs:
            %   obj - The MeshPlotter object
            %   mesh - approx.mesh.Mesh object with nDims <= 3
            %   graph - approx.mesh.Graph object
            %   options - Structure containing plotting options

            % Check if mesh is triangular/tetrahedral for built-in functions
            isTriangular = isTriangularMesh(mesh);

            % Plot elements using fast methods
            if options.ShowElements && ~isequal(options.ElementColor, 'none') && ...
                    mesh.nElements > 0
                plotElementsFast(mesh, options, isTriangular);
            end

            % Plot faces using fast methods
            if options.ShowFaces && ~isequal(options.FaceColor, 'none') && ...
                    mesh.nFaces > 0
                plotFacesFast(mesh, options);
            end

            % Plot edges using fast methods
            if options.ShowEdges && ~isequal(options.EdgeColor, 'none') && ...
                    graph.nEdges > 0
                plotEdgesFast(mesh, graph.edges, options, isTriangular);
            end

            % Plot vertices using fast methods
            if options.ShowVertices && ~isequal(options.VertexColor, 'none') && ...
                    mesh.nVertices > 0
                plotVerticesFast(mesh, options);
            end

            % Plot boundary using fast methods
            if options.ShowBoundary && ~isequal(options.BoundaryColor, 'none') && ...
                    mesh.nBoundaryFaces > 0
                plotBoundaryFast(mesh, options);
            end
        end

        function plotMeshOriginal(obj, mesh, graph, options)
            % PLOTMESHORIGINAL Use original plotting method for compatibility.

            % Plot vertices and edges using GraphPlotter if enabled
            if options.ShowVertices || options.ShowEdges
                graphOptions = {};

                if options.ShowVertices
                    graphOptions = [graphOptions, {'VertexColor', options.VertexColor, ...
                        'VertexSize', options.VertexSize, 'ShowVertices', true}];
                else
                    graphOptions = [graphOptions, {'ShowVertices', false}];
                end

                if options.ShowEdges
                    graphOptions = [graphOptions, {'EdgeColor', options.EdgeColor, ...
                        'EdgeWidth', options.EdgeWidth, 'ShowEdges', true}];
                else
                    graphOptions = [graphOptions, {'ShowEdges', false}];
                end

                obj.graphPlotter.plot(graph, graphOptions{:});
            end

            % Plot faces if enabled
            if options.ShowFaces && ~isequal(options.FaceColor, 'none') && ...
                    mesh.nFaces > 0
                plotFacesSlow(mesh, options.FaceColor, options.Alpha);
            end

            % Plot elements if enabled
            if options.ShowElements && ~isequal(options.ElementColor, 'none') && ...
                    mesh.nElements > 0
                plotElementsSlow(mesh, options.ElementColor, options.Alpha);
            end

            % Plot boundary if enabled
            if options.ShowBoundary && ~isequal(options.BoundaryColor, 'none') && ...
                    mesh.nBoundaryFaces > 0
                plotBoundarySlow(mesh, options.BoundaryColor, options.EdgeWidth*1.5);
            end
        end

        function plotHighDimensional(obj, mesh, graph, options)
            % PLOTHIGHDIMENSIONAL Visualize a mesh with nDims > 3.

            % Get vertex coordinates
            vertices = mesh.vertices;

            % Perform PCA to project to 3D
            if size(vertices, 1) > 1
                [~, score] = pca(vertices);
                projectedVertices = score(:, 1:min(3, size(score, 2)));
            else
                projectedVertices = vertices(:, 1:min(3, size(vertices, 2)));
            end

            % Pad with zeros if necessary
            if size(projectedVertices, 2) < 3
                projectedVertices = [projectedVertices, ...
                    zeros(size(projectedVertices, 1), 3-size(projectedVertices, 2))];
            end

            % Create a 3D graph for visualization
            projectedGraph = approx.mesh.Graph(projectedVertices, graph.edges);

            % Use graph plotter to visualize
            graphOptions = {};

            if options.ShowVertices
                graphOptions = [graphOptions, {'VertexColor', options.VertexColor, ...
                    'VertexSize', options.VertexSize, 'ShowVertices', true}];
            else
                graphOptions = [graphOptions, {'ShowVertices', false}];
            end

            if options.ShowEdges
                graphOptions = [graphOptions, {'EdgeColor', options.EdgeColor, ...
                    'EdgeWidth', options.EdgeWidth, 'ShowEdges', true}];
            else
                graphOptions = [graphOptions, {'ShowEdges', false}];
            end

            obj.graphPlotter.plot(projectedGraph, graphOptions{:});

            % Add warning about projection
            core.except.verify(0, 'DimensionReduction', ...
                'Visualizing %dD mesh using PCA projection to 3D.', mesh.nDims);

            % Set 3D view
            view(3);
        end
    end
end

function isTriangular = isTriangularMesh(mesh)
% ISTRIANGULARMESH Check if mesh consists of triangular/tetrahedral elements.

nDims = mesh.nDims;
nVerticesPerElement = mesh.nVerticesPerElement;

% Check if elements are simplices (triangles or tetrahedra)
expectedVertices = nDims + 1;
isTriangular = (nVerticesPerElement == expectedVertices);
end

function surfaceTriangles = getTetrahedralSurface(mesh)
% GETTETRAHEDRALSURFACE Extract surface triangles from tetrahedral mesh.

if mesh.nDims ~= 3 || size(mesh.elements, 2) ~= 4
    surfaceTriangles = [];
    return;
end

% Get boundary faces which should be triangles
if mesh.nBoundaryFaces > 0 && size(mesh.faces, 2) == 3
    surfaceTriangles = mesh.faces(mesh.boundary, :);
else
    surfaceTriangles = [];
end
end

function plotElementsFast(mesh, options, isTriangularMesh)
% PLOTELEMENTSFAST Plot elements using optimized functions.

nDims = mesh.nDims;
vertices = mesh.vertices;
elements = mesh.elements;

switch nDims
    case 1
        % Plot line segments efficiently
        x = reshape(vertices(elements', 1), 2, []);
        y = zeros(size(x));
        plot(x, y, '-', 'Color', options.ElementColor, 'LineWidth', 3);

    case 2
        if isTriangularMesh
            % Use trisurf for filled triangles
            trisurf(elements, vertices(:, 1), vertices(:, 2), ...
                zeros(size(vertices, 1), 1), ...
                'FaceColor', options.ElementColor, ...
                'FaceAlpha', options.Alpha, ...
                'EdgeColor', 'none');
            view(2); % Set 2D view
        else
            % Use patch for general polygons
            patch('Faces', elements, 'Vertices', vertices, ...
                'FaceColor', options.ElementColor, ...
                'FaceAlpha', options.Alpha, ...
                'EdgeColor', 'none');
        end

    case 3
        if isTriangularMesh && size(elements, 2) == 4
            % Tetrahedral mesh - show surface using boundary
            boundaryTriangles = getTetrahedralSurface(mesh);
            if ~isempty(boundaryTriangles)
                trisurf(boundaryTriangles, vertices(:, 1), ...
                    vertices(:, 2), vertices(:, 3), ...
                    'FaceColor', options.ElementColor, ...
                    'FaceAlpha', options.Alpha, ...
                    'EdgeColor', 'none');
            end
        else
            % Use patch for general 3D elements
            plotElementsSlow(mesh, options.ElementColor, options.Alpha);
        end
end
end

function plotFacesFast(mesh, options)
% PLOTFACESFAST Plot faces using optimized functions.

nDims = mesh.nDims;
vertices = mesh.vertices;
faces = mesh.faces;

if nDims == 3 && size(faces, 2) == 3
    % Use trisurf for triangular faces in 3D
    trisurf(faces, vertices(:, 1), vertices(:, 2), vertices(:, 3), ...
        'FaceColor', options.FaceColor, ...
        'FaceAlpha', options.Alpha, ...
        'EdgeColor', 'none');
elseif nDims >= 2
    % Use patch for general faces
    patch('Faces', faces, 'Vertices', vertices, ...
        'FaceColor', options.FaceColor, ...
        'FaceAlpha', options.Alpha, ...
        'EdgeColor', 'none');
end
end

function plotEdgesFast(mesh, edges, options, isTriangularMesh)
% PLOTEDGESFAST Plot edges using optimized functions.

nDims = mesh.nDims;
vertices = mesh.vertices;
elements = mesh.elements;

switch nDims
    case 2
        if isTriangularMesh
            % Use triplot for triangular mesh edges
            triplot(elements, vertices(:, 1), vertices(:, 2), ...
                'Color', options.EdgeColor, ...
                'LineWidth', options.EdgeWidth);
        else
            % Plot edges directly using vectorized plotting
            plotEdgesVectorized(vertices, edges, options, nDims);
        end

    case 3
        if isTriangularMesh
            % Use trimesh for 3D triangular mesh edges
            if size(elements, 2) == 3
                trimesh(elements, vertices(:, 1), vertices(:, 2), ...
                    vertices(:, 3), 'EdgeColor', options.EdgeColor, ...
                    'LineWidth', options.EdgeWidth, 'FaceColor', 'none');
            else
                % For tetrahedra, plot edges directly
                plotEdgesVectorized(vertices, edges, options, nDims);
            end
        else
            % Plot edges directly
            plotEdgesVectorized(vertices, edges, options, nDims);
        end

    otherwise
        % Plot edges directly for other dimensions
        self.plotEdgesVectorized(vertices, edges, options, nDims);
end
end

function plotEdgesVectorized(vertices, edges, options, nDims)
% PLOTEDGESVECTORIZED Plot edges using vectorized plotting.

if isempty(edges)
    return;
end

% Create vectorized edge coordinates
x = [vertices(edges(:, 1), 1)'; vertices(edges(:, 2), 1)'; NaN(1, size(edges, 1))];
y = [vertices(edges(:, 1), 2)'; vertices(edges(:, 2), 2)'; NaN(1, size(edges, 1))];

if nDims >= 3
    z = [vertices(edges(:, 1), 3)'; vertices(edges(:, 2), 3)'; NaN(1, size(edges, 1))];
    plot3(x(:), y(:), z(:), '-', 'Color', options.EdgeColor, ...
        'LineWidth', options.EdgeWidth);
else
    plot(x(:), y(:), '-', 'Color', options.EdgeColor, ...
        'LineWidth', options.EdgeWidth);
end
end

function plotVerticesFast(mesh, options)
% PLOTVERTICESFAST Plot vertices using optimized functions.

vertices = mesh.vertices;
nDims = mesh.nDims;

switch nDims
    case 1
        scatter(vertices(:, 1), zeros(size(vertices, 1), 1), ...
            options.VertexSize, options.VertexColor, 'filled');
    case 2
        scatter(vertices(:, 1), vertices(:, 2), ...
            options.VertexSize, options.VertexColor, 'filled');
    case 3
        scatter3(vertices(:, 1), vertices(:, 2), vertices(:, 3), ...
            options.VertexSize, options.VertexColor, 'filled');
    otherwise
        % Use first 3 dimensions for visualization
        nPlotDims = min(nDims, 3);
        if nPlotDims == 3
            scatter3(vertices(:, 1), vertices(:, 2), vertices(:, 3), ...
                options.VertexSize, options.VertexColor, 'filled');
        else
            scatter(vertices(:, 1), vertices(:, 2), ...
                options.VertexSize, options.VertexColor, 'filled');
        end
end
end

function plotBoundaryFast(mesh, options)
% PLOTBOUNDARYFAST Plot boundary using optimized functions.

nDims = mesh.nDims;
vertices = mesh.vertices;
faces = mesh.faces;
boundary = mesh.boundary;

if isempty(boundary)
    return;
end

% Get boundary face vertices
boundaryFaces = faces(boundary, :);

switch nDims
    case 1
        % Boundary points in 1D
        boundaryVertices = vertices(unique(boundaryFaces(:)), :);
        scatter(boundaryVertices(:, 1), zeros(size(boundaryVertices, 1), 1), ...
            100, options.BoundaryColor, 'filled', 'MarkerEdgeColor', 'k');

    case 2
        % Boundary edges in 2D - plot as line segments
        if size(boundaryFaces, 2) == 2
            plotEdgesVectorized(vertices, boundaryFaces, ...
                struct('EdgeColor', options.BoundaryColor, ...
                'EdgeWidth', options.EdgeWidth*1.5), nDims);
        end

    case 3
        % Boundary faces in 3D
        if size(boundaryFaces, 2) == 3
            % Use trisurf to plot boundary triangles as wireframe
            trimesh(boundaryFaces, vertices(:, 1), vertices(:, 2), ...
                vertices(:, 3), 'EdgeColor', options.BoundaryColor, ...
                'LineWidth', options.EdgeWidth*1.5, 'FaceColor', 'none');
        else
            % General boundary faces
            patch('Faces', boundaryFaces, 'Vertices', vertices, ...
                'FaceColor', 'none', 'EdgeColor', options.BoundaryColor, ...
                'LineWidth', options.EdgeWidth*1.5);
        end
end
end

function plotFacesSlow(mesh, faceColor, alpha)
% PLOTFACES Plot the faces of a mesh.
%
%   plotFacesSlow(mesh, faceColor, alpha) renders mesh faces as filled polygons
%   with the specified color and transparency.
%
% Inputs:
%   mesh - approx.mesh.Mesh object
%   faceColor - Color for faces
%   alpha - Transparency level (0-1)

nDims = mesh.nDims;

if nDims <= 2
    % In 1D and 2D, faces are lower-dimensional and handled differently
    return;
end

% In 3D, faces are 2D polygons
vertices = mesh.vertices;
faces = mesh.faces;

for iFace = 1:mesh.nFaces
    faceVertices = vertices(faces(iFace, :), :);

    % Handle different face types
    nFaceVertices = size(faceVertices, 1);

    if nFaceVertices >= 3
        try
            % Use convhull for proper ordering if needed
            if nFaceVertices > 3
                k = convhull(faceVertices(:, 1), faceVertices(:, 2));
                faceVertices = faceVertices(k(1:end-1), :);
            end

            fill3(faceVertices(:, 1), faceVertices(:, 2), faceVertices(:, 3), ...
                faceColor, 'FaceAlpha', alpha, 'EdgeColor', 'none');
        catch ME
            % Skip problematic faces
            warning('approx:mesh:MeshPlotter:FacePlotFailed', ...
                'Failed to plot face %d: %s', iFace, ME.message);
        end
    end
end
end

function plotElementsSlow(mesh, elementColor, alpha)
% PLOTELEMENTS Plot the elements of a mesh.
%
%   plotElementsSlow(mesh, elementColor, alpha) renders mesh elements as filled
%   shapes with the specified color and transparency.
%
% Inputs:
%   mesh - approx.mesh.Mesh object
%   elementColor - Color for elements
%   alpha - Transparency level (0-1)

nDims = mesh.nDims;
vertices = mesh.vertices;
elements = mesh.elements;

for iElement = 1:mesh.nElements
    elementVertices = vertices(elements(iElement, :), :);

    try
        switch nDims
            case 1
                % Line segments
                plot(elementVertices(:, 1), zeros(size(elementVertices, 1), 1), ...
                    '-', 'Color', elementColor, 'LineWidth', 3);

            case 2
                % Polygons in 2D
                if size(elementVertices, 1) >= 3
                    fill(elementVertices(:, 1), elementVertices(:, 2), ...
                        elementColor, 'FaceAlpha', alpha, 'EdgeColor', 'none');
                end

            case 3
                % 3D elements (tetrahedra, etc.)
                if size(elementVertices, 1) >= 4
                    % Use convex hull for 3D visualization
                    k = convhulln(elementVertices);
                    patch('Faces', k, 'Vertices', elementVertices, ...
                        'FaceColor', elementColor, 'FaceAlpha', alpha, ...
                        'EdgeColor', 'none');
                end
        end
    catch ME
        % Skip problematic elements
        core.except.verify(0, 'ElementPlotFailed', ...
            'Failed to plot element %d: %s', iElement, ME.message);
    end
end
end

function plotBoundarySlow(mesh, boundaryColor, lineWidth)
% PLOTBOUNDARY Plot the boundary faces of a mesh.
%
%   plotBoundarySlow(mesh, boundaryColor, lineWidth) renders boundary faces as
%   lines or highlighted regions.
%
% Inputs:
%   mesh - approx.mesh.Mesh object
%   boundaryColor - Color for boundary
%   lineWidth - Width of boundary lines

nDims = mesh.nDims;
vertices = mesh.vertices;
faces = mesh.faces;
boundary = mesh.boundary;

for iBoundary = 1:length(boundary)
    iFace = boundary(iBoundary);
    faceVertices = vertices(faces(iFace, :), :);

    try
        switch nDims
            case 1
                % Boundary points in 1D
                scatter(faceVertices, zeros(size(faceVertices)), ...
                    100, boundaryColor, 'filled', 'MarkerEdgeColor', 'k');

            case 2
                % Boundary edges in 2D
                plot(faceVertices(:, 1), faceVertices(:, 2), '-', ...
                    'Color', boundaryColor, 'LineWidth', lineWidth);

            case 3
                % Boundary faces in 3D
                if size(faceVertices, 1) >= 3
                    % Plot face outline
                    faceVertices = [faceVertices; faceVertices(1, :)]; % Close the loop
                    plot3(faceVertices(:, 1), faceVertices(:, 2), ...
                        faceVertices(:, 3), '-', 'Color', boundaryColor, ...
                        'LineWidth', lineWidth);
                end
        end
    catch ME
        % Skip problematic boundary faces
        core.except.verify(0, 'BoundaryPlotFailed', ...
            'Failed to plot boundary face %d: %s', iFace, ME.message);
    end
end
end

function addLegend(options)
% ADDLEGEND Add legend to the plot based on visible elements.
%
%   addLegend(options) creates a legend showing only the mesh components
%   that are currently visible based on the plotting options.
%
% Inputs:
%   options - Structure containing plotting options

legendItems = {};
legendLabels = {};

% Add legend entries for visible components
if options.ShowVertices && ~isequal(options.VertexColor, 'none')
    legendItems{end+1} = scatter(NaN, NaN, options.VertexSize, ...
        options.VertexColor, 'filled');
    legendLabels{end+1} = 'Vertices';
end

if options.ShowEdges && ~isequal(options.EdgeColor, 'none')
    legendItems{end+1} = plot(NaN, NaN, '-', 'Color', options.EdgeColor, ...
        'LineWidth', options.EdgeWidth);
    legendLabels{end+1} = 'Edges';
end

if options.ShowFaces && ~isequal(options.FaceColor, 'none')
    legendItems{end+1} = fill(NaN, NaN, options.FaceColor, ...
        'FaceAlpha', options.Alpha, 'EdgeColor', 'none');
    legendLabels{end+1} = 'Faces';
end

if options.ShowElements && ~isequal(options.ElementColor, 'none')
    legendItems{end+1} = fill(NaN, NaN, options.ElementColor, ...
        'FaceAlpha', options.Alpha, 'EdgeColor', 'none');
    legendLabels{end+1} = 'Elements';
end

if options.ShowBoundary && ~isequal(options.BoundaryColor, 'none')
    legendItems{end+1} = plot(NaN, NaN, '-', 'Color', options.BoundaryColor, ...
        'LineWidth', options.EdgeWidth*1.5);
    legendLabels{end+1} = 'Boundary';
end

% Create legend if there are items to show
if ~isempty(legendItems)
    legend([legendItems{:}], legendLabels);
end
end
