%% Graph Visualization Examples
%
% This script demonstrates comprehensive usage of the GraphPlotter class
% to visualize various types of Graph objects in 1D, 2D, and 3D spaces.
% Examples include basic graphs, custom styling, complex geometries, and
% interactive visualization techniques.
%
% See also:
%   approx.mesh.Graph, approx.mesh.GraphPlotter

% Clear workspace and close figures
clc; clear; close all;

%% Basic 1D Graph Example

% Create simple 1D graph with evenly spaced vertices
vertices1D = (0:5)';  % Column vector of x-coordinates
edges1D = [1, 2; 2, 3; 3, 4; 4, 5; 5, 6]; % Connect adjacent vertices

% Create Graph object
basicGraph1D = approx.mesh.Graph(vertices1D, edges1D);

% Visualize basic 1D graph
figure('Name', 'Basic 1D Graph');
plotter = approx.mesh.GraphPlotter();
plotter.plot(basicGraph1D);
title('Basic 1D Graph with Uniform Spacing');
xlabel('X Position');

%% Custom Styled 1D Graph

% Create 1D graph with non-uniform vertex spacing
vertices1DCustom = [0; 0.5; 1.5; 3; 5; 8.5; 10];
edges1DCustom = [1, 2; 2, 3; 3, 4; 4, 5; 5, 6; 6, 7];

% Create Graph object
customGraph1D = approx.mesh.Graph(vertices1DCustom, edges1DCustom);

% Plot with custom styling options
figure('Name', 'Custom Styled 1D Graph');
plotter.plot(customGraph1D, ...
    'VertexColor', 'red', ...
    'VertexSize', 50, ...
    'EdgeColor', 'blue', ...
    'EdgeWidth', 2);
title('1D Graph with Custom Styling and Non-uniform Spacing');
xlabel('X Position');

%% 1D Graph with Disconnected Components

% Create 1D graph with gaps between connected components
vertices1DGaps = [0; 1; 2; 3; 5; 6; 7; 9; 10];
edges1DGaps = [1, 2; 2, 3; 3, 4; 5, 6; 6, 7; 8, 9];

% Create Graph object with disconnected segments
gapGraph1D = approx.mesh.Graph(vertices1DGaps, edges1DGaps);

% Visualize graph with disconnected components
figure('Name', '1D Graph with Disconnected Components');
plotter.plot(gapGraph1D, ...
    'VertexColor', 'green', ...
    'VertexSize', 36);
title('1D Graph with Disconnected Segments');
xlabel('X Position');

%% Comparison of Multiple 1D Graph Types

figure('Name', 'Comparison of 1D Graph Types');

% Uniform spacing subplot
subplot(3, 1, 1);
plotter.plot(basicGraph1D, ...
    'VertexColor', 'blue', ...
    'VertexSize', 30);
title('Uniform Vertex Spacing');
xlabel('X Position');

% Non-uniform spacing subplot
subplot(3, 1, 2);
plotter.plot(customGraph1D, ...
    'VertexColor', 'red', ...
    'VertexSize', 30);
title('Non-uniform Vertex Spacing');
xlabel('X Position');

% Disconnected components subplot
subplot(3, 1, 3);
plotter.plot(gapGraph1D, ...
    'VertexColor', 'green', ...
    'VertexSize', 30);
title('Disconnected Components');
xlabel('X Position');

%% Enhanced 1D Graph with Custom Annotations

figure('Name', '1D Graph with Custom Labels and Annotations');

% Plot base graph
plotter.plot(customGraph1D, ...
    'VertexColor', 'black', ...
    'VertexSize', 40, ...
    'EdgeColor', [0.6, 0.6, 0.8], ...
    'EdgeWidth', 2);

% Add vertex labels
nVertices1D = size(vertices1DCustom, 1);
for iVertex = 1:nVertices1D
    text(vertices1DCustom(iVertex), 0.1, ...
        sprintf('V%d', iVertex), ...
        'FontSize', 8, ...
        'HorizontalAlignment', 'center');
end

% Add edge distance annotations
nEdges1D = size(edges1DCustom, 1);
for iEdge = 1:nEdges1D
    vertex1Index = edges1DCustom(iEdge, 1);
    vertex2Index = edges1DCustom(iEdge, 2);
    midpointX = (vertices1DCustom(vertex1Index) + vertices1DCustom(vertex2Index)) / 2;
    edgeLength = abs(vertices1DCustom(vertex2Index) - vertices1DCustom(vertex1Index));
    text(midpointX, -0.1, ...
        sprintf('d=%.1f', edgeLength), ...
        'FontSize', 7, ...
        'HorizontalAlignment', 'center');
end

% Enhance plot appearance
title('1D Graph with Vertex Labels and Edge Distances');
xlabel('Position');
grid('minor');

%% Basic 2D Graph Examples

% Create simple square graph
vertices2DSquare = [0, 0; 1, 0; 1, 1; 0, 1];
edges2DSquare = [1, 2; 2, 3; 3, 4; 4, 1];
squareGraph = approx.mesh.Graph(vertices2DSquare, edges2DSquare);

% Plot basic square
figure('Name', 'Basic 2D Square Graph');
plotter.plot(squareGraph);
title('Simple Square Graph');

% Plot square with custom styling
figure('Name', 'Custom Styled 2D Square Graph');
plotter.plot(squareGraph, ...
    'VertexColor', 'red', ...
    'VertexSize', 50, ...
    'EdgeColor', 'blue', ...
    'EdgeWidth', 2);
title('Square Graph with Custom Styling');

%% Complex 2D Grid Graph

% Generate grid coordinates
gridSize = 5;
gridSpacing = 0.5;
[gridX, gridY] = meshgrid(0:gridSpacing:gridSize-gridSpacing, ...
                          0:gridSpacing:gridSize-gridSpacing);
vertices2DGrid = [gridX(:), gridY(:)];

% Calculate grid dimensions
nRowsGrid = size(gridX, 1);
nColsGrid = size(gridX, 2);

% Generate horizontal edges
edgesHorizontal = [];
for iRow = 1:nRowsGrid
    for jCol = 1:(nColsGrid-1)
        vertex1Index = (iRow-1)*nColsGrid + jCol;
        vertex2Index = vertex1Index + 1;
        edgesHorizontal = [edgesHorizontal; vertex1Index, vertex2Index];
    end
end

% Generate vertical edges
edgesVertical = [];
for iRow = 1:(nRowsGrid-1)
    for jCol = 1:nColsGrid
        vertex1Index = (iRow-1)*nColsGrid + jCol;
        vertex2Index = vertex1Index + nColsGrid;
        edgesVertical = [edgesVertical; vertex1Index, vertex2Index];
    end
end

% Combine all edges
edges2DGrid = [edgesHorizontal; edgesVertical];
gridGraph2D = approx.mesh.Graph(vertices2DGrid, edges2DGrid);

% Visualize 2D grid
figure('Name', '2D Grid Graph');
plotter.plot(gridGraph2D, 'VertexSize', 25);
title(sprintf('2D Grid Graph (%dx%d)', nRowsGrid, nColsGrid));

%% 3D Graph Examples

% Create cube vertices and edges
vertices3DCube = [
    0, 0, 0;  % Vertex 1
    1, 0, 0;  % Vertex 2
    1, 1, 0;  % Vertex 3
    0, 1, 0;  % Vertex 4
    0, 0, 1;  % Vertex 5
    1, 0, 1;  % Vertex 6
    1, 1, 1;  % Vertex 7
    0, 1, 1   % Vertex 8
];

edges3DCube = [
    % Bottom face edges
    1, 2; 2, 3; 3, 4; 4, 1;
    % Top face edges
    5, 6; 6, 7; 7, 8; 8, 5;
    % Vertical edges
    1, 5; 2, 6; 3, 7; 4, 8
];

cubeGraph = approx.mesh.Graph(vertices3DCube, edges3DCube);

% Visualize 3D cube
figure('Name', '3D Cube Graph');
plotter.plot(cubeGraph, ...
    'VertexColor', 'green', ...
    'VertexSize', 50);
title('3D Cube Graph');
view(30, 30);

%% Random 3D Graph with Nearest Neighbor Connectivity

% Generate random vertices in 3D space
nVerticesRandom = 30;
vertices3DRandom = rand(nVerticesRandom, 3);

% Connect each vertex to its nearest neighbors
nNearestNeighbors = 3;
edges3DRandom = [];

for iVertex = 1:nVerticesRandom
    % Calculate distances to all other vertices
    distances = inf(nVerticesRandom, 1);
    for jVertex = 1:nVerticesRandom
        if iVertex ~= jVertex
            distances(jVertex) = norm(vertices3DRandom(iVertex, :) - ...
                                    vertices3DRandom(jVertex, :));
        end
    end
    
    % Find nearest neighbors
    [~, nearestIndices] = sort(distances);
    for kNeighbor = 1:nNearestNeighbors
        neighborIndex = nearestIndices(kNeighbor);
        % Avoid duplicate edges by ensuring consistent ordering
        if iVertex < neighborIndex
            edges3DRandom = [edges3DRandom; iVertex, neighborIndex];
        end
    end
end

randomGraph3D = approx.mesh.Graph(vertices3DRandom, edges3DRandom);

% Visualize random 3D graph
figure('Name', 'Random 3D Graph with Nearest Neighbors');
plotter.plot(randomGraph3D, ...
    'VertexColor', 'blue', ...
    'EdgeWidth', 0.5);
title(sprintf('Random 3D Graph - %d Nearest Neighbors', nNearestNeighbors));
view(40, 35);

%% Side-by-Side 3D Graph Comparison

figure('Name', 'Comparison of 3D Graphs');

% Cube graph subplot
subplot(1, 2, 1);
plotter.plot(cubeGraph, ...
    'VertexColor', 'red', ...
    'EdgeColor', 'black');
title('Structured Cube Graph');
view(30, 30);

% Random graph subplot
subplot(1, 2, 2);
plotter.plot(randomGraph3D, ...
    'VertexColor', [0, 0.7, 0.3], ...
    'EdgeColor', [0.7, 0.7, 0.7]);
title('Random 3D Graph');
view(40, 35);

%% Spherical Graph Example

% Generate vertices on sphere surface
nThetaPoints = 10;
nPhiPoints = 16;
thetaValues = linspace(0, pi, nThetaPoints);
phiValues = linspace(0, 2*pi, nPhiPoints);
[thetaMesh, phiMesh] = meshgrid(thetaValues, phiValues);

% Convert spherical to Cartesian coordinates
sphereX = sin(thetaMesh) .* cos(phiMesh);
sphereY = sin(thetaMesh) .* sin(phiMesh);
sphereZ = cos(thetaMesh);

% Reshape coordinate arrays to vertex list
vertices3DSphere = [sphereX(:), sphereY(:), sphereZ(:)];

% Create connectivity along latitude lines
edgesLatitude = [];
for iThetaPoint = 1:nThetaPoints
    for jPhiPoint = 1:(nPhiPoints-1)
        vertex1Index = (jPhiPoint-1)*nThetaPoints + iThetaPoint;
        vertex2Index = jPhiPoint*nThetaPoints + iThetaPoint;
        edgesLatitude = [edgesLatitude; vertex1Index, vertex2Index];
    end
    % Connect last point to first in each latitude circle
    vertex1Index = (nPhiPoints-1)*nThetaPoints + iThetaPoint;
    vertex2Index = iThetaPoint;
    edgesLatitude = [edgesLatitude; vertex1Index, vertex2Index];
end

% Create connectivity along longitude lines
edgesLongitude = [];
for iThetaPoint = 1:(nThetaPoints-1)
    for jPhiPoint = 1:nPhiPoints
        vertex1Index = (jPhiPoint-1)*nThetaPoints + iThetaPoint;
        vertex2Index = vertex1Index + 1;
        edgesLongitude = [edgesLongitude; vertex1Index, vertex2Index];
    end
end

% Combine all spherical edges
edges3DSphere = [edgesLatitude; edgesLongitude];
sphereGraph = approx.mesh.Graph(vertices3DSphere, edges3DSphere);

% Visualize spherical graph
figure('Name', 'Spherical Graph');
plotter.plot(sphereGraph, ...
    'VertexSize', 20, ...
    'EdgeColor', [0.5, 0.5, 0.8], ...
    'EdgeWidth', 1);
title('Sphere-like Graph Structure');
axis('equal', 'tight');
view(30, 20);

%% Visualization Options Demonstration

figure('Name', 'Different Visualization Options');

% Vertices only
subplot(2, 2, 1);
plotter.plot(gridGraph2D, 'EdgeColor', 'none');
title('Vertices Only');

% Edges only
subplot(2, 2, 2);
plotter.plot(gridGraph2D, 'VertexColor', 'none');
title('Edges Only');

% Color-coded vertices by distance from origin
subplot(2, 2, 3);
distancesFromOrigin = sqrt(sum(gridGraph2D.Vertices.^2, 2));
scatter(gridGraph2D.Vertices(:, 1), gridGraph2D.Vertices(:, 2), ...
        50, distancesFromOrigin, 'filled');
hold('on');
plotter.plot(gridGraph2D, ...
    'VertexColor', 'none', ...
    'EdgeColor', [0.7, 0.7, 0.7]);
title('Color-coded by Distance from Origin');
colorbar;

% Custom edge styling based on edge length
subplot(2, 2, 4);
hold('on');
nEdgesGrid = size(gridGraph2D.Edges, 1);
for iEdge = 1:nEdgesGrid
    vertex1Index = gridGraph2D.Edges(iEdge, 1);
    vertex2Index = gridGraph2D.Edges(iEdge, 2);
    edgeStartX = gridGraph2D.Vertices(vertex1Index, 1);
    edgeEndX = gridGraph2D.Vertices(vertex2Index, 1);
    edgeStartY = gridGraph2D.Vertices(vertex1Index, 2);
    edgeEndY = gridGraph2D.Vertices(vertex2Index, 2);
    
    % Calculate edge length and apply color coding
    edgeLength = norm(gridGraph2D.Vertices(vertex1Index, :) - ...
                     gridGraph2D.Vertices(vertex2Index, :));
    plot([edgeStartX, edgeEndX], [edgeStartY, edgeEndY], '-', ...
        'Color', [0, 0.5, edgeLength*2], ...
        'LineWidth', 1 + edgeLength*3);
end
scatter(gridGraph2D.Vertices(:, 1), gridGraph2D.Vertices(:, 2), ...
        25, 'black', 'filled');
title('Custom Edge Styling by Length');
axis('equal');

%% Interactive 3D Graph Exploration

figure('Name', 'Interactive 3D Graph Exploration', ...
       'Position', [100, 100, 800, 600]);
plotter.plot(randomGraph3D, ...
    'VertexColor', 'red', ...
    'VertexSize', 40, ...
    'EdgeColor', 'blue', ...
    'EdgeWidth', 1.5);
title('Interactive 3D Graph - Use Mouse to Rotate, Pan, and Zoom');
view(40, 35);
axis('equal');
grid('on');
rotate3d('on');

% Display usage instructions
fprintf('\n=== Interactive Graph Exploration ===\n');
fprintf('The "Interactive 3D Graph Exploration" figure supports:\n');
fprintf('- Click and drag to rotate the view\n');
fprintf('- Scroll to zoom in/out\n');
fprintf('- Right-click and drag to pan\n');
fprintf('- Double-click to reset view\n\n');

fprintf('=== Graph Statistics Summary ===\n');
fprintf('1D Basic Graph: %d vertices, %d edges\n', ...
        basicGraph1D.NVertices, basicGraph1D.NEdges);
fprintf('2D Grid Graph: %d vertices, %d edges\n', ...
        gridGraph2D.NVertices, gridGraph2D.NEdges);
fprintf('3D Cube Graph: %d vertices, %d edges\n', ...
        cubeGraph.NVertices, cubeGraph.NEdges);
fprintf('3D Random Graph: %d vertices, %d edges\n', ...
        randomGraph3D.NVertices, randomGraph3D.NEdges);
fprintf('3D Sphere Graph: %d vertices, %d edges\n', ...
        sphereGraph.NVertices, sphereGraph.NEdges);