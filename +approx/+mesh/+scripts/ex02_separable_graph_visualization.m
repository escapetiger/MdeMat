%% SeparableGraph Visualization Examples
%
% This script demonstrates comprehensive usage of the SeparableGraphPlotter
% class to visualize SeparableGraph objects across different dimensions.
% Examples include memory-efficient visualization, subsampling techniques,
% projection methods for high-dimensional cases, and performance comparisons
% with full Graph representations.
%
% See also:
%   approx.mesh.SeparableGraph, approx.mesh.SeparableGraphPlotter,
%   approx.mesh.Graph, approx.mesh.GraphPlotter

% Clear workspace and close figures
clc; clear; close all;

%% Basic 1D SeparableGraph Example

% Create 1D separable graph with uniform spacing
nVertices1D = 20;
vertices1D = {linspace(0, 10, nVertices1D)};
edges1D = {[(1:nVertices1D-1)', (2:nVertices1D)']};
separableGraph1D = approx.mesh.SeparableGraph(vertices1D, edges1D);

% Visualize 1D separable graph
figure('Name', '1D SeparableGraph Visualization');
plotter = approx.mesh.SeparableGraphPlotter();
plotter.plot(separableGraph1D, ...
    'VertexColor', 'red', ...
    'VertexSize', 30);
title('1D SeparableGraph with Uniform Vertex Distribution');

% Display graph properties
fprintf('=== 1D SeparableGraph Properties ===\n');
fprintf('Dimensions: %d\n', separableGraph1D.NDims);
fprintf('Vertices per dimension: %s\n', mat2str(separableGraph1D.NVertices));
fprintf('Total vertices: %d\n', separableGraph1D.NTotalVertices);
fprintf('Total edges: %d\n', separableGraph1D.NTotalEdges);

%% 2D SeparableGraph Grid Example

% Create 2D separable graph representing a structured grid
nVerticesX = 10;
nVerticesY = 8;
vertices2D = {linspace(0, 5, nVerticesX), linspace(0, 4, nVerticesY)};
edges2D = {[(1:nVerticesX-1)', (2:nVerticesX)'], ...
           [(1:nVerticesY-1)', (2:nVerticesY)']};
separableGraph2D = approx.mesh.SeparableGraph(vertices2D, edges2D);

% Visualize 2D separable graph
figure('Name', '2D SeparableGraph Grid Visualization');
plotter.plot(separableGraph2D, ...
    'VertexColor', 'blue', ...
    'EdgeColor', [0.5, 0.5, 0.5], ...
    'VertexSize', 25);
title('2D SeparableGraph Grid Structure');

% Display graph properties
fprintf('\n=== 2D SeparableGraph Properties ===\n');
fprintf('Dimensions: %d\n', separableGraph2D.NDims);
fprintf('Vertices per dimension: %s\n', mat2str(separableGraph2D.NVertices));
fprintf('Total vertices: %d\n', separableGraph2D.NTotalVertices);
fprintf('Total edges: %d\n', separableGraph2D.NTotalEdges);

%% 3D SeparableGraph Cube Example

% Create 3D separable graph representing a structured cube
nVerticesCube = 4;
vertices3D = {linspace(0, 3, nVerticesCube), ...
              linspace(0, 3, nVerticesCube), ...
              linspace(0, 3, nVerticesCube)};
edges3D = {[(1:nVerticesCube-1)', (2:nVerticesCube)'], ...
           [(1:nVerticesCube-1)', (2:nVerticesCube)'], ...
           [(1:nVerticesCube-1)', (2:nVerticesCube)']};
separableGraph3D = approx.mesh.SeparableGraph(vertices3D, edges3D);

% Visualize 3D separable graph
figure('Name', '3D SeparableGraph Cube Visualization');
plotter.plot(separableGraph3D, ...
    'VertexSize', 40, ...
    'EdgeWidth', 1.5);
title('3D SeparableGraph Cube Structure');
view(45, 30);

% Display graph properties
fprintf('\n=== 3D SeparableGraph Properties ===\n');
fprintf('Dimensions: %d\n', separableGraph3D.NDims);
fprintf('Vertices per dimension: %s\n', mat2str(separableGraph3D.NVertices));
fprintf('Total vertices: %d\n', separableGraph3D.NTotalVertices);
fprintf('Total edges: %d\n', separableGraph3D.NTotalEdges);

%% Large 3D SeparableGraph with Subsampling

% Create larger 3D grid that demonstrates subsampling functionality
nVerticesLarge = 20;
vertices3DLarge = {linspace(0, 3, nVerticesLarge), ...
                   linspace(0, 3, nVerticesLarge), ...
                   linspace(0, 3, nVerticesLarge)};

% Generate edge connectivity for each dimension
edges3DLarge = cell(1, 3);
for iDim = 1:3
    edges3DLarge{iDim} = [(1:nVerticesLarge-1)', (2:nVerticesLarge)'];
end

separableGraph3DLarge = approx.mesh.SeparableGraph(vertices3DLarge, edges3DLarge);

% Visualize with automatic subsampling for large graphs
figure('Name', 'Large 3D SeparableGraph with Subsampling');
plotter.plot(separableGraph3DLarge, ...
    'VertexSize', 30, ...
    'MaxVertices', 1000);  % Trigger subsampling
title('Large 3D SeparableGraph (Automatically Subsampled)');

% Display large graph properties
fprintf('\n=== Large 3D SeparableGraph Properties ===\n');
fprintf('Dimensions: %d\n', separableGraph3DLarge.NDims);
fprintf('Vertices per dimension: %s\n', mat2str(separableGraph3DLarge.NVertices));
fprintf('Total vertices: %d\n', separableGraph3DLarge.NTotalVertices);
fprintf('Total edges: %d\n', separableGraph3DLarge.NTotalEdges);

%% High-Dimensional SeparableGraph with Projection

% Create 4D separable graph for projection demonstration
nVertices4D = 3;
vertices4D = {linspace(0, 2, nVertices4D), ...
              linspace(0, 2, nVertices4D), ...
              linspace(0, 2, nVertices4D), ...
              linspace(0, 2, nVertices4D)};

% Generate connectivity for 4D hypercube
edges4D = cell(1, 4);
for iDim = 1:4
    edges4D{iDim} = [1, 2; 2, 3];
end

separableGraph4D = approx.mesh.SeparableGraph(vertices4D, edges4D);

% Visualize with default projection (first 3 dimensions)
figure('Name', '4D SeparableGraph - Default Projection');
plotter.plot(separableGraph4D, ...
    'VertexColor', 'green', ...
    'VertexSize', 50);
title('4D SeparableGraph - Default Projection (Dims 1,2,3)');

% Visualize with alternative projection
figure('Name', '4D SeparableGraph - Alternative Projection');
plotter.plot(separableGraph4D, ...
    'VertexColor', 'magenta', ...
    'VertexSize', 50, ...
    'ProjectionDims', [2, 3, 4]);
title('4D SeparableGraph - Alternative Projection (Dims 2,3,4)');

% Display 4D graph properties
fprintf('\n=== 4D SeparableGraph Properties ===\n');
fprintf('Dimensions: %d\n', separableGraph4D.NDims);
fprintf('Vertices per dimension: %s\n', mat2str(separableGraph4D.NVertices));
fprintf('Total vertices: %d\n', separableGraph4D.NTotalVertices);
fprintf('Total edges: %d\n', separableGraph4D.NTotalEdges);

%% Comparison: SeparableGraph vs Full Graph Representation

% Convert 2D separable graph to full representation
fullGraph2D = separableGraph2D.full();

% Visualize both representations side by side
figure('Name', 'SeparableGraph vs Full Graph Comparison');

% SeparableGraph visualization
subplot(1, 2, 1);
plotter.plot(separableGraph2D, ...
    'VertexColor', 'blue', ...
    'EdgeColor', 'red');
title('SeparableGraph Representation');

% Full Graph visualization
subplot(1, 2, 2);
regularPlotter = approx.mesh.GraphPlotter();
regularPlotter.plot(fullGraph2D, ...
    'VertexColor', 'blue', ...
    'EdgeColor', 'red');
title('Full Graph Representation');

%% Memory Usage Analysis and Comparison

% Create medium-sized 3D grid for memory analysis
nVerticesMedium = 50;
vertices3DMedium = {linspace(0, 3, nVerticesMedium), ...
                    linspace(0, 3, nVerticesMedium), ...
                    linspace(0, 3, nVerticesMedium)};

% Generate connectivity for medium-sized graph
edges3DMedium = cell(1, 3);
for iDim = 1:3
    edges3DMedium{iDim} = [(1:nVerticesMedium-1)', (2:nVerticesMedium)'];
end

separableGraphMedium = approx.mesh.SeparableGraph(vertices3DMedium, edges3DMedium);

% Convert to full graph representation for comparison
fprintf('\n=== Converting SeparableGraph to Full Graph ===\n');
fprintf('Converting %dx%dx%d grid to full representation...\n', ...
        nVerticesMedium, nVerticesMedium, nVerticesMedium);
fullGraphMedium = separableGraphMedium.full();

% Analyze memory usage of both representations
separableInfo = whos('separableGraphMedium');
fullInfo = whos('fullGraphMedium');

% Analyze component memory usage
separableVerticesData = separableGraphMedium.Vertices;
separableEdgesData = separableGraphMedium.Edges;
fullVerticesData = fullGraphMedium.Vertices;
fullEdgesData = fullGraphMedium.Edges;

separableVerticesInfo = whos('separableVerticesData');
separableEdgesInfo = whos('separableEdgesData');
fullVerticesInfo = whos('fullVerticesData');
fullEdgesInfo = whos('fullEdgesData');

% Display detailed memory comparison
fprintf('\n=== Memory Usage Analysis ===\n');
fprintf('SeparableGraph object total: %.2f KB\n', separableInfo.bytes/1024);
fprintf('Full Graph object total: %.2f KB\n', fullInfo.bytes/1024);
fprintf('\nDetailed breakdown:\n');
fprintf('  SeparableGraph vertices (factored): %.2f KB\n', separableVerticesInfo.bytes/1024);
fprintf('  SeparableGraph edges (factored): %.2f KB\n', separableEdgesInfo.bytes/1024);
fprintf('  Full Graph vertices (expanded): %.2f KB\n', fullVerticesInfo.bytes/1024);
fprintf('  Full Graph edges (expanded): %.2f KB\n', fullEdgesInfo.bytes/1024);

% Calculate memory efficiency
separableDataTotal = separableVerticesInfo.bytes + separableEdgesInfo.bytes;
fullDataTotal = fullVerticesInfo.bytes + fullEdgesInfo.bytes;
memoryReductionFactor = fullDataTotal / separableDataTotal;

fprintf('\nMemory efficiency:\n');
fprintf('  SeparableGraph data total: %.2f KB\n', separableDataTotal/1024);
fprintf('  Full Graph data total: %.2f KB\n', fullDataTotal/1024);
fprintf('  Memory reduction factor: %.1fx\n', memoryReductionFactor);

%% Theoretical Memory Scaling Analysis

fprintf('\n=== Theoretical Scaling Analysis ===\n');

% Analyze theoretical memory requirements for larger grids
gridSizes = [10, 50, 100, 200];
nGridSizes = length(gridSizes);

fprintf('Grid Size | SeparableGraph (MB) | Full Graph (MB) | Reduction Factor\n');
fprintf('----------|---------------------|-----------------|------------------\n');

for iSize = 1:nGridSizes
    gridSize = gridSizes(iSize);
    
    % Calculate theoretical memory for SeparableGraph
    % Each dimension stores n vertices and (n-1) edges
    separableVerticesMemory = 3 * gridSize * 8; % 3 dims, 8 bytes per double
    separableEdgesMemory = 3 * (gridSize-1) * 2 * 8; % 3 dims, 2 indices per edge
    separableTotal = separableVerticesMemory + separableEdgesMemory;
    
    % Calculate theoretical memory for full Graph
    % Total vertices: n^3, each with 3 coordinates
    fullVerticesMemory = (gridSize^3) * 3 * 8;
    % Total edges: 3*n^2*(n-1), each with 2 vertex indices
    fullEdgesMemory = 3 * (gridSize^2) * (gridSize-1) * 2 * 8;
    fullTotal = fullVerticesMemory + fullEdgesMemory;
    
    % Calculate reduction factor
    reductionFactor = fullTotal / separableTotal;
    
    fprintf('%8d  |%19.2f |%15.2f |%16.1f\n', ...
            gridSize, separableTotal/(1024*1024), fullTotal/(1024*1024), reductionFactor);
end

%% Performance Visualization Demo

% Create different sized grids to demonstrate performance characteristics
performanceGridSizes = [5, 10, 15];
nPerformanceTests = length(performanceGridSizes);

figure('Name', 'SeparableGraph Performance Scaling');

for iTest = 1:nPerformanceTests
    gridSize = performanceGridSizes(iTest);
    
    % Create separable graph
    verticesPerf = {linspace(0, 1, gridSize), ...
                    linspace(0, 1, gridSize), ...
                    linspace(0, 1, gridSize)};
    edgesPerf = cell(1, 3);
    for iDim = 1:3
        edgesPerf{iDim} = [(1:gridSize-1)', (2:gridSize)'];
    end
    
    separableGraphPerf = approx.mesh.SeparableGraph(verticesPerf, edgesPerf);
    
    % Plot in subplot
    subplot(1, nPerformanceTests, iTest);
    plotter.plot(separableGraphPerf, ...
        'VertexSize', max(5, 30-gridSize), ...
        'EdgeWidth', max(0.5, 2-gridSize/10));
    title(sprintf('%dx%dx%d Grid\n%d vertices, %d edges', ...
          gridSize, gridSize, gridSize, ...
          separableGraphPerf.NTotalVertices, ...
          separableGraphPerf.NTotalEdges));
    
    if iTest == 1
        view(45, 30);
    elseif iTest == 2
        view(60, 25);
    else
        view(30, 40);
    end
end

%% Non-Uniform Grid Demonstration

% Create separable graph with non-uniform spacing
verticesNonUniform = {[0, 0.1, 0.5, 1.2, 2.5, 4.0], ...
                      [0, 0.3, 1.0, 2.0], ...
                      [0, 1.5, 3.0]};

% Generate connectivity maintaining separable structure
edgesNonUniform = cell(1, 3);
edgesNonUniform{1} = [(1:5)', (2:6)'];  % 5 edges in first dimension
edgesNonUniform{2} = [(1:3)', (2:4)'];  % 3 edges in second dimension
edgesNonUniform{3} = [1, 2; 2, 3];      % 2 edges in third dimension

separableGraphNonUniform = approx.mesh.SeparableGraph(verticesNonUniform, edgesNonUniform);

figure('Name', 'Non-Uniform SeparableGraph');
plotter.plot(separableGraphNonUniform, ...
    'VertexColor', [0.8, 0.2, 0.6], ...
    'VertexSize', 60, ...
    'EdgeColor', [0.2, 0.2, 0.8], ...
    'EdgeWidth', 2);
title('SeparableGraph with Non-Uniform Vertex Spacing');
view(50, 25);

%% Comparison with Different Visualization Techniques

% Create moderately sized graph for visualization comparison
nVerticesComparison = 8;
verticesComparison = {linspace(0, 2, nVerticesComparison), ...
                      linspace(0, 2, nVerticesComparison), ...
                      linspace(0, 2, nVerticesComparison)};
edgesComparison = cell(1, 3);
for iDim = 1:3
    edgesComparison{iDim} = [(1:nVerticesComparison-1)', (2:nVerticesComparison)'];
end
separableGraphComparison = approx.mesh.SeparableGraph(verticesComparison, edgesComparison);

figure('Name', 'SeparableGraph Visualization Techniques');

% Standard visualization
subplot(2, 2, 1);
plotter.plot(separableGraphComparison);
title('Standard Visualization');
view(45, 30);

% Vertices only
subplot(2, 2, 2);
plotter.plot(separableGraphComparison, ...
    'EdgeColor', 'none', ...
    'VertexSize', 40);
title('Vertices Only');
view(45, 30);

% Edges only
subplot(2, 2, 3);
plotter.plot(separableGraphComparison, ...
    'VertexColor', 'none', ...
    'EdgeWidth', 2);
title('Edges Only');
view(45, 30);

% Custom styling
subplot(2, 2, 4);
plotter.plot(separableGraphComparison, ...
    'VertexColor', [1, 0.5, 0], ...
    'VertexSize', 25, ...
    'EdgeColor', [0, 0.5, 1], ...
    'EdgeWidth', 1.5);
title('Custom Styling');
view(45, 30);

%% Summary and Recommendations

fprintf('\n=== SeparableGraph Advantages Summary ===\n');
fprintf('1. Memory Efficiency:\n');
fprintf('   - Scales as O(sum of factor sizes) vs O(product)\n');
fprintf('   - Significant savings for structured grids\n');
fprintf('   - Example: %dx%dx%d grid uses %.1fx less memory\n', ...
        nVerticesMedium, nVerticesMedium, nVerticesMedium, memoryReductionFactor);

fprintf('\n2. Computational Benefits:\n');
fprintf('   - Direct tensor-product operations\n');
fprintf('   - No need to store full vertex/edge arrays\n');
fprintf('   - Efficient subsampling for visualization\n');

fprintf('\n3. Recommended Use Cases:\n');
fprintf('   - Structured grids and meshes\n');
fprintf('   - High-dimensional tensor-product domains\n');
fprintf('   - Large-scale problems with separable structure\n');

fprintf('\n4. Visualization Features:\n');
fprintf('   - Automatic subsampling for large graphs\n');
fprintf('   - Projection for high-dimensional cases\n');
fprintf('   - Memory-efficient rendering\n');

fprintf('\n=== Performance Notes ===\n');
fprintf('- Use SeparableGraph for structured, tensor-product grids\n');
fprintf('- Use regular Graph for irregular connectivity\n');
fprintf('- Consider subsampling for visualization of large grids\n');
fprintf('- Project high-dimensional graphs to 3D or lower\n\n');