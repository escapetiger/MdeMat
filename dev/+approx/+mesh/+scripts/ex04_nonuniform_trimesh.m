clear; close all; clc;

fprintf('Starting DistMesh Triangulation Example...\n');
addpath('distmesh');

%% Define domain and mesh parameters
fprintf('Defining domain and mesh parameters...\n');

% Circle domain parameters
radius = 1.0;
center = [0.0, 0.0];
elementSize = 0.1;

% Define distance function for circle
distanceFunction = @(p) dcircle(p, center(1), center(2), radius);

% Define boundary constraints (none for simple circle)
boundaryConstraints = [];

% Initial point distribution for better mesh quality
boundingBox = [-radius, -radius; radius, radius];

%% Generate mesh using distmesh
fprintf('Generating mesh using ..\n');

[meshVertices, meshElements] = distmesh2d(distanceFunction, ...
    @huniform, elementSize, boundingBox, boundaryConstraints);

fprintf('Mesh generation completed successfully.\n');
fprintf('Generated %d vertices and %d elements.\n', ...
    size(meshVertices, 1), size(meshElements, 1));

%% Create Triangulation object
fprintf('Creating Triangulation object...\n');

% Create Triangulation from distmesh output
simplexMesh = approx.mesh.Triangulation(meshVertices, meshElements);

fprintf('Triangulation created successfully.\n');
fprintf('Mesh properties:\n');
fprintf('  Dimensions: %d\n', simplexMesh.NDims);
fprintf('  Vertices: %d\n', simplexMesh.NVertices);
fprintf('  Elements: %d\n', simplexMesh.NElements);
fprintf('  Faces: %d\n', simplexMesh.NFaces);
fprintf('  Boundary faces: %d\n', simplexMesh.NBoundaryFaces);
    

%% Visualize the original mesh
fprintf('Visualizing the original mesh...\n');

% Create mesh plotter
meshPlotter = approx.mesh.PolytopalMeshPlotter();

% Create figure for original mesh
figure('Name', 'DistMesh Generated Triangulation', 'Position', [100, 100, 800, 600]);

% Plot original mesh with different views
subplot(2, 2, 1);
meshPlotter.plot(simplexMesh, 'ShowElements', true, 'ElementColor', 'cyan', ...
    'ShowBoundary', true, 'BoundaryColor', 'red', 'Alpha', 0.3);
title('Original Mesh: Elements and Boundary');

subplot(2, 2, 2);
meshPlotter.plot(simplexMesh, 'ShowVertices', true, 'VertexColor', 'blue', ...
    'ShowEdges', true, 'EdgeColor', 'black', 'VertexSize', 25);
title('Original Mesh: Vertices and Edges');

subplot(2, 2, 3);
meshPlotter.plot(simplexMesh, 'ShowFaces', true, 'FaceColor', 'yellow', ...
    'ShowBoundary', true, 'BoundaryColor', 'red', 'Alpha', 0.4);
title('Original Mesh: Faces and Boundary');

subplot(2, 2, 4);
meshPlotter.plot(simplexMesh, 'ShowElements', true, 'ElementColor', 'green', ...
    'ShowVertices', true, 'VertexColor', 'red', 'VertexSize', 30, ...
    'Alpha', 0.5);
title('Original Mesh: Combined View');

%% Demonstrate mesh refinement
fprintf('Demonstrating mesh refinement...\n');


% Refine the mesh once
refinedMesh = simplexMesh.refine(1);

fprintf('Mesh refinement completed.\n');
fprintf('Refined mesh properties:\n');
fprintf('  Vertices: %d (was %d)\n', refinedMesh.NVertices, simplexMesh.NVertices);
fprintf('  Elements: %d (was %d)\n', refinedMesh.NElements, simplexMesh.NElements);
fprintf('  Faces: %d (was %d)\n', refinedMesh.NFaces, simplexMesh.NFaces);

% Create figure for refined mesh comparison
figure('Name', 'Mesh Refinement Comparison', 'Position', [150, 150, 1000, 400]);

% Plot original mesh
subplot(1, 2, 1);
meshPlotter.plot(simplexMesh, 'ShowElements', true, 'ElementColor', 'cyan', ...
    'ShowBoundary', true, 'BoundaryColor', 'red', 'Alpha', 0.3);
title(sprintf('Original Mesh (%d elements)', simplexMesh.NElements));

% Plot refined mesh
subplot(1, 2, 2);
meshPlotter.plot(refinedMesh, 'ShowElements', true, 'ElementColor', 'blue', ...
    'ShowBoundary', true, 'BoundaryColor', 'red', 'Alpha', 0.3);
title(sprintf('Refined Mesh (%d elements)', refinedMesh.NElements));


%% Demonstrate coordinate transformation
fprintf('Demonstrating coordinate transformation...\n');


% Define reference points in Cartesian coordinates
nTestPoints = 4;
referencePoints = [0.5 0 0.5 1/3; 0 0.5 0.5 1/3];

% Select a few elements for transformation
targetElements = 1:min(3, simplexMesh.NElements);

% Transform reference points to physical coordinates
physicalPoints = simplexMesh.collocate(referencePoints, targetElements);

fprintf('Coordinate transformation completed.\n');
fprintf('Transformed %d reference points for %d elements.\n', ...
    nTestPoints, length(targetElements));

% Visualize the transformation
figure('Name', 'Coordinate Transformation Demo', 'Position', [200, 200, 600, 500]);

% Plot mesh
meshPlotter.plot(simplexMesh, 'ShowElements', true, 'ElementColor', 'yellow', ...
    'ShowBoundary', true, 'BoundaryColor', 'black', 'Alpha', 0.2);
hold on;

% Highlight target elements
for iElement = targetElements
    elementVertices = simplexMesh.Vertices(simplexMesh.Elements(iElement, :), :);
    fill(elementVertices(:, 1), elementVertices(:, 2), 'yellow', ...
        'FaceAlpha', 0.5, 'EdgeColor', 'red', 'LineWidth', 2);
end

% Plot transformed points
nTargetElements = length(targetElements);
reshapedPoints = reshape(physicalPoints, simplexMesh.NDims, nTestPoints, nTargetElements);

for iElement = 1:nTargetElements
    elementPoints = reshapedPoints(:, :, iElement);
    scatter(elementPoints(1, :), elementPoints(2, :), 10, 'g', ...
        'filled', 'MarkerEdgeColor', 'black');
end

title('Coordinate Transformation: Reference to Physical');
legend('Mesh Elements', 'Target Elements', 'Transformed Points', ...
    'Location', 'best');
    

%% Display mesh quality information
fprintf('Analyzing mesh quality...\n');

meshQuality = analyzeMeshQuality(simplexMesh);

fprintf('Mesh Quality Analysis:\n');
fprintf('  Average element area: %.6f\n', meshQuality.avgArea);
fprintf('  Min element area: %.6f\n', meshQuality.minArea);
fprintf('  Max element area: %.6f\n', meshQuality.maxArea);
fprintf('  Area ratio (max/min): %.2f\n', meshQuality.areaRatio);
fprintf('  Average element quality: %.4f\n', meshQuality.avgQuality);
fprintf('  Min element quality: %.4f\n', meshQuality.minQuality);

fprintf('DistMesh Triangulation Example completed successfully!\n');
rmpath('distmesh');

function quality = analyzeMeshQuality(mesh)
% ANALYZEMESHQUALITY Compute basic mesh quality metrics.
%
%   quality = analyzeMeshQuality(mesh) computes various quality metrics
%   for a triangular mesh including element areas and shape quality.

nElements = mesh.NElements;
areas = zeros(nElements, 1);
qualities = zeros(nElements, 1);

vertices = mesh.Vertices;
elements = mesh.Elements;

for iElement = 1:nElements
    % Get element vertices
    elementVertices = vertices(elements(iElement, :), :);
    
    if size(elementVertices, 1) == 3
        % Triangle area using cross product
        v1 = elementVertices(2, :) - elementVertices(1, :);
        v2 = elementVertices(3, :) - elementVertices(1, :);
        area = 0.5 * abs(v1(1)*v2(2) - v1(2)*v2(1));
        areas(iElement) = area;
        
        % Triangle quality (ratio of inscribed to circumscribed circle radii)
        % Higher values indicate better quality (equilateral triangle = 1)
        sideLength1 = norm(elementVertices(2, :) - elementVertices(1, :));
        sideLength2 = norm(elementVertices(3, :) - elementVertices(2, :));
        sideLength3 = norm(elementVertices(1, :) - elementVertices(3, :));
        
        semiPerimeter = (sideLength1 + sideLength2 + sideLength3) / 2;
        
        if semiPerimeter > 0 && area > 0
            inRadius = area / semiPerimeter;
            circumRadius = (sideLength1 * sideLength2 * sideLength3) / (4 * area);
            qualities(iElement) = 2 * inRadius / circumRadius;
        else
            qualities(iElement) = 0;
        end
    end
end

% Remove zero entries
validAreas = areas(areas > 0);
validQualities = qualities(qualities > 0);

quality.avgArea = mean(validAreas);
quality.minArea = min(validAreas);
quality.maxArea = max(validAreas);
quality.areaRatio = quality.maxArea / quality.minArea;
quality.avgQuality = mean(validQualities);
quality.minQuality = min(validQualities);
end