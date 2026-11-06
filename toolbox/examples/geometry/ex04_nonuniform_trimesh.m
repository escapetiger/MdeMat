function result = ex04_nonuniform_trimesh(options)
%EX04_NONUNIFORM_TRIMESH Generate and inspect a bundled DistMesh circle mesh.
%   R = ex04_nonuniform_trimesh(ElementSize=0.25, MakePlots=false, Seed=0)
%   returns original/refined triangulations, physical area and mesh quality.
%   The complete default study uses ElementSize=0.1 and visible figures.
%   Upstream progress rendering uses an owned figure. Headless mode hides and
%   deletes it; path, RNG and the caller current figure are restored on errors.
arguments
    options.ElementSize (1, 1) double {mustBeFinite, mustBePositive} = 0.1
    options.MakePlots (1, 1) logical = true
    options.Seed (1, 1) double {mustBeFinite, mustBeInteger, mustBeNonnegative, ...
        mustBeLessThanOrEqual(options.Seed, 4294967295)} = 0
end
originalPath = path;
originalRng = rng;
originalFigure = get(groot, 'CurrentFigure');
restoreCaller = onCleanup(@() restoreContext(originalPath, originalRng, originalFigure));
vendorRoot = fullfile(mdemat.tools.environment.toolboxRoot(), 'third_party', 'distmesh', 'src');
if ~isfolder(vendorRoot)
    vendorRoot = fullfile(fileparts(mdemat.tools.environment.toolboxRoot()), 'third_party', 'distmesh', 'src');
end
if ~isfile(fullfile(vendorRoot, 'distmesh2d.m'))
    error('mdemat:examples:DistMeshUnavailable', 'The bundled DistMesh source is missing.');
end
addpath(vendorRoot, '-begin');
rng(options.Seed, 'twister');
visibility = 'off';
if options.MakePlots, visibility = get(groot, 'DefaultFigureVisible'); end
progress = figure('Name', 'DistMesh Mesh Generation', 'Visible', visibility);
if ~options.MakePlots
    closeProgress = onCleanup(@() deleteOwned(progress));
end

fprintf('Starting DistMesh Triangulation Example...\n');

%% Define domain and mesh parameters
fprintf('Defining domain and mesh parameters...\n');

% Circle domain parameters
radius = 1.0;
center = [0.0, 0.0];
elementSize = options.ElementSize;

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
simplexMesh = mdemat.data.mesh.Triangulation(meshVertices, meshElements);

fprintf('Triangulation created successfully.\n');
fprintf('Mesh properties:\n');
fprintf('  Dimensions: %d\n', simplexMesh.NDims);
fprintf('  Vertices: %d\n', simplexMesh.NVertices);
fprintf('  Elements: %d\n', simplexMesh.NElements);
fprintf('  Faces: %d\n', simplexMesh.NFaces);
fprintf('  Boundary faces: %d\n', simplexMesh.NBoundaryFaces);
    

%% Visualize the original mesh
if options.MakePlots
    fprintf('Visualizing the original mesh...\n');

    % Draw through the mesh method
    % Create figure for original mesh
    figure('Name', 'DistMesh Generated Triangulation', 'Position', [100, 100, 800, 600]);

    % Plot original mesh with different views
    subplot(2, 2, 1);
    plot(simplexMesh, 'showElements', true, 'elementColor', 'cyan', ...
        'showBoundary', true, 'boundaryColor', 'red', 'Alpha', 0.3);
    title('Original Mesh: Elements and Boundary');

    subplot(2, 2, 2);
    plot(simplexMesh, 'showVertices', true, 'vertexColor', 'blue', ...
        'showEdges', true, 'edgeColor', 'black', 'vertexSize', 25);
    title('Original Mesh: Vertices and Edges');

    subplot(2, 2, 3);
    plot(simplexMesh, 'showFaces', true, 'faceColor', 'yellow', ...
        'showBoundary', true, 'boundaryColor', 'red', 'Alpha', 0.4);
    title('Original Mesh: Faces and Boundary');

    subplot(2, 2, 4);
    plot(simplexMesh, 'showElements', true, 'elementColor', 'green', ...
        'showVertices', true, 'vertexColor', 'red', 'vertexSize', 30, ...
        'Alpha', 0.5);
    title('Original Mesh: Combined View');

end

%% Demonstrate mesh refinement
fprintf('Demonstrating mesh refinement...\n');


% Refine the mesh once
refinedMesh = simplexMesh.refine(1);

fprintf('Mesh refinement completed.\n');
fprintf('Refined mesh properties:\n');
fprintf('  Vertices: %d (was %d)\n', refinedMesh.NVertices, simplexMesh.NVertices);
fprintf('  Elements: %d (was %d)\n', refinedMesh.NElements, simplexMesh.NElements);
fprintf('  Faces: %d (was %d)\n', refinedMesh.NFaces, simplexMesh.NFaces);

if options.MakePlots
    % Create figure for refined mesh comparison
    figure('Name', 'Mesh Refinement Comparison', 'Position', [150, 150, 1000, 400]);

    % Plot original mesh
    subplot(1, 2, 1);
    plot(simplexMesh, 'showElements', true, 'elementColor', 'cyan', ...
        'showBoundary', true, 'boundaryColor', 'red', 'Alpha', 0.3);
    title(sprintf('Original Mesh (%d elements)', simplexMesh.NElements));

    % Plot refined mesh
    subplot(1, 2, 2);
    plot(refinedMesh, 'showElements', true, 'elementColor', 'blue', ...
        'showBoundary', true, 'boundaryColor', 'red', 'Alpha', 0.3);
    title(sprintf('Refined Mesh (%d elements)', refinedMesh.NElements));


end

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

if options.MakePlots
    % Visualize the transformation
    figure('Name', 'Coordinate Transformation Demo', 'Position', [200, 200, 600, 500]);

    % Plot mesh
    plot(simplexMesh, 'showElements', true, 'elementColor', 'yellow', ...
        'showBoundary', true, 'boundaryColor', 'black', 'Alpha', 0.2);
    hold on;

    % Highlight target elements
    for iElement = targetElements
        elementVertices = simplexMesh.Vertices(simplexMesh.Elements(iElement, :), :);
        fill(elementVertices(:, 1), elementVertices(:, 2), 'yellow', ...
            'FaceAlpha', 0.5, 'edgeColor', 'red', 'LineWidth', 2);
    end

    % Plot transformed points
    nTargetElements = length(targetElements);
    reshapedPoints = reshape(physicalPoints, simplexMesh.NDims, nTestPoints, nTargetElements);

    for iElement = 1:nTargetElements
        elementPoints = reshapedPoints(:, :, iElement);
        scatter(elementPoints(1, :), elementPoints(2, :), 10, 'g', ...
            'filled', 'MarkeredgeColor', 'black');
    end

    title('Coordinate Transformation: Reference to Physical');
    legend('Mesh Elements', 'Target Elements', 'Transformed Points', ...
        'Location', 'best');
    

end

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

result = struct('schemaVersion', 1, 'mesh', simplexMesh, 'refined', refinedMesh, ...
    'physicalArea', sum(simplexMesh.computeElementJacobianDeterminants())/2, ...
    'quality', meshQuality, 'elementSize', elementSize, 'seed', options.Seed, ...
    'referencePoints', referencePoints, 'physicalPoints', physicalPoints, ...
    'targetElements', targetElements);
fprintf('DistMesh Triangulation Example completed successfully!\n');



end

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
function restoreContext(originalPath, originalRng, originalFigure)
path(originalPath);
rng(originalRng);
if isempty(originalFigure) || isgraphics(originalFigure, 'figure')
    set(groot, 'CurrentFigure', originalFigure);
end
end

function deleteOwned(figureHandle)
if isgraphics(figureHandle, 'figure'), delete(figureHandle); end
end
