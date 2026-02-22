clear;
close all;
clc;

fprintf('Starting Mesh Validation Example...\n');
fprintf('=======================================\n');

%% Generate uniform 4x4 Triangulation on unit box
fprintf('Step 1: Generating uniform 4x4 Triangulation on [0,1]x[0,1]...\n');

[vertices, elements] = generateUniformBoxMesh(4, 4, [0, 1], [0, 1]);

fprintf('Generated mesh with %d vertices and %d elements.\n', ...
    size(vertices, 1), size(elements, 1));

% Create Triangulation object
mesh = approx.mesh.Triangulation(vertices, elements);

fprintf('Triangulation properties:\n');
fprintf('  Dimensions: %d\n', mesh.NDims);
fprintf('  Vertices: %d\n', mesh.NVertices);
fprintf('  Elements: %d\n', mesh.NElements);
fprintf('  Faces: %d\n', mesh.NFaces);
fprintf('  Boundary faces: %d\n', mesh.NBoundaryFaces);

%% Set up periodic boundaries
fprintf('\nStep 2: Setting up periodic boundaries...\n');

% Define bounding box for unit square [xMin, xMax, yMin, yMax]
boundingBox = [0, 1, 0, 1];
tolerance = 1e-10;

% Set periodic boundaries using bounding box
mesh = mesh.setPeriodic(boundingBox, tolerance);

fprintf('Periodic mapping established with tolerance = %.2e\n', tolerance);

% Count total vertex-to-vertex connections
totalV2VConnections = 0;
if ~isempty(mesh.VertexToVertexTables)
    for iTable = 1:length(mesh.VertexToVertexTables)
        totalV2VConnections = totalV2VConnections + nnz(mesh.VertexToVertexTables{iTable});
    end
end

fprintf('Vertex-to-vertex connections: %d\n', totalV2VConnections);
fprintf('Face-to-face connections: %d\n', nnz(mesh.FaceToFaceTable));

%% Validate periodic connections
fprintf('\nStep 3: Validating periodic connections...\n');

validatePeriodicConnections(mesh, boundingBox, tolerance);

%% Test find functionality
fprintf('\nStep 4: Testing find functionality...\n');

testFindFunctionality(mesh);

%% Validate Jacobian computations
fprintf('\nStep 5: Validating Jacobian computations...\n');

validateJacobianComputations(mesh);

%% Create visualizations
fprintf('\nStep 6: Creating visualizations...\n');

createMeshStructureFigure(mesh);
createValidationResultsFigure(mesh);

fprintf('\nMesh validation completed successfully!\n');
fprintf('=====================================\n');

%% Helper Functions

function [vertices, elements] = generateUniformBoxMesh(nx, ny, xRange, yRange)
% GENERATEUNIFORMBOXMESH Generate uniform triangular mesh on rectangular box.
%
%   [vertices, elements] = generateUniformBoxMesh(nx, ny, xRange, yRange)
%   creates a uniform triangular mesh on a rectangular domain.

% Generate grid points
xPoints = linspace(xRange(1), xRange(2), nx+1);
yPoints = linspace(yRange(1), yRange(2), ny+1);
[X, Y] = meshgrid(xPoints, yPoints);

% Convert to vertex list
vertices = [X(:), Y(:)];
nVertices = size(vertices, 1);

% Generate triangular elements
% Each grid cell (i,j) creates two triangles
elements = [];
for j = 1:ny
    for i = 1:nx
        % Grid indices for current cell
        v1 = (j - 1) * (nx + 1) + i; % Bottom-left
        v2 = (j - 1) * (nx + 1) + i + 1; % Bottom-right
        v3 = j * (nx + 1) + i; % Top-left
        v4 = j * (nx + 1) + i + 1; % Top-right

        % Two triangles per cell
        elements = [elements; v1, v2, v3]; % Lower triangle
        elements = [elements; v2, v4, v3]; % Upper triangle
    end
end

fprintf('Generated %d vertices and %d triangular elements.\n', ...
    nVertices, size(elements, 1));

end

function validatePeriodicConnections(mesh, boundingBox, tolerance)
% VALIDATEPERIODICCONNECTIONS Validate periodic vertex and face mappings.
%
%   validatePeriodicConnections(mesh, boundingBox, tolerance) checks that
%   periodic connections are correctly established.

fprintf('Validating periodic vertex connections...\n');

vertices = mesh.Vertices;
nConnections = 0;

% Define expected periodic shifts based on bounding box
xRange = [boundingBox(1), boundingBox(2)];
yRange = [boundingBox(3), boundingBox(4)];
periodicShifts = [; ...
    diff(xRange), 0; ... % Right boundary to left
    -diff(xRange), 0; ... % Left boundary to right
    0, diff(yRange); ... % Top boundary to bottom
    0, -diff(yRange); ... % Bottom boundary to top
    diff(xRange), diff(yRange); ... % Top-right to bottom-left
    -diff(xRange), -diff(yRange); ... % Bottom-left to top-right
    diff(xRange), -diff(yRange); ... % Bottom-right to top-left
    -diff(xRange), diff(yRange); ... % Top-left to bottom-right
    ];

% Test vertex-to-vertex mapping
if ~isempty(mesh.VertexToVertexTables)
    for iTable = 1:length(mesh.VertexToVertexTables)
        V2V = mesh.VertexToVertexTables{iTable};
        [iVertices, jVertices] = find(V2V);

        for k = 1:length(iVertices)
            i = iVertices(k);
            j = jVertices(k);

            if i ~= j
                % Check if connection corresponds to a periodic shift
                displacement = vertices(j, :) - vertices(i, :);

                % Find matching shift
                shiftFound = false;
                for s = 1:size(periodicShifts, 1)
                    if norm(displacement-periodicShifts(s, :)) < tolerance
                        shiftFound = true;
                        nConnections = nConnections + 1;
                        break;
                    end
                end

                if ~shiftFound && norm(displacement) > tolerance
                    fprintf('Warning: Unexpected connection between vertices %d and %d\n', i, j);
                    fprintf('  Displacement: [%.6f, %.6f]\n', displacement);
                end
            end
        end
    end
end

fprintf('Validated %d periodic vertex connections.\n', nConnections);

% Test face-to-face mapping
if ~isempty(mesh.FaceToFaceTable)
    F2F = mesh.FaceToFaceTable;
    faces = mesh.Faces;
    nFaceConnections = 0;

    fprintf('Validating periodic face connections...\n');

    [iFaces, jFaces] = find(F2F);

    for k = 1:length(iFaces)
        i = iFaces(k);
        j = jFaces(k);

        if i ~= j
            % Check if face connection is consistent with vertex connections
            face1Vertices = faces(i, :);
            face2Vertices = faces(j, :);

            % Check if all vertices in face1 have corresponding vertices in face2
            mappingConsistent = true;
            if ~isempty(mesh.VertexToVertexTables)
                for v1 = face1Vertices
                    mapped = false;
                    for iTable = 1:length(mesh.VertexToVertexTables)
                        V2V = mesh.VertexToVertexTables{iTable};
                        for v2 = face2Vertices
                            if V2V(v1, v2)
                                mapped = true;
                                break;
                            end
                        end
                        if mapped
                            break;
                        end
                    end
                    if ~mapped
                        mappingConsistent = false;
                        break;
                    end
                end
            end

            if mappingConsistent
                nFaceConnections = nFaceConnections + 1;
            else
                fprintf('Warning: Inconsistent face mapping between faces %d and %d\n', i, j);
            end
        end
    end

    fprintf('Validated %d periodic face connections.\n', nFaceConnections);
else
    fprintf('No face-to-face mapping table found.\n');
end

end

function testFindFunctionality(mesh)
% TESTFINDFUNCTIONALITY Test element finding methods.
%
%   testFindFunctionality(mesh) tests various element finding methods
%   including internal elements, boundary elements, and neighbors.

fprintf('Testing getInternalElements...\n');
interiorElements = mesh.getInternalElements();
fprintf('Found %d internal elements (out of %d total).\n', ...
    length(interiorElements), mesh.NElements);

fprintf('Testing getBoundaryElements...\n');
nFacesPerElement = mesh.NVerticesPerElement;
totalBoundaryElements = [];

for iFace = 1:nFacesPerElement
    boundaryElements = mesh.getBoundaryElements(iFace);
    totalBoundaryElements = [totalBoundaryElements; boundaryElements(:)];
    fprintf('  Face %d: %d boundary elements\n', iFace, length(boundaryElements));
end

uniqueBoundaryElements = unique(totalBoundaryElements);
fprintf('Total unique boundary elements: %d\n', length(uniqueBoundaryElements));

% Verify interior + boundary = total
allElements = union(interiorElements, uniqueBoundaryElements);
if length(allElements) == mesh.NElements
    fprintf('✓ Internal and boundary elements cover all elements.\n');
else
    fprintf('✗ Element coverage mismatch!\n');
end

fprintf('Testing findNeighborElements...\n');

% Test with a sample of elements
testElements = [1, ceil(mesh.NElements/2), mesh.NElements];
testElements = testElements(testElements <= mesh.NElements);

for iElement = testElements
    fprintf('  Element %d neighbors:\n', iElement);

    for iFace = 1:nFacesPerElement
        try
            % Test dirichlet boundary conditions
            [neighborsDirichlet, facesDirichlet] = mesh.findNeighborElements(iElement, iFace, 'dirichlet');

            % Test periodic boundary conditions
            [neighborsPeriodic, facesPeriodic] = mesh.findNeighborElements(iElement, iFace, 'periodic');

            fprintf('    Face %d: %d dirichlet, %d periodic\n', ...
                iFace, length(neighborsDirichlet), length(neighborsPeriodic));
        catch ME
            fprintf('    Face %d: Error - %s\n', iFace, ME.message);
        end
    end
end

end

function validateJacobianComputations(mesh)
% VALIDATEJACOBIANCOMPUTATIONS Validate Jacobian computation methods.
%
%   validateJacobianComputations(mesh) tests and validates various
%   Jacobian-related computations for geometric consistency.

fprintf('Computing element Jacobian determinants...\n');
detJ = mesh.computeElementJacobianDeterminants();
fprintf('Jacobian determinants: min=%.6f, max=%.6f, mean=%.6f\n', ...
    min(detJ), max(detJ), mean(detJ));

% Check for positive determinants (proper orientation)
negativeDetJ = sum(detJ < 0);
if negativeDetJ == 0
    fprintf('✓ All elements have positive orientation.\n');
else
    fprintf('✗ Warning: %d elements have negative orientation!\n', negativeDetJ);
end

fprintf('Computing element inverse Jacobians...\n');
invJ = mesh.computeElementInverseJacobians();
fprintf('Inverse Jacobians computed for %d elements.\n', size(invJ, 3));

% Validate J * inv(J) = I for a few elements
fprintf('Validating Jacobian inversion...\n');
J = mesh.computeElementJacobians();
testElements = [1, ceil(mesh.NElements/2), mesh.NElements];
maxError = 0;

for iElement = testElements
    if iElement <= mesh.NElements
        if mesh.NDims == 1
            Ji = J(iElement);
            invJi = invJ(iElement);
            product = Ji * invJi;
            error = abs(product-1);
        else
            Ji = J(:, :, iElement);
            invJi = invJ(:, :, iElement);
            product = Ji * invJi;
            identity = eye(size(product));
            error = norm(product-identity, 'fro');
        end
        maxError = max(maxError, error);
    end
end

fprintf('Maximum Jacobian inversion error: %.2e\n', maxError);
if maxError < 1e-12
    fprintf('✓ Jacobian inversion is accurate.\n');
else
    fprintf('✗ Warning: Jacobian inversion has significant errors!\n');
end

fprintf('Computing face Jacobian determinants...\n');
nFacesPerElement = mesh.NVerticesPerElement;
for iFace = 1:nFacesPerElement
    detJFace = mesh.computeFaceJacobianDeterminants([], iFace);
    fprintf('Face %d: min=%.6f, max=%.6f, mean=%.6f\n', ...
        iFace, min(detJFace), max(detJFace), mean(detJFace));
end

fprintf('Computing outward normal vectors...\n');
for iFace = 1:nFacesPerElement
    normals = mesh.computeOutwardNormals([], iFace);

    % Check unit length
    normLengths = sqrt(sum(normals.^2, 1));
    maxNormError = max(abs(normLengths-1));

    fprintf('Face %d normals: max deviation from unit length = %.2e\n', ...
        iFace, maxNormError);
end

fprintf('Computing mesh measures...\n');
detJ = mesh.computeElementJacobianDeterminants();
elementMeasures = sum(detJ) / 2; % Divide by 2 for triangular elements
fprintf('Total mesh measure: %.6f\n', elementMeasures);

% For unit square, total area should be 1
expectedArea = 1.0;
areaError = abs(elementMeasures-expectedArea);
fprintf('Area error compared to unit square: %.6f\n', areaError);

if areaError < 1e-10
    fprintf('✓ Mesh area computation is accurate.\n');
else
    fprintf('✗ Warning: Mesh area computation has errors!\n');
end

end

function createMeshStructureFigure(mesh)
% CREATEMESHSTRUCTUREFIGURE Create figure showing mesh structure and connections.
%
%   createMeshStructureFigure(mesh) creates a comprehensive figure showing
%   mesh structure, periodic connections, and element classification.

figure('Name', 'Mesh Structure and Periodic Connections', ...
    'Position', [100, 100, 1400, 500]);

%% Subplot 1: Basic mesh structure with boundary highlighting
subplot(1, 3, 1);
vertices = mesh.Vertices;
elements = mesh.Elements;

% Plot elements
for iElement = 1:mesh.NElements
    elementVertices = vertices(elements(iElement, :), :);
    fill(elementVertices(:, 1), elementVertices(:, 2), [0.9, 0.9, 0.9], ...
        'EdgeColor', [0.5, 0.5, 0.5], 'LineWidth', 1);
    hold on;
end

% Highlight boundary faces
faces = mesh.Faces;
boundary = mesh.Boundary;
for iBoundaryFace = boundary'
    faceVertices = vertices(faces(iBoundaryFace, :), :);
    plot(faceVertices(:, 1), faceVertices(:, 2), 'r-', 'LineWidth', 3);
end

% Add vertices
plot(vertices(:, 1), vertices(:, 2), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'black');

% Add vertex labels
for iVertex = 1:mesh.NVertices
    text(vertices(iVertex, 1)+0.02, vertices(iVertex, 2)+0.02, ...
        sprintf('%d', iVertex), 'FontSize', 8, 'FontWeight', 'bold');
end

axis equal;
grid on;
xlim([-0.1, 1.1]);
ylim([-0.1, 1.1]);
xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
title('Mesh Structure with Boundary', 'FontSize', 14, 'FontWeight', 'bold');

%% Subplot 2: Periodic connections
subplot(1, 3, 2);

% Plot mesh structure lightly
for iElement = 1:mesh.NElements
    elementVertices = vertices(elements(iElement, :), :);
    fill(elementVertices(:, 1), elementVertices(:, 2), [0.95, 0.95, 0.95], ...
        'EdgeColor', [0.7, 0.7, 0.7], 'LineWidth', 0.5);
    hold on;
end

% Plot vertices
plot(vertices(:, 1), vertices(:, 2), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'black');

% Add periodic connections
if ~isempty(mesh.VertexToVertexTables)
    connectionColors = lines(length(mesh.VertexToVertexTables));

    for iTable = 1:length(mesh.VertexToVertexTables)
        V2V = mesh.VertexToVertexTables{iTable};
        [i, j] = find(V2V);

        for k = 1:length(i)
            if i(k) ~= j(k) % Avoid self-connections
                v1 = vertices(i(k), :);
                v2 = vertices(j(k), :);

                plot([v1(1), v2(1)], [v1(2), v2(2)], '--', ...
                    'Color', connectionColors(iTable, :), 'LineWidth', 2);
            end
        end
    end
end

axis equal;
grid on;
xlim([-0.1, 1.1]);
ylim([-0.1, 1.1]);
xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
title('Periodic Connections', 'FontSize', 14, 'FontWeight', 'bold');

%% Subplot 3: Element classification
subplot(1, 3, 3);

% Find element classifications
interiorElements = mesh.getInternalElements();
allElements = 1:mesh.NElements;
boundaryElements = setdiff(allElements, interiorElements);

% Plot internal elements in blue
interiorHandle = [];
if ~isempty(interiorElements)
    iElement = interiorElements(1);
    elementVertices = vertices(elements(iElement, :), :);
    interiorHandle = fill(elementVertices(:, 1), elementVertices(:, 2), [0.2, 0.6, 1.0], ...
        'EdgeColor', [0.1, 0.3, 0.5], 'LineWidth', 1, 'FaceAlpha', 0.8);
    hold on;

    % Add element number
    centroid = mean(elementVertices, 1);
    text(centroid(1), centroid(2), sprintf('%d', iElement), ...
        'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'Color', 'white');

    % Plot remaining internal elements
    for iElement = interiorElements(2:end)
        elementVertices = vertices(elements(iElement, :), :);
        fill(elementVertices(:, 1), elementVertices(:, 2), [0.2, 0.6, 1.0], ...
            'EdgeColor', [0.1, 0.3, 0.5], 'LineWidth', 1, 'FaceAlpha', 0.8);

        % Add element number
        centroid = mean(elementVertices, 1);
        text(centroid(1), centroid(2), sprintf('%d', iElement), ...
            'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'Color', 'white');
    end
end

% Plot boundary elements in orange
boundaryHandle = [];
if ~isempty(boundaryElements)
    iElement = boundaryElements(1);
    elementVertices = vertices(elements(iElement, :), :);
    boundaryHandle = fill(elementVertices(:, 1), elementVertices(:, 2), [1.0, 0.6, 0.2], ...
        'EdgeColor', [0.5, 0.3, 0.1], 'LineWidth', 1, 'FaceAlpha', 0.8);

    % Add element number
    centroid = mean(elementVertices, 1);
    text(centroid(1), centroid(2), sprintf('%d', iElement), ...
        'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'Color', 'white');

    % Plot remaining boundary elements
    for iElement = boundaryElements(2:end)
        elementVertices = vertices(elements(iElement, :), :);
        fill(elementVertices(:, 1), elementVertices(:, 2), [1.0, 0.6, 0.2], ...
            'EdgeColor', [0.5, 0.3, 0.1], 'LineWidth', 1, 'FaceAlpha', 0.8);

        % Add element number
        centroid = mean(elementVertices, 1);
        text(centroid(1), centroid(2), sprintf('%d', iElement), ...
            'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'Color', 'white');
    end
end

axis equal;
grid on;
xlim([-0.1, 1.1]);
ylim([-0.1, 1.1]);
xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
title('Element Classification', 'FontSize', 14, 'FontWeight', 'bold');

% Create legend with proper handles
legendHandles = [];
legendLabels = {};
if ~isempty(interiorHandle)
    legendHandles = [legendHandles, interiorHandle];
    legendLabels = [legendLabels, {'Internal Elements'}];
end
if ~isempty(boundaryHandle)
    legendHandles = [legendHandles, boundaryHandle];
    legendLabels = [legendLabels, {'Boundary Elements'}];
end
if ~isempty(legendHandles)
    legend(legendHandles, legendLabels, 'Location', 'northeastoutside');
end

end

function createValidationResultsFigure(mesh)
% CREATEVALIDATIONRESULTSFIGURE Create figure showing validation results.
%
%   createValidationResultsFigure(mesh) creates a figure showing Jacobian
%   determinants, normal vectors, and neighbor connectivity.

figure('Name', 'Validation Results: Jacobians and Connectivity', ...
    'Position', [150, 150, 1400, 500]);

vertices = mesh.Vertices;
elements = mesh.Elements;

%% Subplot 1: Jacobian determinants with colormap
subplot(1, 3, 1);
detJ = mesh.computeElementJacobianDeterminants();

% Plot elements colored by Jacobian determinant
for iElement = 1:mesh.NElements
    elementVertices = vertices(elements(iElement, :), :);
    fill(elementVertices(:, 1), elementVertices(:, 2), detJ(iElement), ...
        'EdgeColor', 'black', 'LineWidth', 0.8);
    hold on;

    % Add Jacobian value text
    centroid = mean(elementVertices, 1);
    text(centroid(1), centroid(2), sprintf('%.3f', detJ(iElement)), ...
        'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'Color', 'white', 'BackgroundColor', 'black', 'Margin', 1);
end

axis equal;
grid on;
xlim([-0.05, 1.05]);
ylim([-0.05, 1.05]);
xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
title('Jacobian Determinants', 'FontSize', 14, 'FontWeight', 'bold');
colorbar('FontSize', 10);

%% Subplot 2: Normal vectors on boundary faces
subplot(1, 3, 2);

% Plot mesh structure
for iElement = 1:mesh.NElements
    elementVertices = vertices(elements(iElement, :), :);
    fill(elementVertices(:, 1), elementVertices(:, 2), [0.9, 0.9, 0.9], ...
        'EdgeColor', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
    hold on;
end

% Plot normal vectors for each face type
scaleFactor = 0.1;
nFacesPerElement = mesh.NVerticesPerElement;
normalColors = {'red', 'green', 'blue'};

for iFace = 1:nFacesPerElement
    try
        normals = mesh.computeOutwardNormals([], iFace);
        boundaryElements = mesh.getBoundaryElements(iFace);

        for iElement = boundaryElements'
            if iElement <= mesh.NElements
                elementVertices = vertices(elements(iElement, :), :);

                % Compute face center
                faceVertexIndices = setdiff(1:nFacesPerElement, iFace);
                faceVertices = elementVertices(faceVertexIndices, :);
                faceCenter = mean(faceVertices, 1);

                % Plot normal vector
                if size(normals, 2) >= iElement
                    normal = normals(:, iElement);
                    quiver(faceCenter(1), faceCenter(2), ...
                        normal(1)*scaleFactor, normal(2)*scaleFactor, ...
                        'Color', normalColors{iFace}, 'LineWidth', 3, 'MaxHeadSize', 0.4);
                end
            end
        end
    catch ME
        fprintf('Warning: Could not compute normals for face %d: %s\n', iFace, ME.message);
    end
end

axis equal;
grid on;
xlim([-0.05, 1.05]);
ylim([-0.05, 1.05]);
xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
title('Outward Normal Vectors', 'FontSize', 14, 'FontWeight', 'bold');

%% Subplot 3: Neighbor connectivity demonstration
subplot(1, 3, 3);

% Plot mesh structure lightly
for iElement = 1:mesh.NElements
    elementVertices = vertices(elements(iElement, :), :);
    fill(elementVertices(:, 1), elementVertices(:, 2), [0.95, 0.95, 0.95], ...
        'EdgeColor', [0.7, 0.7, 0.7], 'LineWidth', 0.5);
    hold on;
end

% Select sample elements for neighbor demonstration
sampleElements = [5, 16, 25]; % Strategic selection for good visualization
highlightColors = {[1, 0.3, 0.3], [0.3, 0.3, 1], [0.3, 1, 0.3]};
nFacesPerElement = mesh.NVerticesPerElement;

for iSample = 1:length(sampleElements)
    iElement = sampleElements(iSample);
    if iElement <= mesh.NElements
        color = highlightColors{iSample};

        % Highlight the sample element
        elementVertices = vertices(elements(iElement, :), :);
        fill(elementVertices(:, 1), elementVertices(:, 2), color, ...
            'FaceAlpha', 0.9, 'EdgeColor', 'black', 'LineWidth', 2);

        % Add element number
        centroid = mean(elementVertices, 1);
        text(centroid(1), centroid(2), sprintf('%d', iElement), ...
            'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'Color', 'white');

        % Find and highlight neighbors with periodic boundaries
        allNeighbors = [];
        for iFace = 1:nFacesPerElement
            [neighbors, ~] = mesh.findNeighborElements(iElement, iFace, 'periodic');
            allNeighbors = [allNeighbors, neighbors];

        end

        % Remove duplicates and highlight neighbors
        uniqueNeighbors = unique(allNeighbors);
        for neighbor = uniqueNeighbors
            if neighbor > 0 && neighbor <= mesh.NElements
                neighborVertices = vertices(elements(neighbor, :), :);
                plot(neighborVertices([1:end, 1], 1), neighborVertices([1:end, 1], 2), ...
                    'Color', color, 'LineWidth', 3);
            end
        end
    end
end

axis equal;
grid on;
xlim([-0.05, 1.05]);
ylim([-0.05, 1.05]);
xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
title('Neighbor Connectivity (Periodic)', 'FontSize', 14, 'FontWeight', 'bold');

end