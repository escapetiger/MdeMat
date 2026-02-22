clear; close all; clc;

fprintf('Starting Mesh Validation Example...\n');
fprintf('=======================================\n');

%% Generate uniform 4x4 SimplexMesh on unit box
fprintf('Step 1: Generating uniform 4x4 SimplexMesh on [0,1]x[0,1]...\n');

[vertices, elements] = generateUniformBoxMesh(4, 4, [0, 1], [0, 1]);

fprintf('Generated mesh with %d vertices and %d elements.\n', ...
    size(vertices, 1), size(elements, 1));

% Create SimplexMesh object
mesh = approx.mesh.SimplexMesh(vertices, elements);

fprintf('SimplexMesh properties:\n');
fprintf('  Dimensions: %d\n', mesh.nDims);
fprintf('  Vertices: %d\n', mesh.nVertices);
fprintf('  Elements: %d\n', mesh.nElements);
fprintf('  Faces: %d\n', mesh.nFaces);
fprintf('  Boundary faces: %d\n', mesh.nBoundaryFaces);

%% Set up periodic boundaries
fprintf('\nStep 2: Setting up periodic boundaries...\n');

% Define periodic shifts for unit box
% Left-right: shift by [1, 0]
% Top-bottom: shift by [0, 1]
% Corner connections: shift by [1, 1], [1, -1], [-1, 1], [-1, -1]
periodicShifts = [
    1.0,  0.0;   % Right boundary to left
   -1.0,  0.0;   % Left boundary to right
    0.0,  1.0;   % Top boundary to bottom
    0.0, -1.0;   % Bottom boundary to top
    1.0,  1.0;   % Top-right to bottom-left
   -1.0, -1.0;   % Bottom-left to top-right
    1.0, -1.0;   % Bottom-right to top-left
   -1.0,  1.0    % Top-left to bottom-right
];

tolerance = 1e-10;
mesh = mesh.setPeriodic(periodicShifts, tolerance);

fprintf('Periodic mapping established with tolerance = %.2e\n', tolerance);
fprintf('Vertex-to-vertex connections: %d\n', nnz(mesh.vertexToVertexTable));
fprintf('Face-to-face connections: %d\n', nnz(mesh.faceToFaceTable));

%% Validate periodic connections
fprintf('\nStep 3: Validating periodic connections...\n');

validatePeriodicConnections(mesh, periodicShifts, tolerance);

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
%
% Inputs:
%   nx - Number of divisions in x-direction
%   ny - Number of divisions in y-direction  
%   xRange - [xMin, xMax] range in x-direction
%   yRange - [yMin, yMax] range in y-direction
%
% Outputs:
%   vertices - Vertex coordinates (nVertices × 2)
%   elements - Element connectivity (nElements × 3)

% Generate grid points
xPoints = linspace(xRange(1), xRange(2), nx + 1);
yPoints = linspace(yRange(1), yRange(2), ny + 1);
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
        v1 = (j-1) * (nx+1) + i;       % Bottom-left
        v2 = (j-1) * (nx+1) + i + 1;   % Bottom-right  
        v3 = j * (nx+1) + i;           % Top-left
        v4 = j * (nx+1) + i + 1;       % Top-right
        
        % Two triangles per cell
        elements = [elements; v1, v2, v3];  % Lower triangle
        elements = [elements; v2, v4, v3];  % Upper triangle
    end
end

fprintf('Generated %d vertices and %d triangular elements.\n', ...
    nVertices, size(elements, 1));

end

function validatePeriodicConnections(mesh, shifts, tolerance)
% VALIDATEPERIODICCONNECTIONS Validate periodic vertex and face mappings.
%
%   validatePeriodicConnections(mesh, shifts, tolerance) checks that
%   periodic connections are correctly established.

fprintf('Validating periodic vertex connections...\n');

% Test vertex-to-vertex mapping
V2V = mesh.vertexToVertexTable;
vertices = mesh.vertices;
nConnections = 0;

for i = 1:mesh.nVertices
    % Find connected vertices
    connectedVertices = find(V2V(i, :));
    
    if ~isempty(connectedVertices)
        for j = connectedVertices
            if i ~= j
                % Check if connection corresponds to a periodic shift
                displacement = vertices(j, :) - vertices(i, :);
                
                % Find matching shift
                shiftFound = false;
                for k = 1:size(shifts, 1)
                    if norm(displacement - shifts(k, :)) < tolerance
                        shiftFound = true;
                        nConnections = nConnections + 1;
                        break;
                    end
                end
                
                if ~shiftFound
                    fprintf('Warning: Unexpected connection between vertices %d and %d\n', i, j);
                end
            end
        end
    end
end

fprintf('Validated %d periodic vertex connections.\n', nConnections);

% Test face-to-face mapping
F2F = mesh.faceToFaceTable;
faces = mesh.faces;
nFaceConnections = 0;

fprintf('Validating periodic face connections...\n');

for i = 1:mesh.nFaces
    connectedFaces = find(F2F(i, :));
    
    for j = connectedFaces
        if i ~= j
            % Check if face connection is consistent with vertex connections
            face1Vertices = faces(i, :);
            face2Vertices = faces(j, :);
            
            % All vertices in face1 should map to vertices in face2
            mappingConsistent = true;
            for v1 = face1Vertices
                mapped = false;
                for v2 = face2Vertices
                    if V2V(v1, v2)
                        mapped = true;
                        break;
                    end
                end
                if ~mapped
                    mappingConsistent = false;
                    break;
                end
            end
            
            if mappingConsistent
                nFaceConnections = nFaceConnections + 1;
            else
                fprintf('Warning: Inconsistent face mapping between faces %d and %d\n', i, j);
            end
        end
    end
end

fprintf('Validated %d periodic face connections.\n', nFaceConnections);

end

function testFindFunctionality(mesh)
% TESTFINDFUNCTIONALITY Test element finding methods.
%
%   testFindFunctionality(mesh) tests various element finding methods
%   including interior elements, boundary elements, and neighbors.

fprintf('Testing findInteriorElements...\n');
interiorElements = mesh.findInteriorElements();
fprintf('Found %d interior elements (out of %d total).\n', ...
    length(interiorElements), mesh.nElements);

fprintf('Testing findBoundaryElements...\n');
nFacesPerElement = mesh.nVerticesPerElement;
totalBoundaryElements = [];

for iFace = 1:nFacesPerElement
    boundaryElements = mesh.findBoundaryElements(iFace);
    totalBoundaryElements = [totalBoundaryElements; boundaryElements(:)];
    fprintf('  Face %d: %d boundary elements\n', iFace, length(boundaryElements));
end

uniqueBoundaryElements = unique(totalBoundaryElements);
fprintf('Total unique boundary elements: %d\n', length(uniqueBoundaryElements));

% Verify interior + boundary = total
allElements = union(interiorElements, uniqueBoundaryElements);
if length(allElements) == mesh.nElements
    fprintf('✓ Interior and boundary elements cover all elements.\n');
else
    fprintf('✗ Element coverage mismatch!\n');
end

fprintf('Testing findNeighborElements...\n');

% Test with a sample of elements
testElements = [1, ceil(mesh.nElements/2), mesh.nElements];
testElements = testElements(testElements <= mesh.nElements);

for iElement = testElements
    fprintf('  Element %d neighbors:\n', iElement);
    
    for iFace = 1:nFacesPerElement
        % Test strict boundary conditions
        neighborsStrict = mesh.findNeighborElements(iFace, iElement, 'strict');
        
        % Test periodic boundary conditions  
        neighborsPeriodic = mesh.findNeighborElements(iFace, iElement, 'periodic');
        
        fprintf('    Face %d: %d strict, %d periodic\n', ...
            iFace, length(neighborsStrict), length(neighborsPeriodic));
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
testElements = [1, ceil(mesh.nElements/2), mesh.nElements];
maxError = 0;

for iElement = testElements
    if iElement <= mesh.nElements
        Ji = J(:, :, iElement);
        invJi = invJ(:, :, iElement);
        product = Ji * invJi;
        identity = eye(size(product));
        error = norm(product - identity, 'fro');
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
nFacesPerElement = mesh.nVerticesPerElement;
for iFace = 1:nFacesPerElement
    detJFace = mesh.computeFaceJacobianDeterminants(iFace);
    fprintf('Face %d: min=%.6f, max=%.6f, mean=%.6f\n', ...
        iFace, min(detJFace), max(detJFace), mean(detJFace));
end

fprintf('Computing outward normal vectors...\n');
for iFace = 1:nFacesPerElement
    normals = mesh.computeOutwardNormals(iFace);
    
    % Check unit length
    normLengths = sqrt(sum(normals.^2, 1));
    maxNormError = max(abs(normLengths - 1));
    
    fprintf('Face %d normals: max deviation from unit length = %.2e\n', ...
        iFace, maxNormError);
end

fprintf('Computing absolute Jacobian determinants...\n');
detJ = mesh.computeElementJacobianDeterminants();
elementMeasures = sum(detJ) / 2;
fprintf('Total mesh measure: %.6f\n', elementMeasures);

% For unit square, total area should be 1
expectedArea = 1.0;
areaError = abs(elementMeasures - expectedArea);
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
%   mesh structure, periodic connections, and element classification using
%   MeshPlotter for professional visualization.

% Create mesh plotter
meshPlotter = approx.mesh.MeshPlotter();

figure('Name', 'Mesh Structure and Periodic Connections', ...
    'Position', [100, 100, 1400, 500]);

%% Subplot 1: Basic mesh structure with boundary highlighting
subplot(1, 3, 1);
meshPlotter.plot(mesh, 'ShowElements', true, 'ElementColor', [0.9, 0.9, 0.9], ...
    'ShowBoundary', true, 'BoundaryColor', 'red', 'ShowVertices', true, ...
    'VertexColor', 'black', 'VertexSize', 50, 'EdgeColor', [0.5, 0.5, 0.5], ...
    'Alpha', 0.8);

% Add vertex labels
vertices = mesh.vertices;
for iVertex = 1:mesh.nVertices
    text(vertices(iVertex, 1) + 0.02, vertices(iVertex, 2) + 0.02, ...
        sprintf('%d', iVertex), 'FontSize', 8, 'FontWeight', 'bold');
end

xlim([-0.1, 1.1]);
ylim([-0.1, 1.1]);
title('Mesh Structure with Boundary', 'FontSize', 14, 'FontWeight', 'bold');

%% Subplot 2: Periodic connections
subplot(1, 3, 2);
meshPlotter.plot(mesh, 'ShowElements', true, 'ElementColor', [0.95, 0.95, 0.95], ...
    'ShowVertices', true, 'VertexColor', 'black', 'VertexSize', 40, ...
    'ShowEdges', false, 'ShowBoundary', false, 'Alpha', 0.3);

hold on;

% Add periodic connections with different colors for different types
V2V = mesh.vertexToVertexTable;
[i, j] = find(V2V);
connectionColors = lines(4); % Different colors for different connection types

for k = 1:length(i)
    if i(k) ~= j(k)  % Avoid self-connections
        v1 = vertices(i(k), :);
        v2 = vertices(j(k), :);
        displacement = v2 - v1;
        
        % Determine connection type based on displacement
        if abs(displacement(1)) > 0.5 && abs(displacement(2)) < 0.1
            colorIdx = 1; % Horizontal
        elseif abs(displacement(2)) > 0.5 && abs(displacement(1)) < 0.1
            colorIdx = 2; % Vertical
        elseif abs(displacement(1)) > 0.5 && abs(displacement(2)) > 0.5
            colorIdx = 3; % Diagonal
        else
            colorIdx = 4; % Other
        end
        
        plot([v1(1), v2(1)], [v1(2), v2(2)], '--', ...
            'Color', connectionColors(colorIdx, :), 'LineWidth', 2);
    end
end

xlim([-0.1, 1.1]);
ylim([-0.1, 1.1]);
title('Periodic Connections', 'FontSize', 14, 'FontWeight', 'bold');
legend('Elements', 'Vertices', 'Horizontal', 'Vertical', 'Diagonal', 'Other', ...
    'Location', 'northeastoutside');

%% Subplot 3: Element classification
subplot(1, 3, 3);

% Find element classifications
interiorElements = mesh.findInteriorElements();
allElements = 1:mesh.nElements;
boundaryElements = setdiff(allElements, interiorElements);

% Create custom element visualization
elements = mesh.elements;

% Plot interior elements in blue
for iElement = interiorElements
    elementVertices = vertices(elements(iElement, :), :);
    fill(elementVertices(:, 1), elementVertices(:, 2), [0.2, 0.6, 1.0], ...
        'EdgeColor', [0.1, 0.3, 0.5], 'LineWidth', 1, 'FaceAlpha', 0.8);
    hold on;
    
    % Add element number
    centroid = mean(elementVertices, 1);
    text(centroid(1), centroid(2), sprintf('%d', iElement), ...
        'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'Color', 'white');
end

% Plot boundary elements in orange
for iElement = boundaryElements
    elementVertices = vertices(elements(iElement, :), :);
    fill(elementVertices(:, 1), elementVertices(:, 2), [1.0, 0.6, 0.2], ...
        'EdgeColor', [0.5, 0.3, 0.1], 'LineWidth', 1, 'FaceAlpha', 0.8);
    
    % Add element number
    centroid = mean(elementVertices, 1);
    text(centroid(1), centroid(2), sprintf('%d', iElement), ...
        'FontSize', 10, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
        'Color', 'white');
end

axis equal;
grid on;
xlim([-0.1, 1.1]);
ylim([-0.1, 1.1]);
xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
title('Element Classification', 'FontSize', 14, 'FontWeight', 'bold');
legend('Interior Elements', 'Location', 'northeastoutside');

end

function createValidationResultsFigure(mesh)
% CREATEVALIDATIONRESULTSFIGURE Create figure showing validation results.
%
%   createValidationResultsFigure(mesh) creates a figure showing Jacobian
%   determinants, normal vectors, and neighbor connectivity.

figure('Name', 'Validation Results: Jacobians and Connectivity', ...
    'Position', [150, 150, 1400, 500]);

%% Subplot 1: Jacobian determinants with colormap
subplot(1, 3, 1);
vertices = mesh.vertices;
elements = mesh.elements;
detJ = mesh.computeElementJacobianDeterminants();

% Create smooth colormap
colormap(parula);

% Plot elements colored by Jacobian determinant
for iElement = 1:mesh.nElements
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
faces = mesh.faces;
boundary = mesh.boundary;

% Plot mesh structure
for iElement = 1:mesh.nElements
    elementVertices = vertices(elements(iElement, :), :);
    fill(elementVertices(:, 1), elementVertices(:, 2), [0.9, 0.9, 0.9], ...
        'EdgeColor', [0.5, 0.5, 0.5], 'LineWidth', 0.5);
    hold on;
end

% Plot normal vectors for each face type
scaleFactor = 0.1;
nFacesPerElement = mesh.nVerticesPerElement;
normalColors = {'red', 'green', 'blue'};

for iFace = 1:nFacesPerElement
    normals = mesh.computeOutwardNormals(iFace);
    boundaryElements = mesh.findBoundaryElements(iFace);
    
    for iElement = boundaryElements'
        if iElement <= mesh.nElements
            elementVertices = vertices(elements(iElement, :), :);
            
            % Compute face center
            faceVertexIndices = setdiff(1:nFacesPerElement, iFace);
            faceVertices = elementVertices(faceVertexIndices, :);
            faceCenter = mean(faceVertices, 1);
            
            % Plot normal vector
            normal = normals(:, iElement);
            quiver(faceCenter(1), faceCenter(2), ...
                normal(1)*scaleFactor, normal(2)*scaleFactor, ...
                'Color', normalColors{iFace}, 'LineWidth', 3, 'MaxHeadSize', 0.4);
        end
    end
end

axis equal;
grid on;
xlim([-0.05, 1.05]);
ylim([-0.05, 1.05]);
xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
title('Outward Normal Vectors', 'FontSize', 14, 'FontWeight', 'bold');
legend('Elements', 'Face 1 Normals', 'Face 2 Normals', 'Face 3 Normals', ...
    'Location', 'northeastoutside');

%% Subplot 3: Neighbor connectivity demonstration
subplot(1, 3, 3);

% Plot mesh structure lightly
for iElement = 1:mesh.nElements
    elementVertices = vertices(elements(iElement, :), :);
    fill(elementVertices(:, 1), elementVertices(:, 2), [0.95, 0.95, 0.95], ...
        'EdgeColor', [0.7, 0.7, 0.7], 'LineWidth', 0.5);
    hold on;
end

% Select sample elements for neighbor demonstration
sampleElements = [5, 16, 25]; % Strategic selection for good visualization
highlightColors = {[1, 0.3, 0.3], [0.3, 0.3, 1], [0.3, 1, 0.3]};
nFacesPerElement = mesh.nVerticesPerElement;

for iSample = 1:length(sampleElements)
    iElement = sampleElements(iSample);
    if iElement <= mesh.nElements
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
            neighbors = mesh.findNeighborElements(iFace, iElement, 'periodic');
            allNeighbors = [allNeighbors, neighbors];
        end
        
        % Remove duplicates and highlight neighbors
        uniqueNeighbors = unique(allNeighbors);
        for neighbor = uniqueNeighbors
            if neighbor <= mesh.nElements
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
legend('Mesh Elements', 'Element 5 & Neighbors', 'Element 16 & Neighbors', ...
    'Element 25 & Neighbors', 'Location', 'northeastoutside');

end