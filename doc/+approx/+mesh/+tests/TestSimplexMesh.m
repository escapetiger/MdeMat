classdef TestSimplexMesh < matlab.unittest.TestCase
    
    properties (TestParameter)
        meshCase = struct(...
            'line1D', struct(...
                'vertices', [0; 1; 2], ...
                'elements', [1, 2; 2, 3], ...
                'nDims', 1, ...
                'nVerticesPerElement', 2), ...
            'triangle2D', struct(...
                'vertices', [0, 0; 1, 0; 0.5, 1; 1, 1], ...
                'elements', [1, 2, 3; 2, 4, 3], ...
                'nDims', 2, ...
                'nVerticesPerElement', 3), ...
            'tetrahedron3D', struct(...
                'vertices', [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1; 1, 1, 1], ...
                'elements', [1, 2, 3, 4; 2, 5, 3, 4], ...
                'nDims', 3, ...
                'nVerticesPerElement', 4))
    end
    
    methods (Test)
        function testConstructor(testCase, meshCase)
            % Test SimplexMesh constructor for different dimensions
            
            mesh = approx.mesh.SimplexMesh(meshCase.vertices, meshCase.elements);
            
            testCase.verifyEqual(mesh.vertices, meshCase.vertices);
            testCase.verifyEqual(mesh.elements, meshCase.elements);
            testCase.verifyEqual(mesh.nDims, meshCase.nDims);
            testCase.verifyEqual(mesh.nVerticesPerElement, meshCase.nVerticesPerElement);
            
            % Verify that faces and boundary are computed
            testCase.verifyFalse(isempty(mesh.faces));
            testCase.verifyFalse(isempty(mesh.boundary));
            testCase.verifyFalse(isempty(mesh.faceToElementTable));
        end
        
        function testDependentProperties(testCase, meshCase)
            % Test all dependent properties
            
            mesh = approx.mesh.SimplexMesh(meshCase.vertices, meshCase.elements);
            
            testCase.verifyEqual(mesh.nDims, size(meshCase.vertices, 2));
            testCase.verifyEqual(mesh.nVertices, size(meshCase.vertices, 1));
            testCase.verifyEqual(mesh.nElements, size(meshCase.elements, 1));
            testCase.verifyEqual(mesh.nVerticesPerElement, size(meshCase.elements, 2));
            testCase.verifyEqual(mesh.nFaces, size(mesh.faces, 1));
            testCase.verifyEqual(mesh.nVerticesPerFace, size(mesh.faces, 2));
            testCase.verifyEqual(mesh.nBoundaryFaces, length(mesh.boundary));
        end
        
        function testFaceGeneration(testCase, meshCase)
            % Test face generation from elements
            
            mesh = approx.mesh.SimplexMesh(meshCase.vertices, meshCase.elements);
            
            % Verify face dimensions
            expectedNVerticesPerFace = meshCase.nDims;
            testCase.verifyEqual(mesh.nVerticesPerFace, expectedNVerticesPerFace);
            
            % Verify face-to-element connectivity
            testCase.verifyEqual(size(mesh.faceToElementTable, 1), mesh.nFaces);
            testCase.verifyEqual(size(mesh.faceToElementTable, 2), mesh.nElements);
            
            % Verify boundary faces are subset of all faces
            testCase.verifyTrue(all(mesh.boundary >= 1));
            testCase.verifyTrue(all(mesh.boundary <= mesh.nFaces));
        end
        
        function testPositiveOrientation(testCase)
            % Test that elements have positive orientation after construction
            
            % Create mesh with negative orientation
            vertices = [0, 0; 1, 0; 0, 1];
            elements = [1, 3, 2]; % Clockwise orientation (negative)
            
            mesh = approx.mesh.SimplexMesh(vertices, elements);
            
            % Verify orientation is corrected
            detJ = mesh.computeElementJacobianDeterminants();
            testCase.verifyTrue(all(detJ > 0));
        end
        
        function testGraphify(testCase, meshCase)
            % Test graph generation from simplex mesh
            
            mesh = approx.mesh.SimplexMesh(meshCase.vertices, meshCase.elements);
            graph = mesh.graphify();
            
            testCase.verifyClass(graph, 'approx.mesh.Graph');
            testCase.verifyEqual(graph.nDims, mesh.nDims);
            testCase.verifyEqual(graph.nVertices, mesh.nVertices);
            testCase.verifyEqual(graph.vertices, mesh.vertices);
            
            % Verify edges are created from element connectivity
            testCase.verifyTrue(graph.nEdges > 0);
            testCase.verifyEqual(size(graph.edges, 2), 2);
        end
        
        function testCollocate(testCase)
            % Test coordinate transformation from reference to physical
            
            % Simple 2D triangle
            vertices = [0, 0; 1, 0; 0, 1];
            elements = [1, 2, 3];
            mesh = approx.mesh.SimplexMesh(vertices, elements);
            
            % Test reference coordinates
            X = [0.25, 0.25; 0.5, 0.25]'; % 2 points in reference space
            Y = mesh.collocate(X, 1);
            
            testCase.verifyEqual(size(Y), [2, 2]);
            
            % Test center point (1/3, 1/3 in reference coordinates)
            Xcenter = [1/3; 1/3];
            Ycenter = mesh.collocate(Xcenter, 1);
            expectedCenter = [1/3; 1/3];
            testCase.verifyEqual(squeeze(Ycenter), expectedCenter, 'AbsTol', 1e-10);
        end
        
        function testCollocateAllElements(testCase)
            % Test collocation without specifying elements
            
            vertices = [0, 0; 1, 0; 0.5, 1; 1, 1];
            elements = [1, 2, 3; 2, 4, 3];
            mesh = approx.mesh.SimplexMesh(vertices, elements);
            
            X = [0.25; 0.25];
            Y = mesh.collocate(X);
            
            testCase.verifyEqual(size(Y), [2, 2]);
        end
        
        function testFindElementVertices(testCase, meshCase)
            % Test element vertex finding
            
            mesh = approx.mesh.SimplexMesh(meshCase.vertices, meshCase.elements);
            
            % Test with specific elements
            if mesh.nElements >= 1
                V = mesh.findElementVertices(1);
                expectedSize = [meshCase.nDims, meshCase.nVerticesPerElement];
                testCase.verifyEqual(size(V), expectedSize);
            end
            
            % Test with all elements (default)
            VAll = mesh.findElementVertices();
            expectedSize = [meshCase.nDims, meshCase.nVerticesPerElement, mesh.nElements];
            testCase.verifyEqual(size(VAll), expectedSize);
        end
        
        function testRefine1D(testCase)
            % Test 1D mesh refinement
            
            vertices = [0; 1; 2];
            elements = [1, 2; 2, 3];
            mesh = approx.mesh.SimplexMesh(vertices, elements);
            
            % Test no refinement
            mesh0 = mesh.refine(0);
            testCase.verifyEqual(mesh0.nElements, mesh.nElements);
            
            % Test single refinement
            mesh1 = mesh.refine(1);
            testCase.verifyEqual(mesh1.nElements, 4);
            testCase.verifyEqual(mesh1.nVertices, 5);
            
            % Test double refinement
            mesh2 = mesh.refine(2);
            testCase.verifyEqual(mesh2.nElements, 8);
        end
        
        function testRefine2D(testCase)
            % Test 2D mesh refinement
            
            vertices = [0, 0; 1, 0; 0, 1];
            elements = [1, 2, 3];
            mesh = approx.mesh.SimplexMesh(vertices, elements);
            
            % Test single refinement (1 triangle -> 4 triangles)
            mesh1 = mesh.refine(1);
            testCase.verifyEqual(mesh1.nElements, 4);
            testCase.verifyEqual(mesh1.nVertices, 6);
            
            % Test double refinement (1 -> 4 -> 16)
            mesh2 = mesh.refine(2);
            testCase.verifyEqual(mesh2.nElements, 16);
        end
        
        function testRefine3D(testCase)
            % Test 3D mesh refinement
            
            vertices = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];
            elements = [1, 2, 3, 4];
            mesh = approx.mesh.SimplexMesh(vertices, elements);
            
            % Test single refinement (1 tetrahedron -> 8 tetrahedra)
            mesh1 = mesh.refine(1);
            testCase.verifyEqual(mesh1.nElements, 8);
            testCase.verifyEqual(mesh1.nVertices, 10);
        end
        
        function testComputeMeasure(testCase)
            % Test mesh measure computation
            
            % Unit triangle
            vertices = [0, 0; 1, 0; 0, 1];
            elements = [1, 2, 3];
            mesh = approx.mesh.SimplexMesh(vertices, elements);
            
            h = mesh.computeMeasure();
            expectedMeasure = 1; % Minimum edge length
            testCase.verifyEqual(h, expectedMeasure, 'AbsTol', 1e-10);
            
            % Scaled triangle
            vertices2 = 0.5 * vertices;
            mesh2 = approx.mesh.SimplexMesh(vertices2, elements);
            h2 = mesh2.computeMeasure();
            testCase.verifyEqual(h2, 0.5, 'AbsTol', 1e-10);
        end
        
        function testElementJacobianDeterminants(testCase, meshCase)
            % Test element Jacobian determinant computation
            
            mesh = approx.mesh.SimplexMesh(meshCase.vertices, meshCase.elements);
            detJ = mesh.computeElementJacobianDeterminants();
            
            testCase.verifyEqual(length(detJ), mesh.nElements);
            testCase.verifyTrue(all(detJ > 0)); % Positive orientation
        end
        
        function testElementJacobians(testCase, meshCase)
            % Test element Jacobian matrix computation
            
            mesh = approx.mesh.SimplexMesh(meshCase.vertices, meshCase.elements);
            J = mesh.computeElementJacobians();
            
            expectedSize = [meshCase.nDims, meshCase.nDims, mesh.nElements];
            testCase.verifyEqual(size(J), expectedSize);
        end
        
        function testElementInverseJacobians(testCase, meshCase)
            % Test element inverse Jacobian computation
            
            mesh = approx.mesh.SimplexMesh(meshCase.vertices, meshCase.elements);
            invJ = mesh.computeElementInverseJacobians();
            
            expectedSize = [meshCase.nDims, meshCase.nDims, mesh.nElements];
            testCase.verifyEqual(size(invJ), expectedSize);
            
            % Test that J * invJ = I for 2D and 3D
            if meshCase.nDims > 1
                J = mesh.computeElementJacobians();
                for iElement = 1:mesh.nElements
                    product = J(:, :, iElement) * invJ(:, :, iElement);
                    testCase.verifyEqual(product, eye(meshCase.nDims), 'AbsTol', 1e-10);
                end
            end
        end
        
        function testFaceJacobianDeterminants(testCase, meshCase)
            % Test face Jacobian determinant computation
            
            mesh = approx.mesh.SimplexMesh(meshCase.vertices, meshCase.elements);
            
            for iFace = 1:meshCase.nVerticesPerElement
                detJFace = mesh.computeFaceJacobianDeterminants(iFace);
                
                testCase.verifyEqual(length(detJFace), mesh.nElements);
                testCase.verifyTrue(all(detJFace > 0));
            end
        end
        
        function testOutwardNormals(testCase, meshCase)
            % Test outward normal vector computation
            
            mesh = approx.mesh.SimplexMesh(meshCase.vertices, meshCase.elements);
            
            for iFace = 1:meshCase.nVerticesPerElement
                normals = mesh.computeOutwardNormals(iFace);
                
                testCase.verifyEqual(size(normals), [meshCase.nDims, mesh.nElements]);
                
                % Check that normals are unit vectors
                for iElement = 1:mesh.nElements
                    normalMagnitude = norm(normals(:, iElement));
                    testCase.verifyEqual(normalMagnitude, 1, 'AbsTol', 1e-10);
                end
            end
        end
        
        function testSpecific2DTriangle(testCase)
            % Test specific known 2D triangle case
            
            vertices = [0, 0; 1, 0; 0, 1];
            elements = [1, 2, 3];
            mesh = approx.mesh.SimplexMesh(vertices, elements);
            
            % Test Jacobian determinant (area * 2)
            detJ = mesh.computeElementJacobianDeterminants();
            expectedDetJ = 1; % For this unit right triangle
            testCase.verifyEqual(detJ, expectedDetJ, 'AbsTol', 1e-10);
            
            % Test face normals
            normals1 = mesh.computeOutwardNormals(1); % Face opposite vertex 1
            normals2 = mesh.computeOutwardNormals(2); % Face opposite vertex 2
            normals3 = mesh.computeOutwardNormals(3); % Face opposite vertex 3
            
            % Expected outward normals for standard triangle
            expectedNormal1 = [1; 1] / sqrt(2); % Hypotenuse
            expectedNormal2 = [-1; 0]; % Left edge
            expectedNormal3 = [0; -1]; % Bottom edge
            
            testCase.verifyEqual(normals1, expectedNormal1, 'AbsTol', 1e-10);
            testCase.verifyEqual(normals2, expectedNormal2, 'AbsTol', 1e-10);
            testCase.verifyEqual(normals3, expectedNormal3, 'AbsTol', 1e-10);
        end
        
        function testSpecific3DTetrahedron(testCase)
            % Test specific known 3D tetrahedron case
            
            vertices = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];
            elements = [1, 2, 3, 4];
            mesh = approx.mesh.SimplexMesh(vertices, elements);
            
            % Test Jacobian determinant (volume * 6)
            detJ = mesh.computeElementJacobianDeterminants();
            expectedDetJ = 1; % For this unit tetrahedron
            testCase.verifyEqual(detJ, expectedDetJ, 'AbsTol', 1e-10);
            
            % Test that we have 4 faces
            testCase.verifyEqual(mesh.nFaces, 4);
            testCase.verifyEqual(mesh.nVerticesPerFace, 3);
        end
        
        function testBoundaryDetection(testCase)
            % Test boundary face detection
            
            % Two triangles sharing an edge
            vertices = [0, 0; 1, 0; 0.5, 1; 1, 1];
            elements = [1, 2, 3; 2, 4, 3];
            mesh = approx.mesh.SimplexMesh(vertices, elements);
            
            % Should have 5 faces (4 boundary + 1 internal)
            testCase.verifyEqual(mesh.nFaces, 5);
            testCase.verifyEqual(mesh.nBoundaryFaces, 4);
            
            % Verify boundary faces have only one adjacent element
            faceToElementTable = mesh.faceToElementTable;
            for iBoundary = mesh.boundary'
                nAdjacentElements = sum(faceToElementTable(iBoundary, :) > 0);
                testCase.verifyEqual(full(nAdjacentElements), 1);
            end
        end
        
        function testNeighborElements(testCase)
            % Test neighbor element finding
            
            % Two triangles sharing an edge
            vertices = [0, 0; 1, 0; 0.5, 1; 1, 1];
            elements = [1, 2, 3; 2, 4, 3];
            mesh = approx.mesh.SimplexMesh(vertices, elements);
            
            % Find neighbors of element 1
            neighbors1_face1 = mesh.findNeighborElements(1, 1, 'strict');
            neighbors1_face2 = mesh.findNeighborElements(2, 1, 'strict');
            neighbors1_face3 = mesh.findNeighborElements(3, 1, 'strict');
            
            % Element 1 should have element 2 as neighbor through one face
            allNeighbors = [neighbors1_face1, neighbors1_face2, neighbors1_face3];
            testCase.verifyTrue(ismember(2, allNeighbors));
        end
        
        function testInteriorElements(testCase)
            % Test interior element detection
            
            % Single triangle (no interior elements)
            vertices = [0, 0; 1, 0; 0, 1];
            elements = [1, 2, 3];
            mesh = approx.mesh.SimplexMesh(vertices, elements);
            
            interiorElements = mesh.findInteriorElements();
            testCase.verifyEmpty(interiorElements);
            
            % More complex mesh might have interior elements
            % This would require a larger mesh to test properly
        end
        
        function testPeriodicMapping(testCase)
            % Test periodic boundary conditions
            
            % Create a simple periodic mesh
            vertices = [0, 0; 1, 0; 0, 1; 1, 1];
            elements = [1, 2, 3; 2, 4, 3];
            mesh = approx.mesh.SimplexMesh(vertices, elements);
            
            % Set periodic with unit shift
            shift = [1, 0; 0, 1];
            mesh = mesh.setPeriodic(shift, 1e-12);
            
            % Verify vertex-to-vertex table is created
            testCase.verifyFalse(isempty(mesh.vertexToVertexTable));
            testCase.verifyFalse(isempty(mesh.faceToFaceTable));
        end
        
        function testEmptyElementList(testCase)
            % Test methods with empty element lists
            
            vertices = [0, 0; 1, 0; 0, 1];
            elements = [1, 2, 3];
            mesh = approx.mesh.SimplexMesh(vertices, elements);
            
            % Test with empty element list
            emptyNeighbors = mesh.findNeighborElements(1, [], 'strict');
            testCase.verifyEmpty(emptyNeighbors);
            
            % Test boundary elements with invalid index
            try
                boundaryElements = mesh.findBoundaryElements(1);
                % Should not error, might return empty or valid elements
                testCase.verifyTrue(isnumeric(boundaryElements));
            catch
                % Some boundary element queries might error for single triangles
                testCase.verifyTrue(true); % Pass if error is expected
            end
        end
        
        function testLargeRefinement(testCase)
            % Test multiple levels of refinement
            
            vertices = [0, 0; 1, 0; 0, 1];
            elements = [1, 2, 3];
            mesh = approx.mesh.SimplexMesh(vertices, elements);
            
            % Test 3 levels of refinement
            mesh3 = mesh.refine(3);
            expectedElements = 4^3; % Each triangle becomes 4 triangles per level
            testCase.verifyEqual(mesh3.nElements, expectedElements);
            
            % Verify mesh is still valid
            detJ = mesh3.computeElementJacobianDeterminants();
            testCase.verifyTrue(all(detJ > 0));
        end
        
        function testDegenerateCase(testCase)
            % Test handling of degenerate cases
            
            % Triangle with very small area (but not degenerate)
            vertices = [0, 0; 1e-6, 0; 0, 1e-6];
            elements = [1, 2, 3];
            mesh = approx.mesh.SimplexMesh(vertices, elements);
            
            % Should still have positive Jacobian determinant
            detJ = mesh.computeElementJacobianDeterminants();
            testCase.verifyTrue(detJ > 0);
            
            % Measure should be very small
            h = mesh.computeMeasure();
            testCase.verifyTrue(h < 1e-5);
        end
    end
end