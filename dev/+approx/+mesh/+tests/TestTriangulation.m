classdef TestTriangulation < matlab.unittest.TestCase
    % TESTTRIANGULATION Test suite for Triangulation class.
    %
    %   This test suite validates the functionality of the Triangulation
    %   class following Richard K. Johnson's MATLAB style guidelines,
    %   including proper constructor usage with domain geometries,
    %   coordinate transformations, mesh refinement, and geometric
    %   computations.
    
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
            % Test Triangulation constructor for different dimensions.
            
            mesh = approx.mesh.Triangulation(meshCase.vertices, meshCase.elements);
            
            testCase.verifyEqual(mesh.Vertices, meshCase.vertices);
            testCase.verifyEqual(mesh.Elements, meshCase.elements);
            testCase.verifyEqual(mesh.NDims, meshCase.nDims);
            testCase.verifyEqual(mesh.NVerticesPerElement, meshCase.nVerticesPerElement);
            
            % Verify that faces and boundary are computed
            testCase.verifyFalse(isempty(mesh.Faces));
            testCase.verifyFalse(isempty(mesh.Boundary));
            testCase.verifyFalse(isempty(mesh.FaceToElementTable));
        end
        
        function testDependentProperties(testCase, meshCase)
            % Test all dependent properties.
            
            mesh = approx.mesh.Triangulation(meshCase.vertices, meshCase.elements);
            
            testCase.verifyEqual(mesh.NDims, size(meshCase.vertices, 2));
            testCase.verifyEqual(mesh.NVertices, size(meshCase.vertices, 1));
            testCase.verifyEqual(mesh.NElements, size(meshCase.elements, 1));
            testCase.verifyEqual(mesh.NVerticesPerElement, size(meshCase.elements, 2));
            testCase.verifyEqual(mesh.NFaces, size(mesh.Faces, 1));
            testCase.verifyEqual(mesh.NVerticesPerFace, size(mesh.Faces, 2));
            testCase.verifyEqual(mesh.NBoundaryFaces, length(mesh.Boundary));
        end
        
        function testFaceGeneration(testCase, meshCase)
            % Test face generation from elements.
            
            mesh = approx.mesh.Triangulation(meshCase.vertices, meshCase.elements);
            
            % Verify face dimensions
            expectedNVerticesPerFace = meshCase.nDims;
            testCase.verifyEqual(mesh.NVerticesPerFace, expectedNVerticesPerFace);
            
            % Verify face-to-element connectivity
            testCase.verifyEqual(size(mesh.FaceToElementTable, 1), mesh.NFaces);
            testCase.verifyEqual(size(mesh.FaceToElementTable, 2), mesh.NElements);
            
            % Verify boundary faces are subset of all faces
            testCase.verifyTrue(all(mesh.Boundary >= 1));
            testCase.verifyTrue(all(mesh.Boundary <= mesh.NFaces));
        end
        
        function testPositiveOrientation(testCase)
            % Test that elements have positive orientation after construction.
            
            % Create mesh with negative orientation
            vertices = [0, 0; 1, 0; 0, 1];
            elements = [1, 3, 2]; % Clockwise orientation (negative)
            
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            % Verify orientation is corrected
            detJ = mesh.computeElementJacobianDeterminants();
            testCase.verifyTrue(all(detJ > 0));
        end
        
        function testGraphify(testCase, meshCase)
            % Test graph generation from simplex mesh.
            
            mesh = approx.mesh.Triangulation(meshCase.vertices, meshCase.elements);
            graph = mesh.graphify();
            
            testCase.verifyClass(graph, 'approx.mesh.Graph');
            testCase.verifyEqual(graph.NDims, mesh.NDims);
            testCase.verifyEqual(graph.NVertices, mesh.NVertices);
            testCase.verifyEqual(graph.Vertices, mesh.Vertices);
            
            % Verify edges are created from element connectivity
            testCase.verifyTrue(graph.NEdges > 0);
            testCase.verifyEqual(size(graph.Edges, 2), 2);
        end
        
        function testCollocate(testCase)
            % Test coordinate transformation from reference to physical.
            
            % Simple 2D triangle
            vertices = [0, 0; 1, 0; 0, 1];
            elements = [1, 2, 3];
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            % Test reference coordinates
            X = [0.25; 0.25]; % Point in reference space
            Y = mesh.collocate(X, 1);
            
            testCase.verifyEqual(size(Y), [2, 1]);
            
            % Test center point (1/3, 1/3 in reference coordinates)
            Xcenter = [1/3; 1/3];
            Ycenter = mesh.collocate(Xcenter, 1);
            expectedCenter = [1/3; 1/3];
            testCase.verifyEqual(squeeze(Ycenter), expectedCenter, 'AbsTol', 1e-10);
        end
        
        function testCollocateAllElements(testCase)
            % Test collocation without specifying elements.
            
            vertices = [0, 0; 1, 0; 0.5, 1; 1, 1];
            elements = [1, 2, 3; 2, 4, 3];
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            X = [0.25; 0.25];
            Y = mesh.collocate(X);
            
            testCase.verifyEqual(size(Y), [2, 2]);
        end
        
        function testFindElementVertices(testCase, meshCase)
            % Test element vertex finding.
            
            mesh = approx.mesh.Triangulation(meshCase.vertices, meshCase.elements);
            
            % Test with specific elements
            if mesh.NElements >= 1
                V = mesh.findElementVertices(1);
                expectedSize = [meshCase.nDims, meshCase.nVerticesPerElement];
                testCase.verifyEqual(size(V), expectedSize);
            end
            
            % Test with all elements (default)
            VAll = mesh.findElementVertices();
            expectedSize = [meshCase.nDims, meshCase.nVerticesPerElement, mesh.NElements];
            testCase.verifyEqual(size(VAll), expectedSize);
        end
        
        function testRefine1D(testCase)
            % Test 1D mesh refinement.
            
            vertices = [0; 1; 2];
            elements = [1, 2; 2, 3];
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            % Test no refinement
            mesh0 = mesh.refine(0);
            testCase.verifyEqual(mesh0.NElements, mesh.NElements);
            
            % Test single refinement
            mesh1 = mesh.refine(1);
            testCase.verifyEqual(mesh1.NElements, 4);
            testCase.verifyEqual(mesh1.NVertices, 5);
            
            % Test double refinement
            mesh2 = mesh.refine(2);
            testCase.verifyEqual(mesh2.NElements, 8);
        end
        
        function testRefine2D(testCase)
            % Test 2D mesh refinement.
            
            vertices = [0, 0; 1, 0; 0, 1];
            elements = [1, 2, 3];
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            % Test single refinement (1 triangle -> 4 triangles)
            mesh1 = mesh.refine(1);
            testCase.verifyEqual(mesh1.NElements, 4);
            testCase.verifyEqual(mesh1.NVertices, 6);
            
            % Test double refinement (1 -> 4 -> 16)
            mesh2 = mesh.refine(2);
            testCase.verifyEqual(mesh2.NElements, 16);
        end
        
        function testRefine3D(testCase)
            % Test 3D mesh refinement.
            
            vertices = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];
            elements = [1, 2, 3, 4];
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            % Test single refinement (1 tetrahedron -> 8 tetrahedra)
            mesh1 = mesh.refine(1);
            testCase.verifyEqual(mesh1.NElements, 8);
            testCase.verifyEqual(mesh1.NVertices, 10);
        end
        
        function testComputeMeasure(testCase)
            % Test mesh measure computation.
            
            % Unit triangle
            vertices = [0, 0; 1, 0; 0, 1];
            elements = [1, 2, 3];
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            h = mesh.computeMeasure();
            expectedMeasure = 1; % Minimum edge length
            testCase.verifyEqual(h, expectedMeasure, 'AbsTol', 1e-10);
            
            % Scaled triangle
            vertices2 = 0.5 * vertices;
            mesh2 = approx.mesh.Triangulation(vertices2, elements);
            h2 = mesh2.computeMeasure();
            testCase.verifyEqual(h2, 0.5, 'AbsTol', 1e-10);
        end
        
        function testElementJacobianDeterminants(testCase, meshCase)
            % Test element Jacobian determinant computation.
            
            mesh = approx.mesh.Triangulation(meshCase.vertices, meshCase.elements);
            detJ = mesh.computeElementJacobianDeterminants();
            
            testCase.verifyEqual(length(detJ), mesh.NElements);
            testCase.verifyTrue(all(detJ > 0)); % Positive orientation
        end
        
        function testElementJacobians(testCase, meshCase)
            % Test element Jacobian matrix computation.
            
            mesh = approx.mesh.Triangulation(meshCase.vertices, meshCase.elements);
            J = mesh.computeElementJacobians();
            
            expectedSize = [meshCase.nDims, meshCase.nDims, mesh.NElements];
            testCase.verifyEqual(size(J), expectedSize);
        end
        
        function testElementInverseJacobians(testCase, meshCase)
            % Test element inverse Jacobian computation.
            
            mesh = approx.mesh.Triangulation(meshCase.vertices, meshCase.elements);
            invJ = mesh.computeElementInverseJacobians();
            
            expectedSize = [meshCase.nDims, meshCase.nDims, mesh.NElements];
            testCase.verifyEqual(size(invJ), expectedSize);
            
            % Test that J * invJ = I for 2D and 3D
            if meshCase.nDims > 1
                J = mesh.computeElementJacobians();
                for iElement = 1:mesh.NElements
                    product = J(:, :, iElement) * invJ(:, :, iElement);
                    testCase.verifyEqual(product, eye(meshCase.nDims), 'AbsTol', 1e-10);
                end
            end
        end
        
        function testFaceJacobianDeterminants(testCase, meshCase)
            % Test face Jacobian determinant computation.
            
            mesh = approx.mesh.Triangulation(meshCase.vertices, meshCase.elements);
            
            detJFace = mesh.computeAllFaceJacobianDeterminants();
            
            % Verify we get a cell array with one entry per face type
            testCase.verifyClass(detJFace, 'cell');
            testCase.verifyEqual(length(detJFace), meshCase.nVerticesPerElement);
            
            for iFace = 1:meshCase.nVerticesPerElement
                testCase.verifyEqual(length(detJFace{iFace}), mesh.NElements);
                testCase.verifyTrue(all(detJFace{iFace} > 0));
            end
        end
        
        function testOutwardNormals(testCase, meshCase)
            % Test outward normal vector computation.
            
            mesh = approx.mesh.Triangulation(meshCase.vertices, meshCase.elements);
            
            outN = mesh.computeAllOutwardNormals();
            
            % Verify we get a cell array with one entry per face type
            testCase.verifyClass(outN, 'cell');
            testCase.verifyEqual(length(outN), meshCase.nVerticesPerElement);
            
            for iFace = 1:meshCase.nVerticesPerElement
                testCase.verifyEqual(size(outN{iFace}), [meshCase.nDims, mesh.NElements]);
                
                % Check that normals are unit vectors
                for iElement = 1:mesh.NElements
                    normalMagnitude = norm(outN{iFace}(:, iElement));
                    testCase.verifyEqual(normalMagnitude, 1, 'AbsTol', 1e-10);
                end
            end
        end
        
        function testSpecific2DTriangle(testCase)
            % Test specific known 2D triangle case.
            
            vertices = [0, 0; 1, 0; 0, 1];
            elements = [1, 2, 3];
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            % Test Jacobian determinant (area * 2)
            detJ = mesh.computeElementJacobianDeterminants();
            expectedDetJ = 1; % For this unit right triangle
            testCase.verifyEqual(detJ, expectedDetJ, 'AbsTol', 1e-10);
            
            % Test face normals
            outN = mesh.computeAllOutwardNormals();
            normals1 = outN{1}; % Face opposite vertex 1
            normals2 = outN{2}; % Face opposite vertex 2  
            normals3 = outN{3}; % Face opposite vertex 3
            
            % Expected outward normals for standard triangle
            expectedNormal1 = [1; 1] / sqrt(2); % Hypotenuse
            expectedNormal2 = [-1; 0]; % Left edge
            expectedNormal3 = [0; -1]; % Bottom edge
            
            testCase.verifyEqual(normals1, expectedNormal1, 'AbsTol', 1e-10);
            testCase.verifyEqual(normals2, expectedNormal2, 'AbsTol', 1e-10);
            testCase.verifyEqual(normals3, expectedNormal3, 'AbsTol', 1e-10);
        end
        
        function testSpecific3DTetrahedron(testCase)
            % Test specific known 3D tetrahedron case.
            
            vertices = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];
            elements = [1, 2, 3, 4];
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            % Test Jacobian determinant (volume * 6)
            detJ = mesh.computeElementJacobianDeterminants();
            expectedDetJ = 1; % For this unit tetrahedron
            testCase.verifyEqual(detJ, expectedDetJ, 'AbsTol', 1e-10);
            
            % Test that we have 4 faces
            testCase.verifyEqual(mesh.NFaces, 4);
            testCase.verifyEqual(mesh.NVerticesPerFace, 3);
        end
        
        function testBoundaryDetection(testCase)
            % Test boundary face detection.
            
            % Two triangles sharing an edge
            vertices = [0, 0; 1, 0; 0.5, 1; 1, 1];
            elements = [1, 2, 3; 2, 4, 3];
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            % Should have 5 faces (4 boundary + 1 internal)
            testCase.verifyEqual(mesh.NFaces, 5);
            testCase.verifyEqual(mesh.NBoundaryFaces, 4);
            
            % Verify boundary faces have only one adjacent element
            faceToElementTable = mesh.FaceToElementTable;
            for iBoundary = mesh.Boundary'
                nAdjacentElements = sum(faceToElementTable(iBoundary, :) > 0);
                testCase.verifyEqual(full(nAdjacentElements), 1);
            end
        end
        
        function testNeighborElements(testCase)
            % Test neighbor element finding using PolytopalMesh methods.
            
            % Two triangles sharing an edge
            vertices = [0, 0; 1, 0; 0.5, 1; 1, 1];
            elements = [1, 2, 3; 2, 4, 3];
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            % Find neighbors of element 1 for each local face
            [JE1, JFE1] = mesh.findNeighborElements(1, 1, 'dirichlet');
            [JE2, JFE2] = mesh.findNeighborElements(1, 2, 'dirichlet');
            [JE3, JFE3] = mesh.findNeighborElements(1, 3, 'dirichlet');
            
            % Element 1 should have element 2 as neighbor through one face
            allNeighbors = [JE1, JE2, JE3];
            testCase.verifyTrue(ismember(2, allNeighbors));
        end
        
        function testInternalElements(testCase)
            % Test internal element detection.
            
            % Single triangle (no internal elements)
            vertices = [0, 0; 1, 0; 0, 1];
            elements = [1, 2, 3];
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            internalElements = mesh.getInternalElements();
            testCase.verifyEmpty(internalElements);
        end
        
        function testBoundaryElements(testCase)
            % Test boundary element detection.
            
            % Two triangles sharing an edge
            vertices = [0, 0; 1, 0; 0.5, 1; 1, 1];
            elements = [1, 2, 3; 2, 4, 3];
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            % All elements should be boundary elements for this simple mesh
            boundaryElements = mesh.getBoundaryElements();
            testCase.verifyEqual(sort(boundaryElements), [1, 2]);
        end
        
        function testPeriodicMapping(testCase)
            % Test periodic boundary conditions.
            
            % Create a simple periodic mesh
            vertices = [0, 0; 1, 0; 0, 1; 1, 1];
            elements = [1, 2, 3; 2, 4, 3];
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            % Set periodic with unit shift
            shift = [1, 0; 0, 1];
            mesh = mesh.setPeriodic(shift, 1e-12);
            
            % Verify vertex-to-vertex table is created
            testCase.verifyFalse(isempty(mesh.VertexToVertexTables));
            testCase.verifyFalse(isempty(mesh.FaceToFaceTable));
        end
        
        function testLargeRefinement(testCase)
            % Test multiple levels of refinement.
            
            vertices = [0, 0; 1, 0; 0, 1];
            elements = [1, 2, 3];
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            % Test 3 levels of refinement
            mesh3 = mesh.refine(3);
            expectedElements = 4^3; % Each triangle becomes 4 triangles per level
            testCase.verifyEqual(mesh3.NElements, expectedElements);
            
            % Verify mesh is still valid
            detJ = mesh3.computeElementJacobianDeterminants();
            testCase.verifyTrue(all(detJ > 0));
        end
        
        function testDegenerateCase(testCase)
            % Test handling of degenerate cases.
            
            % Triangle with very small area (but not degenerate)
            vertices = [0, 0; 1e-6, 0; 0, 1e-6];
            elements = [1, 2, 3];
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            % Should still have positive Jacobian determinant
            detJ = mesh.computeElementJacobianDeterminants();
            testCase.verifyTrue(detJ > 0);
            
            % Measure should be very small
            h = mesh.computeMeasure();
            testCase.verifyTrue(h < 1e-5);
        end
        
        function testTransformCache(testCase)
            % Test the transform caching system.
            
            vertices = [0, 0; 1, 0; 0, 1];
            elements = [1, 2, 3];
            mesh = approx.mesh.Triangulation(vertices, elements);
            
            % Initially, transform cache should be empty
            testCase.verifyEqual(mesh.Status, 0b110000);
            
            % Compute Jacobians - should set the appropriate flag
            J = mesh.computeElementJacobians();
            testCase.verifyTrue(bitand(mesh.Status, 0b100000) > 0);
            testCase.verifyFalse(isempty(mesh.Transform.EJac));
            
            % Compute determinants - should set another flag
            detJ = mesh.computeElementJacobianDeterminants();
            testCase.verifyTrue(bitand(mesh.Status, 0b010000) > 0);
            testCase.verifyFalse(isempty(mesh.Transform.detEJac));
            
            % Reset cache
            mesh.resetTransform();
            testCase.verifyEqual(mesh.Status, 0b000000);
            testCase.verifyEmpty(mesh.Transform.EJac);
        end
    end
end