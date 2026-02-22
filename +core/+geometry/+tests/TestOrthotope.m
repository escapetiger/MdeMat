classdef TestOrthotope < matlab.unittest.TestCase
    % TESTORTHOTOPE Unit tests for the Orthotope class.
    %
    % This test class validates the functionality of the Orthotope class,
    % including constructor behavior, geometric operations, and inheritance
    % properties. Tests cover various dimensions and edge cases to ensure
    % robust implementation.
    %
    % See also:
    %   core.geometry.Orthotope
    
    properties (Constant)
        Tolerance = 1e-10 % Numerical tolerance for floating-point comparisons
    end
    
    methods (Test)
        function testConstructor(testCase)
            % TESTCONSTRUCTOR Test orthotope constructor with bbox format.
            %
            % Verifies that the constructor properly handles the bbox input
            % format [a1, b1, a2, b2, ...] and creates correct orthotopes
            % with proper dimensions, bounds, vertices, and faces.
            
            % Test 1D orthotope (interval)
            rect1D = core.geometry.Orthotope([0, 1]);
            testCase.verifyEqual(rect1D.NDims, 1);
            testCase.verifyEqual(rect1D.Lower, [0]);
            testCase.verifyEqual(rect1D.Upper, [1]);
            testCase.verifySize(rect1D.Vertices, [1, 2]);
            testCase.verifySize(rect1D.Faces, [1, 2]);
            
            % Test 2D orthotope (rectangle)
            rect2D = core.geometry.Orthotope([1, 4, 2, 5]);
            testCase.verifyEqual(rect2D.NDims, 2);
            testCase.verifyEqual(rect2D.Lower, [1, 2]);
            testCase.verifyEqual(rect2D.Upper, [4, 5]);
            testCase.verifySize(rect2D.Vertices, [2, 4]);
            testCase.verifySize(rect2D.Faces, [1, 4]);
            
            % Test 3D orthotope (box)
            rect3D = core.geometry.Orthotope([1, 4, 2, 5, 3, 6]);
            testCase.verifyEqual(rect3D.NDims, 3);
            testCase.verifyEqual(rect3D.Lower, [1, 2, 3]);
            testCase.verifyEqual(rect3D.Upper, [4, 5, 6]);
            testCase.verifySize(rect3D.Vertices, [3, 8]);
            testCase.verifySize(rect3D.Faces, [1, 6]);
            
            % Test unit cube
            unitCube = core.geometry.Orthotope([0, 1, 0, 1, 0, 1]);
            testCase.verifyEqual(unitCube.NDims, 3);
            testCase.verifyEqual(unitCube.Lower, [0, 0, 0]);
            testCase.verifyEqual(unitCube.Upper, [1, 1, 1]);
        end
        
        function testConstructorValidation(testCase)
            % TESTCONSTRUCTORVALIDATION Test constructor input validation.
            %
            % Verifies that the constructor properly validates input and
            % throws appropriate errors for invalid inputs.
            
            % Test empty input
            testCase.verifyError(@() core.geometry.Orthotope(), ?MException);
            
            % Test non-vector input
            testCase.verifyError(@() core.geometry.Orthotope([1, 2; 3, 4]), ?MException);
            
            % Test odd number of elements
            testCase.verifyError(@() core.geometry.Orthotope([1, 2, 3]), ?MException);
            
            % Test insufficient elements
            testCase.verifyError(@() core.geometry.Orthotope([1]), ?MException);
            
            % Test lower > upper bounds
            testCase.verifyError(@() core.geometry.Orthotope([3, 2]), ?MException);
            testCase.verifyError(@() core.geometry.Orthotope([1, 2, 4, 3]), ?MException);
            
            % Test valid edge cases
            % Equal bounds (degenerate orthotope)
            degenerateRect = core.geometry.Orthotope([1, 1, 2, 2]);
            testCase.verifyEqual(degenerateRect.Lower, [1, 2]);
            testCase.verifyEqual(degenerateRect.Upper, [1, 2]);
            
            % Negative bounds
            negativeRect = core.geometry.Orthotope([-2, -1, -1, 0]);
            testCase.verifyEqual(negativeRect.Lower, [-2, -1]);
            testCase.verifyEqual(negativeRect.Upper, [-1, 0]);
        end
        
        function testInheritance(testCase)
            % TESTINHERITANCE Test inheritance hierarchy.
            %
            % Verifies that Orthotope correctly inherits from Polytope
            % and Geometry classes.
            
            rect = core.geometry.Orthotope([0, 1, 0, 1]);
            testCase.verifyTrue(isa(rect, 'core.geometry.Polytope'));
            testCase.verifyTrue(isa(rect, 'core.geometry.Geometry'));
        end
        
        function testVerticesGeneration(testCase)
            % TESTVERTICESGENERATION Test vertex generation correctness.
            %
            % Verifies that vertices are correctly generated for orthotopes
            % of various dimensions and that all vertices lie within the
            % expected bounds.
            
            % Test 2D rectangle vertices
            rect2D = core.geometry.Orthotope([0, 1, 0, 1]);
            expectedVertices2D = [0, 0, 1, 1; 0, 1, 0, 1];
            
            for i = 1:4
                vertex = expectedVertices2D(:, i);
                found = false;
                for j = 1:4
                    if norm(rect2D.Vertices(:, j) - vertex) < testCase.Tolerance
                        found = true;
                        break;
                    end
                end
                testCase.verifyTrue(found, ...
                    'Expected vertex not found in generated vertices');
            end
            
            % Test 3D box vertices
            rect3D = core.geometry.Orthotope([1, 4, 2, 5, 3, 6]);
            testCase.verifySize(rect3D.Vertices, [3, 8]);
            
            % Verify all vertices are within bounds
            for i = 1:8
                vertex = rect3D.Vertices(:, i);
                testCase.verifyTrue(all(vertex >= [1, 2, 3]') && ...
                                  all(vertex <= [4, 5, 6]'), ...
                    'Vertex outside expected bounds');
                
                % Verify each coordinate is at a boundary
                testCase.verifyTrue(vertex(1) == 1 || vertex(1) == 4, ...
                    'X-coordinate not at boundary');
                testCase.verifyTrue(vertex(2) == 2 || vertex(2) == 5, ...
                    'Y-coordinate not at boundary');
                testCase.verifyTrue(vertex(3) == 3 || vertex(3) == 6, ...
                    'Z-coordinate not at boundary');
            end
        end
        
        function testMagnitude(testCase)
            % TESTMAGNITUDE Test magnitude (volume) calculation.
            %
            % Verifies that the magnitude method correctly computes the
            % measure of orthotopes in various dimensions.
            
            % Test 1D interval
            rect1D = core.geometry.Orthotope([0, 5]);
            testCase.verifyEqual(rect1D.magnitude(), 5, ...
                'AbsTol', testCase.Tolerance);
            
            % Test 2D rectangle
            rect2D = core.geometry.Orthotope([1, 3, 2, 6]);
            testCase.verifyEqual(rect2D.magnitude(), 8, ...
                'AbsTol', testCase.Tolerance);
            
            % Test 3D box
            rect3D = core.geometry.Orthotope([0, 2, 0, 3, 0, 4]);
            testCase.verifyEqual(rect3D.magnitude(), 24, ...
                'AbsTol', testCase.Tolerance);
            
            % Test degenerate case
            degenerateRect = core.geometry.Orthotope([1, 1, 0, 2]);
            testCase.verifyEqual(degenerateRect.magnitude(), 0, ...
                'AbsTol', testCase.Tolerance);
        end
        
        function testIsInside(testCase)
            % TESTISINSIDE Test interior point detection.
            %
            % Verifies that the isInside method correctly identifies points
            % strictly inside the orthotope (excluding boundaries).
            
            rect = core.geometry.Orthotope([1, 3, 2, 4]);
            
            % Test interior points
            interiorPoints = [2, 1.5, 2.9; 3, 2.5, 3.9];
            testCase.verifyEqual(rect.isInside(interiorPoints), ...
                [true, true, true]);
            
            % Test exterior points
            exteriorPoints = [0, 4, 2; 3, 3, 5];
            testCase.verifyEqual(rect.isInside(exteriorPoints), ...
                [false, false, false]);
            
            % Test boundary points (should be false)
            boundaryPoints = [1, 3, 2; 2, 4, 4];
            testCase.verifyEqual(rect.isInside(boundaryPoints), ...
                [false, false, false]);
            
            % Test dimension mismatch
            testCase.verifyError(@() rect.isInside([1, 2, 3; 2, 3, 4; 5, 6, 7]), ...
                ?MException);
        end
        
        function testIsOnBoundary(testCase)
            % TESTISONBOUNDARY Test boundary point detection.
            %
            % Verifies that the isOnBoundary method correctly identifies
            % points on the orthotope boundary.
            
            rect = core.geometry.Orthotope([1, 3, 2, 4]);
            
            % Test interior points (should be false)
            interiorPoints = [2, 1.5, 2.9; 3, 2.5, 3.9];
            testCase.verifyEqual(rect.isOnBoundary(interiorPoints), ...
                [false, false, false]);
            
            % Test exterior points (should be false)
            exteriorPoints = [0, 4, 2; 3, 3, 5];
            testCase.verifyEqual(rect.isOnBoundary(exteriorPoints), ...
                [false, false, false]);
            
            % Test boundary points
            boundaryPoints = [1, 3, 2, 2; 2, 4, 2, 4];
            testCase.verifyEqual(rect.isOnBoundary(boundaryPoints), ...
                [true, true, true, true]);
            
            % Test mixed boundary points
            mixedBoundaryPoints = [1, 1, 1, 1, 1; 1.9, 2, 3, 4, 4.1];
            expected = [false, true, true, true, false];
            testCase.verifyEqual(rect.isOnBoundary(mixedBoundaryPoints), expected);
            
            % Test corner point
            cornerPoint = [1; 2];
            testCase.verifyEqual(rect.isOnBoundary(cornerPoint), true);
        end
        
        function testDependentProperties(testCase)
            % TESTDEPENDENTPROPERTIES Test dependent property calculations.
            %
            % Verifies that dependent properties like nBoundaries and
            % outerNormals are correctly computed.
            
            % Test 2D rectangle
            rect2D = core.geometry.Orthotope([0, 1, 0, 1]);
            testCase.verifyEqual(rect2D.NFaces, 4);
            testCase.verifySize(rect2D.OutNormals, [2, 4]);
            
            % Test 3D box
            rect3D = core.geometry.Orthotope([0, 1, 0, 1, 0, 1]);
            testCase.verifyEqual(rect3D.NFaces, 6);
            testCase.verifySize(rect3D.OutNormals, [3, 6]);
            
            % Verify outer normals structure
            expectedNormals2D = [-1, 1, 0, 0; 0, 0, -1, 1];
            testCase.verifyEqual(rect2D.OutNormals, expectedNormals2D, ...
                'AbsTol', testCase.Tolerance);
        end
        
        function testPolytopeInheritance(testCase)
            % TESTPOLYTOPEINHERITANCE Test inherited Polytope functionality.
            %
            % Verifies that Orthotope correctly inherits and implements
            % Polytope methods for geometric operations.
            
            rect = core.geometry.Orthotope([0, 1, 0, 1]);
            
            insidePoint = [0.5; 0.5];
            boundaryPoint = [0; 0.5];
            outsidePoint = [2; 2];
            
            % Test point classification
            testCase.verifyTrue(rect.isInside(insidePoint));
            testCase.verifyFalse(rect.isInside(boundaryPoint));
            testCase.verifyFalse(rect.isInside(outsidePoint));
            
            testCase.verifyFalse(rect.isOnBoundary(insidePoint));
            testCase.verifyTrue(rect.isOnBoundary(boundaryPoint));
            testCase.verifyFalse(rect.isOnBoundary(outsidePoint));
        end
        
        function testVectorInputHandling(testCase)
            % TESTVECTORINPUTHANDLING Test bbox vector input variations.
            %
            % Verifies that the constructor correctly handles both row and
            % column vectors for the bbox input.
            
            % Test row vector
            rectRow = core.geometry.Orthotope([1, 2, 3, 4]);
            testCase.verifyEqual(rectRow.Lower, [1, 3]);
            testCase.verifyEqual(rectRow.Upper, [2, 4]);
            
            % Test column vector
            rectCol = core.geometry.Orthotope([1; 2; 3; 4]);
            testCase.verifyEqual(rectCol.Lower, [1, 3]);
            testCase.verifyEqual(rectCol.Upper, [2, 4]);
            
            % Verify both produce identical results
            testCase.verifyEqual(rectRow.Lower, rectCol.Lower);
            testCase.verifyEqual(rectRow.Upper, rectCol.Upper);
            testCase.verifyEqual(rectRow.NDims, rectCol.NDims);
        end
    end
end