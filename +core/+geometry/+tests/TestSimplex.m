classdef TestSimplex < matlab.unittest.TestCase
    % TESTSIMPLEX Unit tests for the Simplex class.
    
    properties (Constant)
        Tolerance = 1e-14
    end
    
    methods (Test)
        function testConstructor(testCase)
            vertices1D = [0, 1];
            simplex1D = core.geometry.Simplex(vertices1D);
            testCase.verifyEqual(simplex1D.NDims, 1);
            testCase.verifyEqual(simplex1D.Vertices, vertices1D);
            testCase.verifySize(simplex1D.Faces, [1, 2]);
            
            vertices2D = [0, 1, 0; 0, 0, 1];
            simplex2D = core.geometry.Simplex(vertices2D);
            testCase.verifyEqual(simplex2D.NDims, 2);
            testCase.verifyEqual(simplex2D.Vertices, vertices2D);
            testCase.verifySize(simplex2D.Faces, [1, 3]);
            
            vertices3D = [0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
            simplex3D = core.geometry.Simplex(vertices3D);
            testCase.verifyEqual(simplex3D.NDims, 3);
            testCase.verifyEqual(simplex3D.Vertices, vertices3D);
            testCase.verifySize(simplex3D.Faces, [1, 4]);
        end
        
        function testConstructorErrors(testCase)
            invalidVertices2D = [0, 1, 0, 0; 0, 0, 1, 0];
            testCase.verifyError(@() core.geometry.Simplex(invalidVertices2D), ...
                'core:geometry:Simplex:InvalidInput');
                
            dependentVertices2D = [0, 1, 2; 0, 1, 2];
            testCase.verifyError(@() core.geometry.Simplex(dependentVertices2D), ...
                'core:geometry:Simplex:InvalidInput');
        end
        
        function testInheritance(testCase)
            vertices2D = [0, 1, 0; 0, 0, 1];
            simplex = core.geometry.Simplex(vertices2D);
            testCase.verifyTrue(isa(simplex, 'core.geometry.Polytope'));
            testCase.verifyTrue(isa(simplex, 'core.geometry.Geometry'));
        end
        
        function testMagnitude(testCase)
            vertices1D = [0, 3];
            simplex1D = core.geometry.Simplex(vertices1D);
            testCase.verifyEqual(simplex1D.magnitude(), 3);
            
            vertices2D = [0, 1, 0; 0, 0, 1];
            simplex2D = core.geometry.Simplex(vertices2D);
            testCase.verifyEqual(simplex2D.magnitude(), 0.5);
            
            vertices2D_2 = [0, 2, 0; 0, 0, 3];
            simplex2D_2 = core.geometry.Simplex(vertices2D_2);
            testCase.verifyEqual(simplex2D_2.magnitude(), 3);
            
            vertices3D = [0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
            simplex3D = core.geometry.Simplex(vertices3D);
            testCase.verifyEqual(simplex3D.magnitude(), 1/6);
        end
        
        function testIsInside(testCase)
            vertices2D = [0, 1, 0; 0, 0, 1];
            simplex = core.geometry.Simplex(vertices2D);
            
            interiorPoints = [0.2, 0.3, 0.4; 0.2, 0.3, 0.1];
            testCase.verifyEqual(simplex.isInside(interiorPoints), [true, true, true]);
            
            exteriorPoints = [-0.1, 1.1, 0.5; 0.1, 0.1, 1.1];
            testCase.verifyEqual(simplex.isInside(exteriorPoints), [false, false, false]);
            
            boundaryPoints = [0, 1, 0, 0.5; 0, 0, 1, 0.5];
            testCase.verifyEqual(simplex.isInside(boundaryPoints), [false, false, false, false]);
            
            testCase.verifyError(@() simplex.isInside([1, 2, 3; 2, 3, 4; 5, 6, 7]), ?MException);
        end
        
        function testIsOnBoundary(testCase)
            vertices2D = [0, 1, 0; 0, 0, 1];
            simplex = core.geometry.Simplex(vertices2D);
            
            interiorPoints = [0.2, 0.3, 0.4; 0.2, 0.3, 0.1];
            testCase.verifyEqual(simplex.isOnBoundary(interiorPoints), [false, false, false]);
            
            exteriorPoints = [-0.1, 1.1, 0.5; 0.1, 0.1, 1.1];
            testCase.verifyEqual(simplex.isOnBoundary(exteriorPoints), [false, false, false]);
            
            boundaryPoints = [0, 1, 0, 0.5; 0, 0, 1, 0.5];
            testCase.verifyEqual(simplex.isOnBoundary(boundaryPoints), [true, true, true, true]);
            
            edgePoints = [0.5, 0.5, 0; 0, 0.5, 0.5];
            testCase.verifyEqual(simplex.isOnBoundary(edgePoints), [true, true, true]);
        end
        
        function testBarycentricCoordinateTransformation(testCase)
            vertices2D = [0, 1, 0; 0, 0, 1];
            simplex = core.geometry.Simplex(vertices2D);
            
            baryCoords = simplex.cartesianToBarycentric(vertices2D);
            expected = eye(3);
            testCase.verifyEqual(baryCoords, expected);
            
            point = [0.25; 0.25];
            baryCoords = simplex.cartesianToBarycentric(point);
            expected = [0.5; 0.25; 0.25];
            testCase.verifyEqual(baryCoords, expected, 'AbsTol', testCase.Tolerance);
            
            testBaryCoords = [0.2; 0.3; 0.5];
            point = simplex.barycentricToCartesian(testBaryCoords);
            baryCoords = simplex.cartesianToBarycentric(point);
            testCase.verifyEqual(baryCoords, testBaryCoords, 'AbsTol', testCase.Tolerance);
            
            points = [0.2, 0.6, 0.3; 0.4, 0.2, 0.5];
            baryCoords = simplex.cartesianToBarycentric(points);
            recoveredPoints = simplex.barycentricToCartesian(baryCoords);
            testCase.verifyEqual(recoveredPoints, points, 'AbsTol', testCase.Tolerance);
        end
        
        function testInvalidBarycentricCoordinates(testCase)
            vertices2D = [0, 1, 0; 0, 0, 1];
            simplex = core.geometry.Simplex(vertices2D);
            
            invalidCoords = [0.2; 0.3; 0.3];
            testCase.verifyError(@() simplex.barycentricToCartesian(invalidCoords), ...
                'core:geometry:Simplex:InvalidInput');
            
            negativeCoords = [0.5; 0.8; -0.3];
            testCase.verifyError(@() simplex.barycentricToCartesian(negativeCoords), ...
                'core:geometry:Simplex:InvalidInput');
            
            invalidSizeCoords = [0.5; 0.5];
            testCase.verifyError(@() simplex.barycentricToCartesian(invalidSizeCoords), ...
                'core:geometry:Simplex:DimensionMismatch');
        end
        
        function testMultidimensionalInput(testCase)
            vertices2D = [0, 1, 0; 0, 0, 1];
            simplex = core.geometry.Simplex(vertices2D);
            
            points = zeros(2, 2, 2);
            points(:, 1, 1) = [0.3; 0.3];
            points(:, 2, 1) = [0.4; 0.2];
            points(:, 1, 2) = [0.1; 0.1];
            points(:, 2, 2) = [0.2; 0.6];
            
            baryCoords = simplex.cartesianToBarycentric(points);
            testCase.verifySize(baryCoords, [3, 2, 2]);
            
            recoveredPoints = simplex.barycentricToCartesian(baryCoords);
            testCase.verifyEqual(recoveredPoints, points, 'AbsTol', testCase.Tolerance);
            
            for i = 1:2
                for j = 1:2
                    point = points(:, i, j);
                    baryCoord = baryCoords(:, i, j);
                    testCase.verifyEqual(sum(baryCoord), 1, 'AbsTol', testCase.Tolerance);
                    testCase.verifyEqual(simplex.Vertices * baryCoord, point, 'AbsTol', testCase.Tolerance);
                end
            end
        end
    end
end