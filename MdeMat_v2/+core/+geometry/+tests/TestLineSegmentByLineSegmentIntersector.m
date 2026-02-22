classdef TestLineSegmentByLineSegmentIntersector < matlab.unittest.TestCase
    % TESTLINESEGMENTBYLINESEGMENTINTERSECTOR Tests for LineSegmentByLineSegmentIntersector.
    %
    % This class tests the functionality of the LineSegmentByLineSegmentIntersector class,
    % which computes intersections between line segments in n-dimensional space.
    %
    % See also: LINESEGMENTBYLINESEGMENTINTERSECTOR
    
    properties
        defaultTol = sqrt(eps);
    end
    
    methods (Test)
        function testConstructor(testCase)
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(2);
            testCase.verifyEqual(intersector.nDims, 2);
            
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(3);
            testCase.verifyEqual(intersector.nDims, 3);
            
            testCase.verifyError(@() core.geometry.LineSegmentByLineSegmentIntersector(1), ...
                'core:geometry:LineSegmentByLineSegmentIntersector:InvalidInput');
        end
        
        function testBasic2DIntersection(testCase)
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(2);
            
            X1 = [0; 0];
            X2 = [1; 1];
            X3 = [0; 1];
            X4 = [1; 0];
            
            [X, TF] = intersector.intersect(X1, X2, X3, X4);
            
            testCase.verifyEqual(size(X, 2), 1);
            testCase.verifyEqual(X, [0.5; 0.5], 'AbsTol', testCase.defaultTol);
            testCase.verifyTrue(TF);
        end
        
        function testNonIntersectingLines(testCase)
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(2);
            
            X1 = [0; 0];
            X2 = [1; 1];
            X3 = [2; 0];
            X4 = [3; 1];
            
            [X, TF] = intersector.intersect(X1, X2, X3, X4);
            
            testCase.verifyEqual(size(X, 2), 0);
            testCase.verifyFalse(TF);
        end
        
        function test3DSkewLines(testCase)
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(3);
            
            X1 = [0; 0; 0];
            X2 = [1; 0; 0];
            X3 = [0; 0; 1];
            X4 = [1; 0; 1];
            
            [X, TF] = intersector.intersect(X1, X2, X3, X4);
            
            testCase.verifyEqual(size(X, 2), 0);
            testCase.verifyFalse(TF);
        end
        
        function test3DIntersectingLines(testCase)
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(3);
            
            X1 = [0; 0; 0];
            X2 = [1; 1; 1];
            X3 = [0; 0; 0];
            X4 = [1; 1; 0];
            
            [X, TF] = intersector.intersect(X1, X2, X3, X4);
            
            testCase.verifyEqual(size(X, 2), 1);
            testCase.verifyEqual(X, [0; 0; 0], 'AbsTol', testCase.defaultTol);
            testCase.verifyTrue(TF);
        end
        
        function testMultiplePairs(testCase)
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(2);
            
            X1 = [0, 0; 0, 0];
            X2 = [1, 1; 1, 2];
            X3 = [0, 2; 1, 1];
            X4 = [1, 3; 0, 3];
            
            [X, TF] = intersector.intersect(X1, X2, X3, X4);
            
            testCase.verifyEqual(size(X, 2), 1);
            testCase.verifyTrue(TF(1));
            testCase.verifyFalse(TF(2));
            testCase.verifyEqual(X, [0.5; 0.5], 'AbsTol', testCase.defaultTol);
        end
        
        function testHigherDimensions(testCase)
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(4);
            
            X1 = [0; 0; 0; 0];
            X2 = [1; 1; 1; 1];
            X3 = [0; 0; 1; 1];
            X4 = [1; 1; 0; 0];
            
            [X, TF] = intersector.intersect(X1, X2, X3, X4);
            
            testCase.verifyEqual(size(X, 2), 1);
            testCase.verifyEqual(X, [0.5; 0.5; 0.5; 0.5], 'AbsTol', testCase.defaultTol);
            testCase.verifyTrue(TF);
        end
        
        function testZeroLengthSegment(testCase)
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(2);
            
            X1 = [1; 1];
            X2 = [1; 1];
            X3 = [0; 0];
            X4 = [2; 2];
            
            [X, TF] = intersector.intersect(X1, X2, X3, X4);
            
            testCase.verifyEqual(size(X, 2), 1);
            testCase.verifyEqual(X, [1; 1], 'AbsTol', testCase.defaultTol);
            testCase.verifyTrue(TF);
        end
        
        function testBothZeroLength(testCase)
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(2);
            
            X1 = [1; 1];
            X2 = [1; 1];
            X3 = [1; 1];
            X4 = [1; 1];
            
            [X, TF] = intersector.intersect(X1, X2, X3, X4);
            
            testCase.verifyEqual(size(X, 2), 1);
            testCase.verifyEqual(X, [1; 1], 'AbsTol', testCase.defaultTol);
            testCase.verifyTrue(TF);
        end
        
        function testCoincidentSegments(testCase)
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(2);
            
            X1 = [0; 0];
            X2 = [1; 1];
            X3 = [0; 0];
            X4 = [1; 1];
            
            testCase.verifyError(@() intersector.intersect(X1, X2, X3, X4), ...
                'core:geometry:LineSegmentByLineSegmentIntersector:CoincidentSegments');
            
            X1 = [0; 0];
            X2 = [1; 1];
            X3 = [1; 1];
            X4 = [0; 0];
            
            testCase.verifyError(@() intersector.intersect(X1, X2, X3, X4), ...
                'core:geometry:LineSegmentByLineSegmentIntersector:CoincidentSegments');
        end
        
        function testCollinearOverlapping(testCase)
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(2);
            
            X1 = [0; 0];
            X2 = [2; 2];
            X3 = [1; 1];
            X4 = [3; 3];
            
            testCase.verifyError(@() intersector.intersect(X1, X2, X3, X4), ...
                'core:geometry:LineSegmentByLineSegmentIntersector:CollinearOverlap');
            
            X1 = [1; 1];
            X2 = [2; 2];
            X3 = [0; 0];
            X4 = [3; 3];
            
            testCase.verifyError(@() intersector.intersect(X1, X2, X3, X4), ...
                'core:geometry:LineSegmentByLineSegmentIntersector:CollinearOverlap');
        end
        
        function testCollinearEndpoint(testCase)
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(2);
            
            X1 = [0; 0];
            X2 = [1; 1];
            X3 = [1; 1];
            X4 = [2; 2];
            
            [X, TF] = intersector.intersect(X1, X2, X3, X4);
            
            testCase.verifyEqual(size(X, 2), 1);
            testCase.verifyEqual(X, [1; 1], 'AbsTol', testCase.defaultTol);
            testCase.verifyTrue(TF);
        end
        
        function testCollinearNoIntersection(testCase)
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(2);
            
            X1 = [0; 0];
            X2 = [1; 1];
            X3 = [2; 2];
            X4 = [3; 3];
            
            [X, TF] = intersector.intersect(X1, X2, X3, X4);
            
            testCase.verifyEqual(size(X, 2), 0);
            testCase.verifyFalse(TF);
        end
        
        function testInvalidDimensions(testCase)
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(2);
            
            X1 = [0; 0];
            X2 = [1; 1];
            X3 = [0; 0; 0];
            X4 = [1; 1; 1];
            
            testCase.verifyError(@() intersector.intersect(X1, X2, X3, X4), ...
                'core:geometry:LineSegmentByLineSegmentIntersector:InvalidInput');
        end
        
        function testMismatchedSegments(testCase)
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(2);
            
            X1 = [0, 0; 0, 0];
            X2 = [1, 1; 1, 1];
            X3 = [0; 1];
            X4 = [1; 0];
            
            testCase.verifyError(@() intersector.intersect(X1, X2, X3, X4), ...
                'core:geometry:LineSegmentByLineSegmentIntersector:InvalidInput');
        end
        
        function testNearlyParallel(testCase)
            intersector = core.geometry.LineSegmentByLineSegmentIntersector(2);
            
            X1 = [0; 0];
            X2 = [1; 1 + 1e-12];
            X3 = [0; 1];
            X4 = [1; 0];
            
            [X, TF] = intersector.intersect(X1, X2, X3, X4);
            
            testCase.verifyEqual(size(X, 2), 1);
            testCase.verifyEqual(X(1), 0.5, 'AbsTol', 1e-10);
            testCase.verifyLessThan(abs(X(2) - 0.5), 1e-10);
            testCase.verifyTrue(TF);
            
            X1 = [0; 0];
            X2 = [1; 1 + 0.1];
            X3 = [0; 1];
            X4 = [1; 0];
            
            [X, TF] = intersector.intersect(X1, X2, X3, X4);
            
            testCase.verifyEqual(size(X, 2), 1);
            testCase.verifyEqual(X(1), 0.5, 'AbsTol', 0.1);
            testCase.verifyTrue(TF);
        end
    end
end