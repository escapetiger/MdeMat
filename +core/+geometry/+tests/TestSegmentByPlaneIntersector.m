classdef TestSegmentByPlaneIntersector < matlab.unittest.TestCase
    % TESTSEGMENTBYPLANEINTERSECTOR Tests for SegmentByPlaneIntersector.
    %
    % This class tests the intersection computation between line segments
    % and planes in n-dimensional space.
    %
    % See also: SEGMENTBYPLANEINTERSECTOR
    
    methods (Test)
        function testConstructor(testCase)
            intersector = core.geometry.SegmentByPlaneIntersector(2);
            testCase.verifyEqual(intersector.NDims, 2);
            
            intersector = core.geometry.SegmentByPlaneIntersector(3);
            testCase.verifyEqual(intersector.NDims, 3);
            
            testCase.verifyError(@() core.geometry.SegmentByPlaneIntersector(1), ...
                'MATLAB:validators:mustBeGreaterThan');
        end
        
        function test2DLineIntersection(testCase)
            % Line segment from (0,0) to (2,2) with line x+y=2
            intersector = core.geometry.SegmentByPlaneIntersector(2);
            
            X1 = [0; 0];
            X2 = [2; 2];
            a = [1; 1];
            b = -2;
            
            [X, TF] = intersector.intersect(X1, X2, a, b);
            
            testCase.verifyEqual(size(X, 2), 1);
            testCase.verifyEqual(X, [1; 1], 'AbsTol', sqrt(eps));
            testCase.verifyTrue(TF);
        end
        
        function test3DLineIntersection(testCase)
            % Line segment from origin to (3,3,3) with plane x+y+z=3
            intersector = core.geometry.SegmentByPlaneIntersector(3);
            
            X1 = [0; 0; 0];
            X2 = [3; 3; 3];
            a = [1; 1; 1];
            b = -3;
            
            [X, TF] = intersector.intersect(X1, X2, a, b);
            
            testCase.verifyEqual(size(X, 2), 1);
            testCase.verifyEqual(X, [1; 1; 1], 'AbsTol', sqrt(eps));
            testCase.verifyTrue(TF);
        end
        
        function testNoIntersection(testCase)
            % Line segment from origin to (1,1,1) with plane x+y+z=4
            intersector = core.geometry.SegmentByPlaneIntersector(3);
            
            X1 = [0; 0; 0];
            X2 = [1; 1; 1];
            a = [1; 1; 1];
            b = -4;
            
            [X, TF] = intersector.intersect(X1, X2, a, b);
            
            testCase.verifyEqual(size(X, 2), 0);
            testCase.verifyFalse(TF);
        end
        
        function testParallelLine(testCase)
            % Line segment along x-axis with yz-plane
            intersector = core.geometry.SegmentByPlaneIntersector(3);
            
            X1 = [0; 0; 0];
            X2 = [1; 0; 0];
            a = [0; 1; 0];
            b = -2;
            
            [X, TF] = intersector.intersect(X1, X2, a, b);
            
            testCase.verifyEqual(size(X, 2), 0);
            testCase.verifyFalse(TF);
        end
        
        function testEndpointIntersection(testCase)
            % Line segment from origin to (1,1,1) with plane through origin
            intersector = core.geometry.SegmentByPlaneIntersector(3);
            
            X1 = [0; 0; 0];
            X2 = [1; 1; 1];
            a = [1; 1; 1];
            b = 0;
            
            [X, TF] = intersector.intersect(X1, X2, a, b);
            
            testCase.verifyEqual(size(X, 2), 1);
            testCase.verifyEqual(X, [0; 0; 0], 'AbsTol', sqrt(eps));
            testCase.verifyTrue(TF);
        end
        
        function testMultipleLineSegments(testCase)
            % Three line segments intersecting with plane x=1
            intersector = core.geometry.SegmentByPlaneIntersector(3);
            
            X1 = [0, 0, 0;
                  0, 1, 2;
                  0, 0, 0];
                  
            X2 = [2, 2, 2;
                  0, 1, 2;
                  0, 0, 0];
                  
            a = [1; 0; 0];
            b = -1;
            
            [X, TF] = intersector.intersect(X1, X2, a, b);
            
            testCase.verifyEqual(size(X, 2), 3);
            testCase.verifyTrue(all(TF));
            testCase.verifyEqual(X(1,:), [1, 1, 1], 'AbsTol', sqrt(eps));
        end
        
        function testMultiplePlanes(testCase)
            % One line segment intersecting with multiple planes
            intersector = core.geometry.SegmentByPlaneIntersector(3);
            
            X1 = [0; 0; 0];
            X2 = [4; 0; 0];
            
            a = [1, 1, 1;
                 0, 0, 0;
                 0, 0, 0];
            b = [-1, -2, -5];
            
            [X, TF] = intersector.intersect(X1, X2, a, b);
            
            testCase.verifyEqual(size(X, 2), 2);
            testCase.verifyEqual(sum(TF), 2);
            testCase.verifyEqual(X(1, TF), [1, 2], 'AbsTol', sqrt(eps));
        end
        
        function testLineInPlane(testCase)
            % Line segment lying in a plane
            intersector = core.geometry.SegmentByPlaneIntersector(3);
            
            X1 = [0; 0; 0];
            X2 = [1; 1; 0];
            a = [0; 0; 1];
            b = 0;
            
            [X, TF] = intersector.intersect(X1, X2, a, b);
            
            testCase.verifyEqual(size(X, 2), 0);
            testCase.verifyFalse(TF);
        end
        
        function testHigherDimension(testCase)
            % 4D line segment and hyperplane
            intersector = core.geometry.SegmentByPlaneIntersector(4);
            
            X1 = [0; 0; 0; 0];
            X2 = [4; 4; 4; 4];
            a = [1; 1; 1; 1];
            b = -4;
            
            [X, TF] = intersector.intersect(X1, X2, a, b);
            
            testCase.verifyEqual(size(X, 2), 1);
            testCase.verifyTrue(TF);
            testCase.verifyEqual(X, [1; 1; 1; 1], 'AbsTol', sqrt(eps));
        end
        
        function testDimensionMismatch(testCase)
            intersector = core.geometry.SegmentByPlaneIntersector(3);
            
            X1 = [0; 0; 0];
            X2 = [1; 1; 1];
            a = [1; 0];
            b = -1;
            
            testCase.verifyError(@() intersector.intersect(X1, X2, a, b), ...
                'core:geometry:SegmentByPlaneIntersector:DimensionMismatch');
        end
        
        function testMatchingSegmentPlaneCount(testCase)
            % Multiple segments and planes (same count)
            intersector = core.geometry.SegmentByPlaneIntersector(2);
            
            X1 = [0, 1;
                  0, 0];
                  
            X2 = [2, 3;
                  0, 0];
                  
            a = [1, 0;
                 0, 1];
            b = [-1, -0];
            
            [X, TF] = intersector.intersect(X1, X2, a, b);
            
            testCase.verifyEqual(size(X, 2), 1);
            testCase.verifyTrue(TF(1));
            testCase.verifyFalse(TF(2));
        end
        
        function testNearParallelCase(testCase)
            % Line segment almost parallel to plane
            intersector = core.geometry.SegmentByPlaneIntersector(3);
            
            X1 = [0; 0; 0];
            X2 = [0; 0; 1e-6];
            a = [0; 0; 1];
            b = -1;
            
            [X, TF] = intersector.intersect(X1, X2, a, b);
            
            testCase.verifyEqual(size(X, 2), 0);
            testCase.verifyFalse(TF);
        end
        
        function testIntersectionBeyondSegment(testCase)
            % Line would intersect plane, but beyond segment end
            intersector = core.geometry.SegmentByPlaneIntersector(2);
            
            X1 = [0; 0];
            X2 = [1; 1];
            a = [1; 1];
            b = -3;
            
            [X, TF] = intersector.intersect(X1, X2, a, b);
            
            testCase.verifyEqual(size(X, 2), 0);
            testCase.verifyFalse(TF);
        end
        
        function testZeroLengthSegment(testCase)
            % Zero-length line segment (point)
            intersector = core.geometry.SegmentByPlaneIntersector(3);
            
            X1 = [1; 1; 1];
            X2 = [1; 1; 1];
            a = [1; 0; 0];
            b = -1;
            
            [X, TF] = intersector.intersect(X1, X2, a, b);
            
            testCase.verifyEqual(size(X, 2), 1);
            testCase.verifyTrue(TF);
            testCase.verifyEqual(X, [1; 1; 1], 'AbsTol', sqrt(eps));
        end
        
        function testZeroNormalVector(testCase)
            % Plane with zero normal vector
            intersector = core.geometry.SegmentByPlaneIntersector(2);
            
            X1 = [0; 0];
            X2 = [1; 1];
            a = [0; 0];
            b = -1;
            
            testCase.verifyError(@() intersector.intersect(X1, X2, a, b), ...
                'core:geometry:SegmentByPlaneIntersector:InvalidInput');
        end
    end
end