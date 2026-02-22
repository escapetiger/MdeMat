classdef TestBall < matlab.unittest.TestCase
    % TESTBALL Unit tests for the Ball class.
    %
    %   TestBall provides comprehensive test coverage for the
    %   Ball class functionality including constructor validation,
    %   magnitude calculation, and point containment testing.
    %
    % See also:
    %   core.geometry.Ball

    properties (Constant)
        Tolerance = 1e-12
    end

    methods (Test)
        function testConstructor(testCase)
            ball1 = core.geometry.Ball([0, 0, 0], 1);
            testCase.verifyEqual(ball1.NDims, 3);
            testCase.verifyEqual(ball1.Center, [0, 0, 0]);
            testCase.verifyEqual(ball1.Radius, 1);

            ball2 = core.geometry.Ball([1, 2, 3], 2);
            testCase.verifyEqual(ball2.NDims, 3);
            testCase.verifyEqual(ball2.Center, [1, 2, 3]);
            testCase.verifyEqual(ball2.Radius, 2);

            testCase.verifyError(@() core.geometry.Ball([1, 2], -1), ?MException);
        end

        function testUnitBall(testCase)
            ball1D = core.geometry.Ball.unit(1);
            testCase.verifyEqual(ball1D.NDims, 1);
            testCase.verifyEqual(ball1D.Center, [0]);
            testCase.verifyEqual(ball1D.Radius, 1);

            ball2D = core.geometry.Ball.unit(2);
            testCase.verifyEqual(ball2D.NDims, 2);
            testCase.verifyEqual(ball2D.Center, [0, 0]);
            testCase.verifyEqual(ball2D.Radius, 1);

            ball3D = core.geometry.Ball.unit(3);
            testCase.verifyEqual(ball3D.NDims, 3);
            testCase.verifyEqual(ball3D.Center, [0, 0, 0]);
            testCase.verifyEqual(ball3D.Radius, 1);

            testCase.verifyError(@() core.geometry.Ball.unit(-1), ?MException);
            testCase.verifyError(@() core.geometry.Ball.unit(0), ?MException);
        end

        function testMagnitude(testCase)
            ball1D = core.geometry.Ball([0], 2);
            testCase.verifyEqual(ball1D.magnitude(), 4);

            ball2D = core.geometry.Ball([0, 0], 3);
            testCase.verifyEqual(ball2D.magnitude(), pi*3^2, 'AbsTol', testCase.Tolerance);

            ball3D = core.geometry.Ball([0, 0, 0], 2);
            testCase.verifyEqual(ball3D.magnitude(), 4/3*pi*2^3, 'AbsTol', testCase.Tolerance);

            ball4D = core.geometry.Ball([0, 0, 0, 0], 2);
            testCase.verifyEqual(ball4D.magnitude(), pi^2*2^4/2, 'AbsTol', testCase.Tolerance);
        end

        function testIsInside(testCase)
            ball = core.geometry.Ball([1, 1], 2);

            interiorPoints = [1, 0, 2; 1, 0, 1];
            testCase.verifyEqual(ball.isInside(interiorPoints), [true, true, true]);

            boundaryPoints = [3, 1, 3; 1, 3, 1];
            testCase.verifyEqual(ball.isInside(boundaryPoints), [false, false, false]);

            exteriorPoints = [4, -2, 1; 1, 4, -2];
            testCase.verifyEqual(ball.isInside(exteriorPoints), [false, false, false]);

            exactBoundaryPoint = [1; 1 + 2];
            testCase.verifyEqual(ball.isInside(exactBoundaryPoint), false);

            testCase.verifyError(@() ball.isInside([1, 2, 3]), ?MException);
        end

        function testIsOnBoundary(testCase)
            ball = core.geometry.Ball([1, 1], 2);

            interiorPoints = [1, 0, 2; 1, 0, 1];
            testCase.verifyEqual(ball.isOnBoundary(interiorPoints), [false, false, false]);

            boundaryPoints = [3, 1, 3, 1 - sqrt(2); 1, 3, 1, 1 + sqrt(2)];
            testCase.verifyEqual(ball.isOnBoundary(boundaryPoints), [true, true, true, true]);

            exteriorPoints = [4, -2, 1; 1, 4, -2];
            testCase.verifyEqual(ball.isOnBoundary(exteriorPoints), [false, false, false]);

            almostBoundaryPoint = [1; 1 + 2 * (1 - eps * 1e1)];
            testCase.verifyEqual(ball.isOnBoundary(almostBoundaryPoint), true);
        end

        function testGetBoundary(testCase)
            ball = core.geometry.Ball([1, 2, 3], 2);
            boundary = ball.Boundary;

            testCase.verifyEqual(boundary.Center, [1, 2, 3]);
            testCase.verifyEqual(boundary.Radius, 2);
            testCase.verifyEqual(boundary.NSphDims, 2);
            testCase.verifyClass(boundary, ?core.geometry.Sphere);
        end

        function testZeroRadius(testCase)
            ball = core.geometry.Ball([1, 2], 0);

            testCase.verifyEqual(ball.Radius, 0);
            testCase.verifyEqual(ball.Center, [1, 2]);
            testCase.verifyEqual(ball.magnitude(), 0);

            points = [1, 0, 2; 2, 2, 1];
            testCase.verifyEqual(ball.isInside(points), [false, false, false]);

            centerPoint = [1; 2];
            otherPoints = [1, 0, 2; 1, 2, 1];

            testCase.verifyEqual(ball.isOnBoundary(centerPoint), true);
            testCase.verifyEqual(ball.isOnBoundary(otherPoints), [false, false, false]);

            nearCenterPoint = [1 + eps(1); 2 + eps(1)];
            testCase.verifyEqual(ball.isOnBoundary(nearCenterPoint), true);

            boundary = ball.Boundary;
            testCase.verifyEqual(boundary.Radius, 0);
            testCase.verifyEqual(boundary.Center, [1, 2]);
        end

        function testSphericalCoordinates2D(testCase)
            ball = core.geometry.Ball([0, 0], 2);

            cartPoints = {; ...
                [2; 0], ...
                [0; 2], ...
                [-2; 0], ...
                [0; -2], ...
                [1; 0], ...
                [0; 1], ...
                [1; 1]; ...
                };

            expectedSpherical = {; ...
                [2; 0], ...
                [2; pi / 2], ...
                [2; pi], ...
                [2; 3 * pi / 2], ...
                [1; 0], ...
                [1; pi / 2], ...
                [sqrt(2); pi / 4]; ...
                };

            for i = 1:numel(cartPoints)
                cartPoint = cartPoints{i};
                expectedSph = expectedSpherical{i};

                spherical = ball.cartesianToSpherical(cartPoint);

                testCase.verifyEqual(spherical(1), expectedSph(1), 'AbsTol', testCase.Tolerance);

                angleDiff = mod(spherical(2)-expectedSph(2), 2*pi);
                testCase.verifyTrue(angleDiff < testCase.Tolerance || ...
                    abs(angleDiff-2*pi) < testCase.Tolerance);

                reconstructed = ball.sphericalToCartesian(spherical);
                testCase.verifyEqual(reconstructed, cartPoint, 'AbsTol', testCase.Tolerance);
            end
        end

        function testSphericalCoordinates3D(testCase)
            ball = core.geometry.Ball([0, 0, 0], 2);

            cartPoints = {; ...
                [2; 0; 0], ...
                [0; 2; 0], ...
                [0; 0; 2], ...
                [-2; 0; 0], ...
                [0; -2; 0], ...
                [0; 0; -2], ...
                [1; 1; 1], ...
                [0; 0; 0]; ...
                };

            expectedSpherical = {; ...
                [2; 0; 0], ...
                [2; pi / 2; 0], ...
                [2; pi / 2; pi / 2], ...
                [2; pi; 0], ...
                [2; pi / 2; pi], ...
                [2; pi / 2; 3 * pi / 2], ...
                [sqrt(3); acos(1/sqrt(3)); pi / 4], ...
                [0; 0; 0]; ...
                };

            for i = 1:numel(cartPoints)
                cartPoint = cartPoints{i};
                expectedSph = expectedSpherical{i};

                skipAngleCheck = (i == 8);

                spherical = ball.cartesianToSpherical(cartPoint);

                testCase.verifyEqual(spherical(1), expectedSph(1), 'AbsTol', testCase.Tolerance);

                if ~skipAngleCheck
                    thetaDiff = mod(spherical(2)-expectedSph(2), pi);
                    testCase.verifyTrue(thetaDiff < testCase.Tolerance || ...
                        abs(thetaDiff-pi) < testCase.Tolerance);

                    if ~(i == 3 || i == 6)
                        phiDiff = mod(spherical(3)-expectedSph(3), 2*pi);
                        testCase.verifyTrue(phiDiff < testCase.Tolerance || ...
                            abs(phiDiff-2*pi) < testCase.Tolerance);
                    end
                end

                reconstructed = ball.sphericalToCartesian(spherical);
                testCase.verifyEqual(reconstructed, cartPoint, 'AbsTol', testCase.Tolerance);
            end
        end

        function testNonOriginCenteredBall(testCase)
            ball = core.geometry.Ball([1, 2], 2);

            cartPoints = {; ...
                [3; 2], ...
                [1; 4], ...
                [-1; 2], ...
                [1; 0], ...
                [1; 2], ...
                [2; 3]; ...
                };

            expectedSpherical = {; ...
                [2; 0], ...
                [2; pi / 2], ...
                [2; pi], ...
                [2; 3 * pi / 2], ...
                [0; 0], ...
                [sqrt(2); pi / 4]; ...
                };

            for i = 1:numel(cartPoints)
                cartPoint = cartPoints{i};
                expectedSph = expectedSpherical{i};

                skipAngleCheck = (i == 5);

                spherical = ball.cartesianToSpherical(cartPoint);

                testCase.verifyEqual(spherical(1), expectedSph(1), 'AbsTol', testCase.Tolerance);

                if ~skipAngleCheck
                    angleDiff = mod(spherical(2)-expectedSph(2), 2*pi);
                    testCase.verifyTrue(angleDiff < testCase.Tolerance || ...
                        abs(angleDiff-2*pi) < testCase.Tolerance);
                end

                reconstructed = ball.sphericalToCartesian(spherical);
                testCase.verifyEqual(reconstructed, cartPoint, 'AbsTol', testCase.Tolerance);
            end
        end

        function testMultidimensionalInput(testCase)
            ball = core.geometry.Ball([0, 0], 1);

            cartesian = zeros(2, 2, 2);
            cartesian(:, 1, 1) = [1; 0];
            cartesian(:, 2, 1) = [0; 1];
            cartesian(:, 1, 2) = [-1; 0];
            cartesian(:, 2, 2) = [0; -1];

            spherical = ball.cartesianToSpherical(cartesian);

            testCase.verifySize(spherical, size(cartesian));
            testCase.verifyEqual(spherical(1, :, :), ones(1, 2, 2), 'AbsTol', testCase.Tolerance);

            testCase.verifyEqual(spherical(2, 1, 1), 0, 'AbsTol', testCase.Tolerance);
            testCase.verifyEqual(spherical(2, 2, 1), pi/2, 'AbsTol', testCase.Tolerance);
            testCase.verifyEqual(spherical(2, 1, 2), pi, 'AbsTol', testCase.Tolerance);
            testCase.verifyEqual(spherical(2, 2, 2), 3*pi/2, 'AbsTol', testCase.Tolerance);

            reconstructed = ball.sphericalToCartesian(spherical);
            testCase.verifyEqual(reconstructed, cartesian, 'AbsTol', testCase.Tolerance);
        end

        function testOutsidePointRejection(testCase)
            ball = core.geometry.Ball([0, 0], 1);
            outsidePoint = [1.5; 0];

            testCase.verifyError(@() ball.cartesianToSpherical(outsidePoint), ...
                'core:geometry:Ball:InvalidInput');
        end
    end
end