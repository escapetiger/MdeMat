classdef TestHyperball < matlab.unittest.TestCase
    % TESTHYPERBALL Unit tests for the Hyperball class.

    properties (Constant)
        Tolerance = 1e-12
    end

    methods (Test)
        function testConstructor(testCase)
            ball1 = core.geometry.Hyperball(3);
            testCase.verifyEqual(ball1.nDims, 3);
            testCase.verifyEqual(ball1.centroid, [0, 0, 0]);
            testCase.verifyEqual(ball1.radius, 1);

            ball2 = core.geometry.Hyperball([1, 2, 3], 2);
            testCase.verifyEqual(ball2.nDims, 3);
            testCase.verifyEqual(ball2.centroid, [1, 2, 3]);
            testCase.verifyEqual(ball2.radius, 2);

            testCase.verifyError(@() core.geometry.Hyperball(-1), ?MException);
            testCase.verifyError(@() core.geometry.Hyperball([1, 2], -1), ?MException);
        end

        function testMagnitude(testCase)
            ball1D = core.geometry.Hyperball(1, 2);
            testCase.verifyEqual(ball1D.magnitude(), 4);

            ball2D = core.geometry.Hyperball([0, 0], 3);
            testCase.verifyEqual(ball2D.magnitude(), pi*3^2, 'AbsTol', testCase.Tolerance);

            ball3D = core.geometry.Hyperball([0, 0, 0], 2);
            testCase.verifyEqual(ball3D.magnitude(), 4/3*pi*2^3, 'AbsTol', testCase.Tolerance);

            ball4D = core.geometry.Hyperball([0, 0, 0, 0], 2);
            testCase.verifyEqual(ball4D.magnitude(), pi^2*2^4/2, 'AbsTol', testCase.Tolerance);
        end

        function testIsInside(testCase)
            ball = core.geometry.Hyperball([1, 1], 2);

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
            ball = core.geometry.Hyperball([1, 1], 2);

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
            ball = core.geometry.Hyperball([1, 2, 3], 2);
            boundary = ball.boundary;

            testCase.verifyEqual(boundary.centroid, [1, 2, 3]);
            testCase.verifyEqual(boundary.radius, 2);
            testCase.verifyEqual(boundary.nSphereDims, 2);
            testCase.verifyClass(boundary, ?core.geometry.Hypersphere);
        end

        function testZeroRadius(testCase)
            ball = core.geometry.Hyperball([1, 2], 0);

            testCase.verifyEqual(ball.radius, 0);
            testCase.verifyEqual(ball.centroid, [1, 2]);
            testCase.verifyEqual(ball.magnitude(), 0);

            points = [1, 0, 2; 2, 2, 1];
            testCase.verifyEqual(ball.isInside(points), [false, false, false]);

            centroidPoint = [1; 2];
            otherPoints = [1, 0, 2; 1, 2, 1];

            testCase.verifyEqual(ball.isOnBoundary(centroidPoint), true);
            testCase.verifyEqual(ball.isOnBoundary(otherPoints), [false, false, false]);

            nearCentroidPoint = [1 + eps(1); 2 + eps(1)];
            testCase.verifyEqual(ball.isOnBoundary(nearCentroidPoint), true);

            boundary = ball.boundary();
            testCase.verifyEqual(boundary.radius, 0);
            testCase.verifyEqual(boundary.centroid, [1, 2]);
        end

        function testHypersphericalCoordinates2D(testCase)
            ball = core.geometry.Hyperball([0, 0], 2);

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

        function testHypersphericalCoordinates3D(testCase)
            ball = core.geometry.Hyperball([0, 0, 0], 2);

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
            ball = core.geometry.Hyperball([1, 2], 2);

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
            ball = core.geometry.Hyperball([0, 0], 1);

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
            ball = core.geometry.Hyperball([0, 0], 1);
            outsidePoint = [1.5; 0];

            testCase.verifyError(@() ball.cartesianToSpherical(outsidePoint), ...
                'core:geometry:HypersphericalGeometry:InvalidInput');
        end
    end
end