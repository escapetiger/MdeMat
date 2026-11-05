classdef TestHypersphere < matlab.unittest.TestCase
    % TESTHYPERSPHERE Unit tests for the Hypersphere class.

    properties (Constant)
        Tolerance = 1e-12
    end

    methods (Test)
        function testConstructor(testCase)
            sphere1 = core.geometry.Hypersphere(2);
            testCase.verifyEqual(sphere1.nDims, 3);
            testCase.verifyEqual(sphere1.nSphereDims, 2);
            testCase.verifyEqual(sphere1.centroid, [0, 0, 0]);
            testCase.verifyEqual(sphere1.radius, 1);

            sphere2 = core.geometry.Hypersphere([1, 2, 3], 2);
            testCase.verifyEqual(sphere2.nDims, 3);
            testCase.verifyEqual(sphere2.nSphereDims, 2);
            testCase.verifyEqual(sphere2.centroid, [1, 2, 3]);
            testCase.verifyEqual(sphere2.radius, 2);

            testCase.verifyError(@() core.geometry.Hypersphere(-1), ?MException);
            testCase.verifyError(@() core.geometry.Hypersphere([1, 2], -1), ?MException);
        end

        function testMagnitude(testCase)
            sphere0D = core.geometry.Hypersphere([0], 1);
            testCase.verifyEqual(sphere0D.magnitude(), 2);

            sphere1D = core.geometry.Hypersphere([0, 0], 3);
            testCase.verifyEqual(sphere1D.magnitude(), 2*pi*3, 'AbsTol', testCase.Tolerance);

            sphere2D = core.geometry.Hypersphere([0, 0, 0], 2);
            testCase.verifyEqual(sphere2D.magnitude(), 4*pi*2^2, 'AbsTol', testCase.Tolerance);

            sphere3D = core.geometry.Hypersphere([0, 0, 0, 0], 2);
            testCase.verifyEqual(sphere3D.magnitude(), 2*pi^2*2^3, 'AbsTol', testCase.Tolerance);
        end

        function testIsInside(testCase)
            sphere = core.geometry.Hypersphere([1, 1, 1], 2);

            interiorPoints = [1, 0, 2; 1, 0, 1; 2, 2, 0];
            testCase.verifyEqual(sphere.isInside(interiorPoints), [false, false, false]);

            testCase.verifyError(@() sphere.isInside([1, 2]), ?MException);
        end

        function testIsOnBoundary(testCase)
            sphere = core.geometry.Hypersphere([1, 1, 1], 2);

            boundaryPoints = [3, 1, 1, 1 + sqrt(8) / 2; 1, 3, 1, 1 - sqrt(8) / 2; 1, 1, 3, 1];
            testCase.verifyEqual(sphere.isOnBoundary(boundaryPoints), [true, true, true, true]);

            nonBoundaryPoints = [4, 0, 1; 1, 0, 1; 1, 1, 1];
            testCase.verifyEqual(sphere.isOnBoundary(nonBoundaryPoints), [false, false, false]);

            almostBoundaryPoint = [1; 1; 1 + 2 * (1 - eps * 1e1)];
            testCase.verifyEqual(sphere.isOnBoundary(almostBoundaryPoint), true);
        end

        function testZeroRadius(testCase)
            sphere = core.geometry.Hypersphere([1, 2, 3], 0);

            testCase.verifyEqual(sphere.radius, 0);
            testCase.verifyEqual(sphere.centroid, [1, 2, 3]);
            testCase.verifyEqual(sphere.magnitude(), 0);

            points = [1, 0, 2; 2, 2, 1; 3, 3, 3];
            testCase.verifyEqual(sphere.isInside(points), [false, false, false]);

            centroidPoint = [1; 2; 3];
            otherPoints = [1, 0, 2; 1, 2, 1; 4, 5, 6];

            testCase.verifyEqual(sphere.isOnBoundary(centroidPoint), true);
            testCase.verifyEqual(sphere.isOnBoundary(otherPoints), [false, false, false]);

            nearCentroidPoint = [1 + eps(1); 2 + eps(1); 3 + eps(1)];
            testCase.verifyEqual(sphere.isOnBoundary(nearCentroidPoint), true);
        end

        function testHypersphericalCoordinates2D(testCase)
            sphere = core.geometry.Hypersphere([0, 0], 2);

            cartPoints = {; ...
                [2; 0], ...
                [0; 2], ...
                [-2; 0], ...
                [0; -2]; ...
                };

            expectedSpherical = {; ...
                [0], ...
                [pi / 2], ...
                [pi], ...
                [3 * pi / 2]; ...
                };

            for i = 1:numel(cartPoints)
                cartPoint = cartPoints{i};
                expectedSph = expectedSpherical{i};

                spherical = sphere.cartesianToSpherical(cartPoint);

                testCase.verifySize(spherical, [1, 1]);

                angleDiff = mod(spherical(1)-expectedSph(1), 2*pi);
                testCase.verifyTrue(angleDiff < testCase.Tolerance || ...
                    abs(angleDiff-2*pi) < testCase.Tolerance);

                reconstructed = sphere.sphericalToCartesian(spherical);
                testCase.verifyEqual(reconstructed, cartPoint, 'AbsTol', testCase.Tolerance);

                testCase.verifyTrue(sphere.isOnBoundary(reconstructed));
            end
        end

        function testHypersphericalCoordinates3D(testCase)
            sphere = core.geometry.Hypersphere([0, 0, 0], 2);

            cartPoints = {; ...
                [2; 0; 0], ...
                [0; 2; 0], ...
                [0; 0; 2], ...
                [-2; 0; 0], ...
                [0; -2; 0], ...
                [0; 0; -2]; ...
                };

            expectedSpherical = {; ...
                [0; 0], ...
                [pi / 2; 0], ...
                [pi / 2; pi / 2], ...
                [pi; 0], ...
                [pi / 2; pi], ...
                [pi / 2; 3 * pi / 2]; ...
                };

            for i = 1:numel(cartPoints)
                cartPoint = cartPoints{i};
                expectedSph = expectedSpherical{i};

                skipPhiCheck = (i == 3 || i == 6);

                spherical = sphere.cartesianToSpherical(cartPoint);

                testCase.verifySize(spherical, [2, 1]);

                thetaDiff = mod(spherical(1)-expectedSph(1), pi);
                testCase.verifyTrue(thetaDiff < testCase.Tolerance || ...
                    abs(thetaDiff-pi) < testCase.Tolerance);

                if ~skipPhiCheck
                    phiDiff = mod(spherical(2)-expectedSph(2), 2*pi);
                    testCase.verifyTrue(phiDiff < testCase.Tolerance || ...
                        abs(phiDiff-2*pi) < testCase.Tolerance);
                end

                reconstructed = sphere.sphericalToCartesian(spherical);
                testCase.verifyEqual(reconstructed, cartPoint, 'AbsTol', testCase.Tolerance);

                testCase.verifyTrue(sphere.isOnBoundary(reconstructed));
            end
        end

        function testNonOriginCenteredSphere(testCase)
            sphere = core.geometry.Hypersphere([1, 2], 2);

            cartPoints = {; ...
                [3; 2], ...
                [1; 4], ...
                [-1; 2], ...
                [1; 0]; ...
                };

            expectedSpherical = {; ...
                [0], ...
                [pi / 2], ...
                [pi], ...
                [3 * pi / 2]; ...
                };

            for i = 1:numel(cartPoints)
                cartPoint = cartPoints{i};
                expectedSph = expectedSpherical{i};

                spherical = sphere.cartesianToSpherical(cartPoint);

                testCase.verifySize(spherical, [1, 1]);

                angleDiff = mod(spherical(1)-expectedSph(1), 2*pi);
                testCase.verifyTrue(angleDiff < testCase.Tolerance || ...
                    abs(angleDiff-2*pi) < testCase.Tolerance);

                reconstructed = sphere.sphericalToCartesian(spherical);
                testCase.verifyEqual(reconstructed, cartPoint, 'AbsTol', testCase.Tolerance);

                testCase.verifyTrue(sphere.isOnBoundary(reconstructed));
            end
        end

        function testMultidimensionalInput(testCase)
            sphere = core.geometry.Hypersphere([0, 0], 1);

            cartesian = zeros(2, 2, 2);
            cartesian(:, 1, 1) = [1; 0];
            cartesian(:, 2, 1) = [0; 1];
            cartesian(:, 1, 2) = [-1; 0];
            cartesian(:, 2, 2) = [0; -1];

            spherical = sphere.cartesianToSpherical(cartesian);

            expectedSize = [1, 2, 2];
            testCase.verifySize(spherical, expectedSize);

            testCase.verifyEqual(spherical(1, 1, 1), 0, 'AbsTol', testCase.Tolerance);
            testCase.verifyEqual(spherical(1, 2, 1), pi/2, 'AbsTol', testCase.Tolerance);
            testCase.verifyEqual(spherical(1, 1, 2), pi, 'AbsTol', testCase.Tolerance);
            testCase.verifyEqual(spherical(1, 2, 2), 3*pi/2, 'AbsTol', testCase.Tolerance);

            reconstructed = sphere.sphericalToCartesian(spherical);
            testCase.verifyEqual(reconstructed, cartesian, 'AbsTol', testCase.Tolerance);

            for i = 1:2
                for j = 1:2
                    point = reconstructed(:, i, j);
                    testCase.verifyTrue(sphere.isOnBoundary(point));
                end
            end
        end

        function testPointsNotOnSphere(testCase)
            sphere = core.geometry.Hypersphere([0, 0, 0], 1);

            nearPoint = [0.8; 0; 0];
            testCase.verifyError(@() sphere.cartesianToSpherical(nearPoint), ...
                'core:geometry:Hypersphere:InvalidInput');

            outsidePoint = [1.5; 0; 0];
            testCase.verifyError(@() sphere.cartesianToSpherical(outsidePoint), ...
                'core:geometry:Hypersphere:InvalidInput');
        end

        function testTransformValidation(testCase)
            sphere = core.geometry.Hypersphere([0, 0, 0], 2);

            spherical = [pi / 4; pi / 4];
            cartesian = sphere.sphericalToCartesian(spherical);

            testCase.verifyTrue(sphere.isOnBoundary(cartesian));

            distance = norm(cartesian);
            testCase.verifyEqual(distance, sphere.radius, 'AbsTol', testCase.Tolerance);
        end
    end
end