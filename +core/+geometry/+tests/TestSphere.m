classdef TestSphere < matlab.unittest.TestCase
    % TESTSPHERE Unit tests for the Sphere class.

    properties (Constant)
        Tolerance = 1e-12
    end

    methods (Test)
        function testConstructor(testCase)
            sphere1 = core.geometry.Sphere([0, 0, 0], 1);
            testCase.verifyEqual(sphere1.NDims, 3);
            testCase.verifyEqual(sphere1.NSphDims, 2);
            testCase.verifyEqual(sphere1.Center, [0, 0, 0]);
            testCase.verifyEqual(sphere1.Radius, 1);

            sphere2 = core.geometry.Sphere([1, 2, 3], 2);
            testCase.verifyEqual(sphere2.NDims, 3);
            testCase.verifyEqual(sphere2.NSphDims, 2);
            testCase.verifyEqual(sphere2.Center, [1, 2, 3]);
            testCase.verifyEqual(sphere2.Radius, 2);

            testCase.verifyError(@() core.geometry.Sphere([1, 2], -1), ?MException);
        end

        function testUnitSphere(testCase)
            sphere1D = core.geometry.Sphere.unit(1);
            testCase.verifyEqual(sphere1D.NDims, 1);
            testCase.verifyEqual(sphere1D.NSphDims, 0);
            testCase.verifyEqual(sphere1D.Center, [0]);
            testCase.verifyEqual(sphere1D.Radius, 1);

            sphere2D = core.geometry.Sphere.unit(2);
            testCase.verifyEqual(sphere2D.NDims, 2);
            testCase.verifyEqual(sphere2D.NSphDims, 1);
            testCase.verifyEqual(sphere2D.Center, [0, 0]);
            testCase.verifyEqual(sphere2D.Radius, 1);

            sphere3D = core.geometry.Sphere.unit(3);
            testCase.verifyEqual(sphere3D.NDims, 3);
            testCase.verifyEqual(sphere3D.NSphDims, 2);
            testCase.verifyEqual(sphere3D.Center, [0, 0, 0]);
            testCase.verifyEqual(sphere3D.Radius, 1);

            testCase.verifyError(@() core.geometry.Sphere.unit(-1), ?MException);
            testCase.verifyError(@() core.geometry.Sphere.unit(0), ?MException);
        end

        function testMagnitude(testCase)
            sphere0D = core.geometry.Sphere([0], 1);
            testCase.verifyEqual(sphere0D.magnitude(), 2);

            sphere1D = core.geometry.Sphere([0, 0], 3);
            testCase.verifyEqual(sphere1D.magnitude(), 2*pi*3, 'AbsTol', testCase.Tolerance);

            sphere2D = core.geometry.Sphere([0, 0, 0], 2);
            testCase.verifyEqual(sphere2D.magnitude(), 4*pi*2^2, 'AbsTol', testCase.Tolerance);

            sphere3D = core.geometry.Sphere([0, 0, 0, 0], 2);
            testCase.verifyEqual(sphere3D.magnitude(), 2*pi^2*2^3, 'AbsTol', testCase.Tolerance);
        end

        function testIsInside(testCase)
            sphere = core.geometry.Sphere([1, 1, 1], 2);

            interiorPoints = [1, 0, 2; 1, 0, 1; 2, 2, 0];
            testCase.verifyEqual(sphere.isInside(interiorPoints), [false, false, false]);

            testCase.verifyError(@() sphere.isInside([1, 2]), ?MException);
        end

        function testIsOnBoundary(testCase)
            sphere = core.geometry.Sphere([1, 1, 1], 2);

            boundaryPoints = [3, 1, 1, 1 + sqrt(8) / 2; 1, 3, 1, 1 - sqrt(8) / 2; 1, 1, 3, 1];
            testCase.verifyEqual(sphere.isOnBoundary(boundaryPoints), [true, true, true, true]);

            nonBoundaryPoints = [4, 0, 1; 1, 0, 1; 1, 1, 1];
            testCase.verifyEqual(sphere.isOnBoundary(nonBoundaryPoints), [false, false, false]);

            almostBoundaryPoint = [1; 1; 1 + 2 * (1 - eps * 1e1)];
            testCase.verifyEqual(sphere.isOnBoundary(almostBoundaryPoint), true);
        end

        function testZeroRadius(testCase)
            sphere = core.geometry.Sphere([1, 2, 3], 0);

            testCase.verifyEqual(sphere.Radius, 0);
            testCase.verifyEqual(sphere.Center, [1, 2, 3]);
            testCase.verifyEqual(sphere.magnitude(), 0);

            points = [1, 0, 2; 2, 2, 1; 3, 3, 3];
            testCase.verifyEqual(sphere.isInside(points), [false, false, false]);

            centerPoint = [1; 2; 3];
            otherPoints = [1, 0, 2; 1, 2, 1; 4, 5, 6];

            testCase.verifyEqual(sphere.isOnBoundary(centerPoint), true);
            testCase.verifyEqual(sphere.isOnBoundary(otherPoints), [false, false, false]);

            nearCenterPoint = [1 + eps(1); 2 + eps(1); 3 + eps(1)];
            testCase.verifyEqual(sphere.isOnBoundary(nearCenterPoint), true);
        end

        function testSphericalCoordinates2D(testCase)
            sphere = core.geometry.Sphere([0, 0], 2);

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

        function testSphericalCoordinates3D(testCase)
            sphere = core.geometry.Sphere([0, 0, 0], 2);

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
            sphere = core.geometry.Sphere([1, 2], 2);

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
            sphere = core.geometry.Sphere([0, 0], 1);

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
            sphere = core.geometry.Sphere([0, 0, 0], 1);

            nearPoint = [0.8; 0; 0];
            testCase.verifyError(@() sphere.cartesianToSpherical(nearPoint), ...
                'core:geometry:Sphere:InvalidInput');

            outsidePoint = [1.5; 0; 0];
            testCase.verifyError(@() sphere.cartesianToSpherical(outsidePoint), ...
                'core:geometry:Sphere:InvalidInput');
        end

        function testTransformValidation(testCase)
            sphere = core.geometry.Sphere([0, 0, 0], 2);

            spherical = [pi / 4; pi / 4];
            cartesian = sphere.sphericalToCartesian(spherical);

            testCase.verifyTrue(sphere.isOnBoundary(cartesian));

            distance = norm(cartesian);
            testCase.verifyEqual(distance, sphere.Radius, 'AbsTol', testCase.Tolerance);
        end
    end
end