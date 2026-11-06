classdef TestQuadratureIdentities < matlab.unittest.TestCase
    methods (Test)
        function polynomialExactnessOnPhysicalIntervals(testCase)
            a = 2;
            b = 5;
            geometry = mdemat.data.geometry.Orthotope([a, b]);
            for count = 1:7
                rules = {mdemat.common.quadrature.GaussLegendreRule(1), ...
                    mdemat.common.quadrature.GaussLobattoRule(1), ...
                    mdemat.common.quadrature.ClosedNewtonCotesRule(1)};
                exactness = [2*count-1, max(1, 2*count-3), count-mod(count+1, 2)];
                for family = 1:3
                    [points, weights] = rules{family}.generate(geometry, count);
                    centered = (2*points-a-b)/(b-a);
                    for degree = 0:exactness(family)
                        expected = (b-a)*(mod(degree, 2) == 0)/(degree+1);
                        actual = weights*centered(:).^degree;
                        testCase.verifyEqual(actual, expected, 'AbsTol', 2e-11);
                    end
                end
            end
        end

        function ambientVolumeAndFaceMeasureAreDistinct(testCase)
            face = mdemat.data.geometry.Orthotope([2, 2, -1, 3]);
            rules = {mdemat.common.quadrature.GaussLegendreRule(2), ...
                mdemat.common.quadrature.GaussLobattoRule(2), ...
                mdemat.common.quadrature.ClosedNewtonCotesRule(2), ...
                mdemat.common.quadrature.PeriodicTrapezoidalRule(2)};
            for family = 1:numel(rules)
                [~, volumeWeights] = rules{family}.generate(face, [4, 5]);
                testCase.verifyEqual(sum(volumeWeights), 0);
                [points, faceWeights] = rules{family}.generateOnBoundary(face, [4, 5]);
                testCase.verifyEqual(points(1, :), 2*ones(1, size(points, 2)));
                testCase.verifyEqual(sum(faceWeights), 4, 'AbsTol', 2e-12);
            end
        end

        function standardGaussianMoments(testCase)
            rule = mdemat.common.quadrature.GaussHermiteRule(1);
            [points, weights] = rule.generate(mdemat.data.geometry.Orthotope([-Inf, Inf]), 6);
            for degree = 0:11
                if mod(degree, 2) == 0
                    expected = prod(1:2:degree-1);
                else
                    expected = 0;
                end
                testCase.verifyEqual(weights*points(:).^degree, expected, ...
                    'AbsTol', 5e-12, 'RelTol', 5e-12);
            end
        end

        function periodicModesIntegrateToZero(testCase)
            a = 0.3;
            b = 4.1;
            rule = mdemat.common.quadrature.PeriodicTrapezoidalRule(1);
            [points, weights] = rule.generate(mdemat.data.geometry.Orthotope([a, b]), 16);
            for frequency = 1:15
                values = exp(2i*pi*frequency*(points-a)/(b-a));
                testCase.verifyEqual(weights*values(:), 0, 'AbsTol', 2e-13);
            end
        end

        function sphericalPhysicalMomentsInAllSupportedDimensions(testCase)
            radius = 1.7;
            pointSphere = mdemat.data.geometry.Sphere(2, radius);
            [points, weights] = mdemat.common.quadrature.GaussTrapezoidalRule(1).generate( ...
                pointSphere, 2);
            testCase.verifyEqual(points, 2+[-radius, radius], 'AbsTol', 1e-14);
            testCase.verifyEqual(sum(weights), 2);
            for dimension = 2:4
                center = (1:dimension)/3;
                sphere = mdemat.data.geometry.Sphere(center, radius);
                counts = [repmat(6, 1, dimension-2), 16];
                rule = mdemat.common.quadrature.GaussTrapezoidalRule(dimension);
                [points, weights] = rule.generate(sphere, counts);
                area = 2*pi^(dimension/2)/gamma(dimension/2)*radius^(dimension-1);
                centered = points-center(:);
                testCase.verifyEqual(sum(weights), area, 'AbsTol', 3e-12);
                testCase.verifyEqual(centered*weights(:), zeros(dimension, 1), ...
                    'AbsTol', 3e-12);
                testCase.verifyEqual((centered.^2)*weights(:), ...
                    area*radius^2/dimension*ones(dimension, 1), 'AbsTol', 8e-12);
            end
        end

        function monteCarloCosineCoordinatesPreserveTheAzimuth(testCase)
            for dimension = 2:4
                sphere = mdemat.data.geometry.Sphere(ones(1, dimension), 1.3);
                car = mdemat.common.quadrature.MonteCarloSphereRule(dimension, seed=17);
                cosine = mdemat.common.quadrature.MonteCarloSphereRule( ...
                    dimension, seed=17, coord='cossph');
                [cartesian, ~] = car.generate(sphere, 31);
                [angles, ~] = cosine.generate(sphere, 31);
                angles(1:dimension-2, :) = acos(-angles(1:dimension-2, :));
                testCase.verifyEqual(sphere.sphericalToCartesian(angles), cartesian, ...
                    'AbsTol', 5e-13);
            end
        end

        function maxwellianFirstDerivativeMatchesIndependentDifference(testCase)
            temperature = 0.7;
            center = [0.3; -0.2];
            distribution = mdemat.data.function.Maxwellian(2, ...
                temperature=temperature, center=center);
            points = [-0.8, 0.1, 0.9; 0.3, -0.7, 0.4];
            h = eps^(1/3);
            for dimension = 1:2
                shift = zeros(2, 1);
                shift(dimension) = h;
                difference = (distribution.eval(points+shift)- ...
                    distribution.eval(points-shift))/(2*h);
                order = double((1:2) == dimension);
                testCase.verifyEqual(distribution.diff(points, order), difference, ...
                    'AbsTol', 3e-10);
            end
            testCase.verifyError(@() distribution.diff(points, [2, 0]), ...
                'MdeMat:Maxwellian:UnsupportedDerivative');
        end
    end
end
