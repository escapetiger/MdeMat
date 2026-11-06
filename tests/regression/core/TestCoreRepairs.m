classdef TestCoreRepairs < matlab.unittest.TestCase
    % Independent regressions for the retained v4 core/approximation contracts.

    methods (Test)
        function chebyshevValuesAndWeightedOrthonormality(testCase)
            theta = ((1:32)-0.5)*pi/32;
            x = cos(theta);
            expected = cos((0:7)'*theta);
            ordinary = mdemat.common.spectral.ChebyshevBasis(7);
            normalized = mdemat.common.spectral.ChebyshevBasis(7,isNormalized=true);
            testCase.verifyEqual(ordinary.eval(x),expected,AbsTol=2e-14);
            factors = [1/sqrt(pi);sqrt(2/pi)*ones(7,1)];
            values = normalized.eval(x);
            testCase.verifyEqual(values,factors.*expected,AbsTol=2e-14);
            testCase.verifyEqual((pi/32)*values*values.',eye(8),AbsTol=3e-14);
        end
        function affineReconstructionAndInverse(testCase)
            element = TestCoreRepairs.linearElement();
            space = mdemat.common.fem.AffineSpace(element);
            points = [-0.4, 0, 0.3];
            coefficients = [2; 3];
            direction = [0.2, -0.1, 0.4];
            expected = 2 + 3 * points + direction;
            testCase.verifyEqual(space.NDims, 1);
            testCase.verifyEqual(space.eval(points, coefficients, direction), expected, ...
                'AbsTol', 1e-13);
            testCase.verifyEqual(space.direction(expected, points, coefficients), direction, ...
                'AbsTol', 1e-13);

            scaled = mdemat.common.fem.ScaledAffineSpace(element, []);
            testCase.verifyEqual(scaled.eval(points, coefficients, direction), expected, ...
                'AbsTol', 1e-13);
            scaled.setScales([2; 4; 0.5]);
            expectedScaled = 4 + 12 * points + 0.5 * direction;
            testCase.verifyEqual(scaled.eval(points, coefficients, direction), ...
                expectedScaled, 'AbsTol', 1e-13);
            testCase.verifyEqual(scaled.direction(expectedScaled, points, coefficients), ...
                direction, 'AbsTol', 1e-13);
        end

        function explicitPhysicalFittingPoints(testCase)
            element = TestCoreRepairs.linearElement();
            mesh = mdemat.data.mesh.UniformGrid([0, 1], 2);
            space = mdemat.common.fem.FiniteElementSpace(mesh, element);
            points = mesh.collocate(element.Volume.Nodes);
            coefficients = space.fit(@(x) 2 + 3 * x, X=points);
            expected = [2.75; 1.5; 4.25; 1.5];
            testCase.verifyEqual(coefficients, expected, 'AbsTol', 2e-13);
            testCase.verifyEqual(space.fit(@(x) 2 + 3 * x), expected, 'AbsTol', 2e-13);
            testCase.verifyEqual(space.fit(@(x) 2 + 3 * x, ...
                EI=2, X=points(:, 4:6)), expected(3:4), 'AbsTol', 2e-13);
        end

        function denseAndSparseCellDerivatives(testCase)
            first = mdemat.data.function.Linear(2, b=1);
            second = mdemat.data.function.Linear(-3, b=4);
            points = {[-0.4, 0.2], [-0.6, 0.1, 0.5]};
            for pattern = {'', 'P', 'Q'}
                basis = mdemat.common.approximation.SeparableFunction( ...
                    factors={first, second}, pattern=pattern{1});
                actual = basis.diff(points, [1, 0]);
                expected = kron(2 * ones(1, 2), 4 - 3 * points{2});
                testCase.verifyEqual(actual, expected, 'AbsTol', 1e-13);
                testCase.verifyEqual(basis.eval(points), ...
                    kron(1 + 2 * points{1}, 4 - 3 * points{2}), 'AbsTol', 1e-13);
            end
        end

        function physicalIntegrationOnNonuniformMesh(testCase)
            vertices = [-2, -1.7, -0.2, 1, 3];
            mesh = mdemat.data.mesh.Grid({vertices});
            space = mdemat.common.fem.FiniteElementSpace(mesh, TestCoreRepairs.linearElement());
            weights = space.integrationWeights();
            testCase.verifySize(weights, [space.NGlobalDofs, 1]);
            testCase.verifyEqual(weights, kron(diff(vertices).', [1; 0]), ...
                'AbsTol', 2e-13);
            constant = space.fit(@(x) ones(size(x)));
            linear = space.fit(@(x) 2 + 3*x);
            testCase.verifyEqual(weights.'*constant, vertices(end)-vertices(1), ...
                'AbsTol', 2e-13);
            exact = 2*(vertices(end)-vertices(1)) + ...
                1.5*(vertices(end)^2-vertices(1)^2);
            testCase.verifyEqual(weights.'*linear, exact, 'AbsTol', 2e-12);
        end

        function complexProjectionUsesHermitianMetric(testCase)
            basis = mdemat.data.function.Linear([0; 1 + 1i], b=[1; 0]);
            approximation = mdemat.common.approximation.ModalApproximator(basis);
            points = [-1, 0, 1];
            weights = [1, 4, 1] / 3;
            approximation.setMass(points, weights);
            approximation.setDesign(points, weights);
            testCase.verifyEqual(approximation.Mass, diag([2, 4/3]), 'AbsTol', 1e-13);
            coefficients = [2 - 1i; 3 + 2i];
            values = coefficients.' * basis.eval(points);
            load = approximation.embed(values, basis.eval(points), weights);
            testCase.verifyEqual(approximation.project(load), coefficients, 'AbsTol', 1e-13);
        end

        function rectangularFitAvoidsNormalEquations(testCase)
            approximation = mdemat.common.approximation.NodalApproximator( ...
                mdemat.data.function.Linear([0; 1], b=[1; 0]));
            points = linspace(-1, 1, 17);
            approximation.setDesign(points);
            coefficients = [2 - 1i; 3 + 2i];
            values = approximation.Design * coefficients;
            testCase.verifyEqual(approximation.fit(values), coefficients, 'AbsTol', 1e-13);
        end

        function multilevelCopyAndVectorExtrapolation(testCase)
            rule = mdemat.common.quadrature.ClosedNewtonCotesRule(1);
            integrator = mdemat.common.quadrature.MultiLevelIntegrator(rule, nLevels=2);
            integrator.setPoints(mdemat.data.geometry.Orthotope([0, 1]), 2);
            integrator.setValues(@(x) [x; x.^2]);
            copy = integrator.copy();
            copy.Integrators(1).Nodes(1) = 0.123;
            testCase.verifyNotEqual(copy.Integrators(1).Nodes, integrator.Integrators(1).Nodes);
            testCase.verifyEqual(integrator.integrate(useExtrapolation=true), [1/2; 1/3], ...
                'AbsTol', 1e-13);
            testCase.verifyEqual(integrator.integrate(@(x) [x; x.^2], ...
                useExtrapolation=true), [1/2; 1/3], 'AbsTol', 1e-13);
            integrator.setPoints(mdemat.data.geometry.Orthotope([0, 2]), 3);
            testCase.verifyEmpty(integrator.Integrators(1).Values);
            testCase.verifyEmpty(integrator.Integrators(2).Values);
        end

        function seededRulesPreserveGlobalState(testCase)
            state = rng;
            cleanup = onCleanup(@() rng(state)); %#ok<NASGU>
            rng(19);
            expectedState = rng;
            orthotope = mdemat.common.quadrature.MonteCarloOrthotopeRule(1, seed=7);
            [points, weights] = orthotope.generate(mdemat.data.geometry.Orthotope([2, 5]), 20);
            testCase.verifyEqual(rng, expectedState);
            testCase.verifyEqual(sum(weights), 3, 'AbsTol', 1e-14);
            [again, ~] = orthotope.generate(mdemat.data.geometry.Orthotope([2, 5]), 20);
            testCase.verifyEqual(points, again);
            sphere = mdemat.common.quadrature.MonteCarloSphereRule(3, seed=11);
            [points, weights] = sphere.generate(mdemat.data.geometry.Sphere([0, 0, 0], 2), 200);
            testCase.verifyEqual(rng, expectedState);
            testCase.verifyEqual(vecnorm(points, 2, 1), 2 * ones(1, 200), 'AbsTol', 1e-13);
            testCase.verifyEqual(sum(weights), 16*pi, 'AbsTol', 1e-11);
        end

        function angularRuleRestoresCallerPath(testCase)
            oldPath = path;
            rule = mdemat.common.quadrature.LebedevRule();
            [points, weights] = rule.generate(mdemat.data.geometry.Sphere([0, 0, 0], 1), 6);
            testCase.verifyEqual(path, oldPath);
            testCase.verifyEqual(weights * (points(1, :).^2).', 4*pi/3, 'AbsTol', 1e-13);
        end

        function compiledFamiliesResolveOutsideRepository(testCase)
            originalDirectory = pwd;
            cleanup = onCleanup(@() cd(originalDirectory)); %#ok<NASGU>
            cd(tempdir);
            basis = mdemat.common.spectral.CompiledOrthogonalBasis(1, 'normal_spherical_harmonic');
            values = basis.eval([0, 0.5; 0, 1]);
            testCase.verifyEqual(values, repmat(1/sqrt(4*pi), 1, 2), 'AbsTol', 1e-13);
            legendre = mdemat.common.spectral.CompiledOrthogonalBasis(1, 'normal_legendre');
            testCase.verifyEqual(legendre.eval([-0.2, 0.3]), repmat(1/sqrt(2), 1, 2), ...
                'AbsTol', 1e-13);
        end
    end

    methods (Static, Access = private)
        function element = linearElement()
            geometry = mdemat.data.geometry.Orthotope([-0.5, 0.5]);
            volume = mdemat.common.quadrature.Integrator(mdemat.common.quadrature.GaussLegendreRule(1));
            volume.setPoints(geometry, 3);
            basis = mdemat.data.function.Linear([0; 1], b=[1; 0]);
            approximation = mdemat.common.approximation.ModalApproximator(basis);
            approximation.setMass(volume.Nodes, volume.Weights);
            approximation.setDesign(volume.Nodes, volume.Weights);
            element = mdemat.common.fem.L2Element(geometry, approximation, volume);
        end
    end
end
