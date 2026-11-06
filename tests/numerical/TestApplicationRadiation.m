classdef TestApplicationRadiation < matlab.unittest.TestCase
    methods (Test)
        function uniformEquilibriumUsesEveryVelocityRepresentation(testCase)
            figures = numel(findall(groot, 'Type', 'figure'));
            configurations = {1, 'macro'; 1, 'micro'; 1, 'combined'; ...
                2, 'combined'; 3, 'combined'};
            for row = 1:size(configurations, 1)
                dimension = configurations{row, 1};
                mode = configurations{row, 2};
                space = TestApplicationRadiation.spatial('nodal', 2, 3);
                velocity = TestApplicationRadiation.velocity(dimension, mode);
                config = struct('scattering', 0.7, 'decomposition', 'full', ...
                    'tOdeIntName', 'sdirk2');
                record = mdemat.application.radiation.problem(space, velocity, 1, ...
                    0.02, 0.01, config=config);
                [result, scheme] = mdemat.application.radiation.solve(record, ...
                    observe=@(state, time) struct('time', time, 'mass', state.History.mass(end)), ...
                    callback=@(event, state) testCase.verifyEqual( ...
                        event.time, state.History.time(end), 'AbsTol', 1e-14));
                if velocity.NRhsDofs > 0
                    element = velocity.Rhs.Element;
                else
                    element = velocity.Lhs.Element;
                end
                points = element.Volume.Nodes;
                actual = result.state.distribution([-0.3, 0.2], points);
                testCase.verifyEqual(actual, ones(size(actual)), 'AbsTol', 2e-10);
                testCase.verifyEqual(result.diagnostics.mass(end), ...
                    element.Geometry.magnitude(), 'AbsTol', 3e-10);
                testCase.verifyEmpty(scheme.Observers.visualizer);
                testCase.verifyEqual(scheme.Observation.Analysis.PointsPerDimension, 0);
                testCase.verifyEmpty(scheme.Observers.timer);
                testCase.verifyEqual(result.times(end), 0.02, 'AbsTol', 1e-14);
                testCase.verifyNumElements(result.observations, numel(result.times));
            end
            testCase.verifyEqual(numel(findall(groot, 'Type', 'figure')), figures);
        end

        function homogeneousScatteringAbsorptionAndSourceMatchExactSolution(testCase)
            space = TestApplicationRadiation.spatial('modal', 2, 3);
            velocity = TestApplicationRadiation.velocity(1, 'combined');
            scattering = 0.7;
            absorption = 0.2;
            source = 0.3;
            initial = @(x, v, t) 1+0.2*v(1, :);
            config = struct('scattering', scattering, 'absorption', absorption, ...
                'source', source, 'decomposition', 'partial', 'tOdeIntName', 'sdirk2');
            final = 0.1;
            record = mdemat.application.radiation.problem(space, velocity, initial, ...
                final, 0.002, config=config);
            result = mdemat.application.radiation.solve(record);
            velocities = velocity.Rhs.Element.Volume.Nodes;
            mean = source/absorption+(1-source/absorption)*exp(-absorption*final);
            exact = mean+0.2*velocities*exp(-(scattering+absorption)*final);
            actual = result.state.distribution(0, velocities);
            testCase.verifyEqual(actual, repmat(exact, space.NMeshElements, 1), ...
                'AbsTol', 8e-7);
        end

        function diffusionScalingHasAnIndependentTelegraphReference(testCase)
            mesh = mdemat.data.mesh.UniformGrid([0, 2*pi], 12);
            element = mdemat.common.fem.BH1OrthotopeElement.modal(1, 3);
            space = mdemat.common.fem.FiniteElementSpace(mesh, element);
            velocity = TestApplicationRadiation.velocity(1, 'combined');
            amplitude = 0.1;
            final = 0.05;
            reference = linspace(-0.5, 0.5, 5);
            positions = mesh.collocate(reference).';
            for epsilon = [0.1, 0.01]
                % Slow eigenmode of rho_t+j_x/epsilon=0,
                % j_t+rho_x/epsilon=-j/epsilon^2, wave number one.
                lambda = -2/(1+sqrt(1-4*epsilon^2));
                initial = @(x, v, t) 1+amplitude*cos(x) - ...
                    epsilon*lambda*amplitude*v(1, :).*sin(x);
                config = struct('epsilon', epsilon, 'timeScale', 1, ...
                    'scatteringScale', -1, 'absorptionScale', 1, ...
                    'scattering', 1, 'decomposition', 'full', 'tOdeIntName', 'sdirk2');
                record = mdemat.application.radiation.problem(space, velocity, initial, ...
                    final, 0.005, config=config);
                result = mdemat.application.radiation.solve(record);
                density = result.state.density(reference);
                exact = 2*(1+amplitude*exp(lambda*final)*cos(positions));
                heatLimit = 2*(1+amplitude*exp(-final)*cos(positions));
                testCase.verifyLessThan(max(abs(density-exact)), 3e-3);
                testCase.verifyLessThan(max(abs(density-heatLimit)), 3e-3+epsilon^2);
                testCase.verifyEqual(result.diagnostics.mass(end), 4*pi, 'AbsTol', 3e-10);
            end
        end

        function pureMicroFreeStreamingAndTimeDependentInflow(testCase)
            velocity = TestApplicationRadiation.velocity(1, 'micro');
            velocities = velocity.Rhs.Element.Volume.Nodes;
            final = 0.03;
            reference = linspace(-0.5, 0.5, 5);
            for boundary = {'periodic', 'inflow'}
                if strcmp(boundary{1}, 'periodic')
                    bounds = [0, 2*pi];
                    initial = @(x, v, t) 1+0.1*cos(x);
                    exact = @(x, v, t) 1+0.1*cos(x-v*t);
                    config = struct('tOdeIntName', 'sdirk2');
                    tolerance = 2e-3; % P1 spatial error on h=2*pi/16.
                else
                    bounds = [0, 1];
                    initial = @(x, v, t) 1+0.2*x;
                    exact = @(x, v, t) 1+0.2*(x-v*t);
                    config = struct('tOdeIntName', 'sdirk2', 'bc', exact);
                    tolerance = 2e-10; % Linear space/time solution is representable.
                end
                mesh = mdemat.data.mesh.UniformGrid(bounds, 16);
                element = mdemat.common.fem.BH1OrthotopeElement.modal(1, 2);
                space = mdemat.common.fem.FiniteElementSpace(mesh, element);
                record = mdemat.application.radiation.problem(space, velocity, initial, ...
                    final, 0.01, config=config);
                result = mdemat.application.radiation.solve(record);
                positions = mesh.collocate(reference).';
                values = result.state.distribution(reference, velocities);
                testCase.verifyEqual(values, exact(positions, velocities, final), ...
                    'AbsTol', tolerance);
            end
        end

        function uniformEquilibriumInTwoAndThreeSpatialDimensions(testCase)
            for dimension = 2:3
                mesh = mdemat.data.mesh.UniformGrid( ...
                    repmat([0, 1], 1, dimension), repmat(2, 1, dimension));
                element = mdemat.common.fem.BH1OrthotopeElement.nodal(dimension, 2);
                space = mdemat.common.fem.FiniteElementSpace(mesh, element);
                velocity = TestApplicationRadiation.velocity(dimension, 'combined');
                config = struct('scattering', 0.5, 'decomposition', 'partial', ...
                    'xBasisOrder', 2, 'tOdeIntName', 'sdirk2');
                record = mdemat.application.radiation.problem(space, velocity, 1, ...
                    0.02, 0.01, config=config);
                result = mdemat.application.radiation.solve(record);
                nodes = velocity.Rhs.Element.Volume.Nodes;
                values = result.state.distribution(zeros(dimension, 1), nodes);
                testCase.verifyEqual(values, ones(size(values)), 'AbsTol', 3e-9);
                testCase.verifyEqual(result.diagnostics.mass(end), ...
                    velocity.Rhs.Element.Geometry.magnitude(), 'AbsTol', 3e-9);
                testCase.verifyEqual(space.NMeshElements, 2^dimension);
            end
        end

        function angularFiltersDampHigherModesAndPreserveMass(testCase)
            space = TestApplicationRadiation.spatial('modal', 2, 2);
            macro = mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2SphereElement.modal(1, 4, reduction='symmetry'));
            velocity = mdemat.common.fem.SumSpace(macro, []);
            initial = @(x, v, t) 1+0.1*(5*v.^3-3*v)/2;
            for filter = {'rot', 'exp'}
                config = struct('useFilter', true, 'filterType', filter{1}, ...
                    'filterStrength', 0.5, 'filterCutoff', 0, 'filterCollisionRate', 0);
                record = mdemat.application.radiation.problem(space, velocity, initial, ...
                    0.1, 0.02, config=config);
                result = mdemat.application.radiation.solve(record);
                if strcmp(filter{1}, 'rot')
                    damping = 0.5^result.execution.steps;
                else
                    damping = 0.5^0.1;
                end
                highest = space.eval(0, result.state.Dofs.U(:, 4));
                testCase.verifyEqual(highest, (0.2/sqrt(7))*damping*ones(2, 1), ...
                    'AbsTol', 2e-12);
                testCase.verifyEqual(result.diagnostics.mass(end), 2, 'AbsTol', 3e-12);
            end
        end

        function densityLimiterUsesPhysicalSamplesAndPreservesAdmissibleMeans(testCase)
            for basis = {'modal', 'nodal'}
                space = TestApplicationRadiation.spatial(basis{1}, 3, 2);
                macro = mdemat.common.fem.SpectralSpace( ...
                    mdemat.common.fem.L2SphereElement.modal(1, 1, reduction='symmetry'));
                velocity = mdemat.common.fem.SumSpace(macro, []);
                initial = @(x, v, t) 1+4*(mod(x, 0.5)/0.5-0.5);
                config = struct('usePositivityLimiter', true, ...
                    'positivityLimiterType', 'zhang_shu');
                record = mdemat.application.radiation.problem(space, velocity, initial, ...
                    0.001, 0.001, config=config);
                result = mdemat.application.radiation.solve(record);
                density = result.state.density(space.Element.Volume.Nodes);
                testCase.verifyGreaterThanOrEqual(min(density), -3e-13);
                testCase.verifyEqual(result.diagnostics.mass(end), 2, 'AbsTol', 3e-12);
                config.positivityLimiterType = 'cutoff';
                record = mdemat.application.radiation.problem(space, velocity, -0.2, ...
                    0.001, 0.001, config=config);
                result = mdemat.application.radiation.solve(record);
                testCase.verifyEqual(result.diagnostics.mass(end), 0, 'AbsTol', 3e-13);
            end
        end

        function problemRecordsRejectUnknownMissingAndIncompatibleSchema(testCase)
            space = TestApplicationRadiation.spatial('nodal', 2, 2);
            velocity = TestApplicationRadiation.velocity(1, 'combined');
            record = mdemat.application.radiation.problem(space, velocity, 1, 0, 0.01);
            extra = record;
            extra.unknown = 1;
            wrongModel = record;
            wrongModel.model = 'other';
            wrongSchema = record;
            wrongSchema.schemaVersion = 2;
            for invalid = {extra, rmfield(record, 'initial'), wrongModel, wrongSchema}
                testCase.verifyError(@() mdemat.application.radiation.solve(invalid{1}), ...
                    'MdeMat:Radiation:InvalidProblem');
            end
        end

        function angularRepresentationBoundsPreserveReducedAndSphericalSpaces(testCase)
            spatial = TestApplicationRadiation.spatial('nodal', 2, 2);
            wrong = mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2OrthotopeElement.hermite(1, 3));
            velocity = mdemat.common.fem.SumSpace([], wrong);
            testCase.verifyError(@() mdemat.application.radiation.problem( ...
                spatial, velocity, 1, 0, 0.01), 'MdeMat:Radiation:UnsupportedRepresentation');

            nodal = mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2SphereElement.nodal(1, 2));
            testCase.verifyError(@() mdemat.application.radiation.problem( ...
                spatial, mdemat.common.fem.SumSpace(nodal, []), 1, 0, 0.01), ...
                'MdeMat:Radiation:UnsupportedRepresentation');
            macro = mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2SphereElement.modal(1, 1));
            otherRadius = mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2SphereElement.nodal(1, 2, radius=2));
            testCase.verifyError(@() mdemat.application.radiation.problem( ...
                spatial, mdemat.common.fem.SumSpace(macro, otherRadius), 1, 0, 0.01), ...
                'MdeMat:Radiation:IncompatibleVelocitySpaces');

            % The interval is the supported symmetry reduction in one dimension.
            macro = mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2SphereElement.modal(1, 1, reduction='symmetry'));
            micro = mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2SphereElement.nodal(1, 3, reduction='symmetry'));
            record = mdemat.application.radiation.problem(spatial, ...
                mdemat.common.fem.SumSpace(macro, micro), 0.7, 0, 0.01);
            result = mdemat.application.radiation.solve(record);
            testCase.verifyEqual(result.state.Dofs.F, ...
                0.7*ones(size(result.state.Dofs.F)), 'AbsTol', 2e-12);

            % In the original 2D symmetry scenarios both angular domains are S^2,
            % even though the two factory calls use physical dimensions 2 and 3.
            mesh = mdemat.data.mesh.UniformGrid([0, 1, -1, 1], [1, 1]);
            spatial2 = mdemat.common.fem.FiniteElementSpace(mesh, ...
                mdemat.common.fem.BH1OrthotopeElement.nodal(2, 2));
            macro = mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2SphereElement.modal(2, 1, reduction='symmetry'));
            micro = mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2SphereElement.nodal(3, [2, 4], reduction='symmetry'));
            record = mdemat.application.radiation.problem(spatial2, ...
                mdemat.common.fem.SumSpace(macro, micro), 0.7, 0, 0.01);
            result = mdemat.application.radiation.solve(record);
            testCase.verifyEqual(result.state.Dofs.F, ...
                0.7*ones(size(result.state.Dofs.F)), 'AbsTol', 2e-12);

            % A reduction label alone does not invalidate the same S^2 chart.
            micro.Element.DimReduction = 'topology';
            record = mdemat.application.radiation.problem(spatial2, ...
                mdemat.common.fem.SumSpace(macro, micro), 0.7, 0, 0.01);
            testCase.verifyEqual(record.velocitySpace.Rhs.Element.Geometry.BBox, ...
                record.velocitySpace.Lhs.Element.Geometry.BBox);
        end

        function radialVelocityKernelGradientsHaveTheCorrectDirection(testCase)
            nodes = [-0.4, 0.7; 0.2, -0.3];
            points = [-0.8, 0.1, 0.9; 0.3, -0.7, 0.4];
            for kernel = {'gaussian', 'multiquadric', 'inverse_multiquadric'}
                basis = mdemat.common.interpolation.RbfKernelBasis(nodes, ...
                    kernel=kernel{1}, epsilon=[1; 1.7]);
                for axis = 1:2
                    shift = zeros(2, 1);
                    shift(axis) = 1e-6;
                    difference = (basis.eval(points+shift)-basis.eval(points-shift))/(2e-6);
                    testCase.verifyEqual(basis.diff(points, double((1:2) == axis)), ...
                        difference, 'AbsTol', 3e-9);
                end
            end
        end
    end

    methods (Static, Access = private)
        function space = spatial(type, order, count)
            mesh = mdemat.data.mesh.UniformGrid([0, 1], count);
            if strcmp(type, 'modal')
                element = mdemat.common.fem.BH1OrthotopeElement.modal(1, order);
            else
                element = mdemat.common.fem.BH1OrthotopeElement.nodal(1, order);
            end
            space = mdemat.common.fem.FiniteElementSpace(mesh, element);
        end

        function space = velocity(dimension, mode)
            macro = [];
            micro = [];
            if ~strcmp(mode, 'micro')
                if dimension == 1
                    count = 1;
                    reduction = 'topology';
                else
                    count = 2;
                    reduction = 'topology';
                end
                macro = mdemat.common.fem.SpectralSpace( ...
                    mdemat.common.fem.L2SphereElement.modal(dimension, count, reduction=reduction));
            end
            if ~strcmp(mode, 'macro')
                counts = {2, 8, [3, 6]};
                micro = mdemat.common.fem.SpectralSpace( ...
                    mdemat.common.fem.L2SphereElement.nodal(dimension, counts{dimension}, ...
                        reduction='topology'));
            end
            space = mdemat.common.fem.SumSpace(macro, micro);
        end
    end
end
