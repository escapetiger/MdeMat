classdef TestApplicationPlasma < matlab.unittest.TestCase
    methods (Test)
        function homogeneousMaxwellianIsStationary(testCase)
            for method = {'ap', 'dghs'}
                for temperature = [0.7, 1.5]
                    [space, velocity] = TestApplicationPlasma.spaces(3, temperature, 4);
                    initial = @(x, nh) [ones(1, size(x, 2)); zeros(nh-1, size(x, 2))];
                    config = struct('tau', 0.4, 'tau0', 0.4);
                    problem = mdemat.application.plasma.problem(space, velocity, initial, ...
                        0.02, 0.01, method=method{1}, config=config);
                    [result, scheme] = mdemat.application.plasma.solve(problem);
                    v = linspace(-4, 4, 31)*sqrt(temperature);
                    exact = exp(-v.^2/(2*temperature))/sqrt(2*pi*temperature);
                    actual = result.state.distribution(0, v);
                    testCase.verifyEqual(actual, repmat(exact, 4, 1), 'AbsTol', 5e-11);
                    testCase.verifyEqual(result.diagnostics.mass(end), 2*pi, 'AbsTol', 3e-11);
                    testCase.verifyEqual(result.diagnostics.kineticEnergy(end), ...
                        pi*temperature, 'AbsTol', 3e-11);
                    testCase.verifyEqual(result.diagnostics.totalEnergy(end), ...
                        result.diagnostics.kineticEnergy(end)+ ...
                        result.diagnostics.potentialEnergy(end), 'AbsTol', 1e-13);
                    testCase.verifyEmpty(scheme.Observers.timer);
                    testCase.verifyEmpty(scheme.Observers.visualizer);
                    testCase.verifyEqual(scheme.Observation.Analysis.PointsPerDimension, 0);
                    testCase.verifyEqual(result.state.VBBox, [-Inf; Inf]);
                end
            end
        end

        function collisionRelaxationHasIndependentSecondOrderTimeReference(testCase)
            tau = 0.3;
            final = 0.1;
            for method = {'ap', 'dghs'}
                [space, velocity] = TestApplicationPlasma.spaces(4, 1, 3);
                initial = @(x, nh) repmat([1; 0; 0.1; zeros(nh-3, 1)], 1, size(x, 2));
                errors = zeros(1, 2);
                for index = 1:2
                    steps = [0.02, 0.01];
                    config = struct('tau', tau, 'tau0', tau);
                    problem = mdemat.application.plasma.problem(space, velocity, initial, ...
                        final, steps(index), method=method{1}, config=config);
                    result = mdemat.application.plasma.solve(problem);
                    actual = space.eval(0, result.state.Dofs.D(:, 3));
                    exact = 0.1*exp(-2*final/tau);
                    errors(index) = max(abs(actual-exact));
                    testCase.verifyEqual(result.diagnostics.mass(end), 2*pi, 'AbsTol', 3e-11);
                end
                testCase.verifyLessThan(errors(2), errors(1)/3);
                testCase.verifyLessThan(errors(2), 1e-4);
            end
        end

        function shiftedGaussianConvergesWithHermiteTruncation(testCase)
            shift = 0.4;
            errors = zeros(1, 3);
            counts = [3, 5, 8];
            for index = 1:numel(counts)
                [space, velocity] = TestApplicationPlasma.spaces(counts(index), 1, 3);
                initial = @(x, nh) repmat((shift.^(0:nh-1)./sqrt(factorial(0:nh-1))).', ...
                    1, size(x, 2));
                problem = mdemat.application.plasma.problem(space, velocity, initial, ...
                    0.02, 0.01, method='dghs', config=struct('tau0', Inf));
                result = mdemat.application.plasma.solve(problem);
                v = linspace(-4, 4, 101);
                exact = exp(-(v-shift).^2/2)/sqrt(2*pi);
                actual = result.state.distribution(0, v);
                errors(index) = max(abs(actual-exact), [], 'all');
                testCase.verifyEqual(result.diagnostics.mass(end), 2*pi, 'AbsTol', 3e-11);
                testCase.verifyEqual(result.diagnostics.momentum(end), 2*pi*shift, ...
                    'AbsTol', 3e-11);
                testCase.verifyEqual(result.diagnostics.kineticEnergy(end), pi*(1+shift^2), ...
                    'AbsTol', 3e-11);
            end
            testCase.verifyLessThan(errors(2), errors(1)/10);
            testCase.verifyLessThan(errors(3), errors(2)/10);
        end

        function poissonCouplingMatchesContinuousFieldAndPreservesMass(testCase)
            amplitude = 0.02;
            for method = {'ap', 'dghs'}
                [space, velocity] = TestApplicationPlasma.spaces(5, 1, 12);
                initial = @(x, nh) [1+amplitude*cos(x); zeros(nh-1, size(x, 2))];
                config = struct('tau', Inf, 'tau0', Inf);
                problem = mdemat.application.plasma.problem(space, velocity, initial, ...
                    0, 0.001, method=method{1}, config=config);
                initialResult = mdemat.application.plasma.solve(problem);
                points = space.Mesh.collocate(0).';
                if strcmp(method{1}, 'ap')
                    field = space.eval(0, initialResult.state.Dofs.Q);
                else
                    field = space.eval(0, initialResult.state.Dofs.E);
                end
                testCase.verifyLessThan(max(abs(field-amplitude*sin(points))), 2e-3);
                problem.finalTime = 0.01;
                result = mdemat.application.plasma.solve(problem);
                testCase.verifyEqual(result.diagnostics.mass(end), 2*pi, 'AbsTol', 2e-9);
                testCase.verifyTrue(all(isfinite(result.state.Dofs.D), 'all'));
                testCase.verifyLessThan(max(abs(result.diagnostics.mass-2*pi)), 2e-9);
            end
        end

        function nonuniformContinuousEquilibriumConvergesUnderSpatialRefinement(testCase)
            % T*rhoInf'=EInf*rhoInf and EInf'=rhoInf-rhoi hold continuously.
            amplitude = 0.1;
            errors = zeros(1, 2);
            counts = [8, 16];
            for index = 1:2
                [space, velocity] = TestApplicationPlasma.spaces(5, 1, counts(index));
                config = struct('epsilon', 1, 'tau', 0.4, ...
                    'phiInf', @(x) -amplitude*cos(x), ...
                    'rhoInf', @(x) exp(amplitude*cos(x)), ...
                    'sqrtRhoInf', @(x) exp(amplitude*cos(x)/2), ...
                    'EInf', @(x) -amplitude*sin(x), ...
                    'rhoi', @(x) exp(amplitude*cos(x))+amplitude*cos(x));
                initial = @(x, nh) [exp(amplitude*cos(x)/2); zeros(nh-1, size(x, 2))];
                problem = mdemat.application.plasma.problem(space, velocity, initial, ...
                    0.02, 0.002, method='ap', config=config);
                result = mdemat.application.plasma.solve(problem);
                v = [-1, 0, 1];
                positions = space.Mesh.collocate(0).';
                exact = exp(amplitude*cos(positions))*exp(-v.^2/2)/sqrt(2*pi);
                values = result.state.distribution(0, v);
                errors(index) = max(abs(values-exact), [], 'all');
                testCase.verifyEqual(result.diagnostics.mass(end), ...
                    2*pi*besseli(0, amplitude), 'AbsTol', 2e-3);
            end
            testCase.verifyLessThan(errors(2), errors(1)/2.5);
            testCase.verifyLessThan(errors(2), 2e-3);
        end

        function velocityRepresentationAndTemperatureAreConsistent(testCase)
            temperature = 0.7;
            [space, velocity] = TestApplicationPlasma.spaces(3, temperature, 2);
            initial = @(x, nh) [ones(1, size(x, 2)); zeros(nh-1, size(x, 2))];
            record = mdemat.application.plasma.problem(space, velocity, initial, ...
                0, 0.01, config=struct('T0', temperature));
            result = mdemat.application.plasma.solve(record);
            testCase.verifyEqual(result.state.distribution(0, 0), ...
                ones(2, 1)/sqrt(2*pi*temperature), 'AbsTol', 3e-12);
            testCase.verifyError(@() mdemat.application.plasma.problem( ...
                space, velocity, initial, 0, 0.01, config=struct('T0', 1)), ...
                'MdeMat:Plasma:TemperatureMismatch');
            for value = {0, -1, NaN, Inf, [temperature, temperature], @(x) temperature}
                testCase.verifyError(@() mdemat.application.plasma.problem( ...
                    space, velocity, initial, 0, 0.01, config=struct('T0', value{1})), ...
                    'MdeMat:Plasma:InvalidTemperature');
            end
            otherBasis = mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.BH1OrthotopeElement.modal(1, 3));
            testCase.verifyError(@() mdemat.application.plasma.problem( ...
                space, otherBasis, initial, 0, 0.01), ...
                'MdeMat:Plasma:UnsupportedRepresentation');
            factor = velocity.Element.Approximator.Basis.Lhs.getFactor(1);
            factor.Lhs.IsNormalized = false;
            testCase.verifyError(@() mdemat.application.plasma.problem( ...
                space, velocity, initial, 0, 0.01), ...
                'MdeMat:Plasma:UnsupportedRepresentation');
            factor.Lhs.IsNormalized = true;
            velocity.Element.Approximator.Basis.Rhs.Center = 0.1;
            testCase.verifyError(@() mdemat.application.plasma.problem( ...
                space, velocity, initial, 0, 0.01), ...
                'MdeMat:Plasma:UnsupportedRepresentation');
            velocity.Element.Approximator.Basis.Rhs.Center = 0;
            factor.Rhs.A = 2/sqrt(temperature);
            testCase.verifyError(@() mdemat.application.plasma.problem( ...
                space, velocity, initial, 0, 0.01), ...
                'MdeMat:Plasma:UnsupportedRepresentation');
        end

        function problemRecordsRejectUnknownMissingAndIncompatibleSchema(testCase)
            [space, velocity] = TestApplicationPlasma.spaces(3, 1, 2);
            initial = @(x, nh) [ones(1, size(x, 2)); zeros(nh-1, size(x, 2))];
            record = mdemat.application.plasma.problem(space, velocity, initial, 0, 0.01);
            extra = record;
            extra.unknown = 1;
            wrongModel = record;
            wrongModel.model = 'other';
            wrongSchema = record;
            wrongSchema.schemaVersion = 2;
            for invalid = {extra, rmfield(record, 'method'), wrongModel, wrongSchema}
                testCase.verifyError(@() mdemat.application.plasma.solve(invalid{1}), ...
                    'MdeMat:Plasma:InvalidProblem');
            end
        end

        function hermiteFiltersDampOnlyUnresolvedModes(testCase)
            for method = {'ap', 'dghs'}
                for filter = {'hou-li', 'tc', 'quasi-tc'}
                    [space, velocity] = TestApplicationPlasma.spaces(5, 1, 3);
                    initial = @(x, nh) repmat([1; zeros(nh-2, 1); 0.01], 1, size(x, 2));
                    config = struct('tau', Inf, 'tau0', Inf, 'useFilter', true, ...
                        'dealiasingType', filter{1}, 'dealiasingBeta', 4, ...
                        'dealiasingGamma', 4, 'dealiasingCutoff', 0.5);
                    problem = mdemat.application.plasma.problem(space, velocity, initial, ...
                        0.01, 0.01, method=method{1}, config=config);
                    result = mdemat.application.plasma.solve(problem);
                    highest = space.eval(0, result.state.Dofs.D(:, end));
                    testCase.verifyGreaterThanOrEqual(min(highest), -1e-12);
                    testCase.verifyLessThan(max(highest), 0.01);
                    testCase.verifyEqual(result.diagnostics.mass(end), 2*pi, 'AbsTol', 3e-11);
                end
            end
        end
    end

    methods (Static, Access = private)
        function [space, velocity] = spaces(count, temperature, cells)
            mesh = mdemat.data.mesh.UniformGrid([0, 2*pi], cells);
            element = mdemat.common.fem.BH1OrthotopeElement.nodal(1, 2);
            space = mdemat.common.fem.FiniteElementSpace(mesh, element);
            velocity = mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2OrthotopeElement.hermite(1, count, T=temperature));
        end
    end
end
