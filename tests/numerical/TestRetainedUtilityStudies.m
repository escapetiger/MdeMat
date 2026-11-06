classdef TestRetainedUtilityStudies < matlab.unittest.TestCase
    properties (TestParameter)
        entry = {'ode/ex01_prothero_robinson', 'ode/ex02_heat_equation', ...
            'ode/ex03_convergence_analysis', 'ode/ex04_van_der_pol', ...
            'ode/ex05_stability_region', 'ode/ex06_custom_step_functions', ...
            'geometry/ex01_graph_visualization', ...
            'geometry/ex02_separable_graph_visualization', 'geometry/ex03_uniform_trimesh', ...
            'geometry/ex04_nonuniform_trimesh'}
    end

    methods (Test)
        function actualRetainedEntryExecutesItsNumericsAndRendering(testCase, entry)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            source = fullfile(root, 'toolbox', 'examples', entry);
            [folder, name] = fileparts(source);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(folder));
            factory = str2func(name);
            testCase.assertEqual(which(name), [source, '.m']);
            visible = get(groot, 'DefaultFigureVisible');
            testCase.addTeardown(@() set(groot, 'DefaultFigureVisible', visible));
            set(groot, 'DefaultFigureVisible', 'off');
            original = findall(groot, 'Type', 'figure');
            testCase.addTeardown(@() TestRetainedUtilityStudies.deleteCreatedFigures(original));
            randomState = rng;
            if startsWith(entry, 'ode/')
                result = factory(Smoke=true, MakePlots=false);
                testCase.verifyEqual(findall(groot, 'Type', 'figure'), original);
            else
                sentinel = figure('Visible', 'off', 'Name', 'Caller figure');
                result = factory();
                testCase.verifyTrue(isgraphics(sentinel));
                testCase.verifyEqual(get(sentinel, 'Name'), 'Caller figure');
                testCase.verifyGreaterThan(numel(findall(groot, 'Type', 'figure')), ...
                    numel(original)+1);
            end
            testCase.verifyEqual(result.schemaVersion, 1);
            testCase.verifyEqual(rng, randomState);
            switch name
                case 'ex01_prothero_robinson'
                    testCase.verifyGreaterThan(numel(result.names), 15);
                    testCase.verifyTrue(all(isfinite(result.finalErrors)));
                    testCase.verifyLessThan(max(result.finalErrors), 0.1);
                    testCase.verifyEqual(cellfun(@(t) t(end), result.times), ...
                        result.finalTime*ones(size(result.times)), 'AbsTol', 1e-14);
                case 'ex02_heat_equation'
                    error = max(abs(result.solutions-result.semiDiscreteExact.'), [], 2);
                    testCase.verifyLessThan(max(error), 1e-3);
                    testCase.verifyEqual(result.x, (1:8).'/9, 'AbsTol', 1e-15);
                case 'ex03_convergence_analysis'
                    testCase.verifyTrue(all(cellfun(@(x) all(isfinite(x)), result.errors)));
                    index = find(strcmp(result.names, 'Be'));
                    testCase.verifyGreaterThan(result.observedOrders(index), 0.8);
                    testCase.verifyLessThan(result.observedOrders(index), 1.2);
                case 'ex04_van_der_pol'
                    derivative = @(t, u) [u(2); result.mu*(1-u(1)^2)*u(2)-u(1)];
                    [~, reference] = ode45(derivative, [0, result.finalTime], [2; 0], ...
                        odeset('RelTol', 1e-11, 'AbsTol', 1e-13));
                    for index = 1:numel(result.solutions)
                        testCase.verifyTrue(all(isfinite(result.solutions{index}), 'all'));
                        testCase.verifyLessThan(norm(result.solutions{index}(:, end) ...
                            -reference(end, :).', Inf), 1e-2);
                    end
                case 'ex05_stability_region'
                    column = find(abs(result.z+5) < eps, 1);
                    testCase.assertNotEmpty(column);
                    names = {'Be', 'Fe', 'Exrk4', 'Bdf2'};
                    expected = [1/6, -4, 1-5+25/2-125/6+625/24, 1/sqrt(13)];
                    for index = 1:numel(names)
                        method = find(strcmp(result.names, names{index}));
                        testCase.verifyEqual(result.amplification{method}(column), ...
                            expected(index), 'AbsTol', 2e-12);
                    end
                    rk = find(strcmp(result.testNames, 'Exrk4'));
                    testCase.verifyEqual(result.solutions(rk, :), ...
                        expected(3).^(0:2), 'AbsTol', 2e-11);
                case 'ex06_custom_step_functions'
                    h = result.stepSize;
                    testCase.verifyEqual(result.basic, (1+h*sin(h))/(1-result.lambda*h), ...
                        'AbsTol', 2e-14);
                    testCase.verifyEqual(result.custom, result.basic, 'AbsTol', 2e-14);
                    testCase.verifyEqual(result.values(3, :), ...
                        1+result.lambda*result.timeSteps, 'AbsTol', 2e-14);
                case 'ex01_graph_visualization'
                    testCase.verifyEqual(result.basic.NVertices, 6);
                    testCase.verifyEqual(result.basic.NEdges, 5);
                    testCase.verifyEqual(result.cube.NVertices, 8);
                    testCase.verifyEqual(result.cube.NEdges, 12);
                    testCase.verifyEqual(vecnorm(result.sphere.Vertices, 2, 2), ...
                        ones(result.sphere.NVertices, 1), 'AbsTol', 2e-14);
                case 'ex02_separable_graph_visualization'
                    full = result.grid.full();
                    testCase.verifyEqual(full.NVertices, prod(result.grid.NVertices));
                    testCase.verifyEqual(full.NEdges, result.grid.NTotalEdges);
                    testCase.verifyEqual(result.projected.NDims, 4);
                    testCase.verifyLessThan(result.separableVerticesBytes, result.fullVerticesBytes);
                case 'ex03_uniform_trimesh'
                    testCase.verifyEqual(result.mesh.NElements, 32);
                    testCase.verifyEqual(result.physicalArea, 1, 'AbsTol', 2e-14);
                    testCase.verifyEqual(nnz(result.mesh.FaceToFaceTable), 16);
                case 'ex04_nonuniform_trimesh'
                    testCase.verifyEqual(result.refined.NElements, 4*result.mesh.NElements);
                    testCase.verifyLessThan(abs(result.physicalArea-pi), 0.03);
                    testCase.verifyGreaterThan(result.quality.minQuality, 0.2);
            end
        end
    end

    methods (Static, Access = private)
        function deleteCreatedFigures(original)
            figures = findall(groot, 'Type', 'figure');
            delete(setdiff(figures, original));
        end
    end
end
