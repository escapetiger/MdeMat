classdef TestAssemblyIdentities < matlab.unittest.TestCase
    methods (Test)
        function componentAxesAndElementSubsets(testCase)
            mesh = mdemat.data.mesh.Grid({[-1, -0.7, 0.2, 1.5]});
            element = mdemat.common.fem.BH1OrthotopeElement.modal(1, 3);
            space = mdemat.common.fem.FiniteElementSpace(mesh, element);
            field = @(x) [2+3*x; -1+0.5*x];
            coefficients = space.fit(field);
            reference = [-0.4, 0.1, 0.35];
            points = mesh.collocate(reference);
            testCase.verifyEqual(space.eval(reference, coefficients), field(points).', ...
                'AbsTol', 3e-13);
            testCase.verifyEqual(space.feval(field, reference), field(points).', 'AbsTol', 2e-14);
            packed = @(x) [2+3*x, -1+0.5*x];
            legacy = @(x) [(2+3*x).', (-1+0.5*x).'];
            testCase.verifyEqual(space.fit(packed), coefficients, 'AbsTol', 3e-13);
            testCase.verifyEqual(space.fit(legacy), coefficients, 'AbsTol', 3e-13);
            selected = [3, 1];
            selectedPoints = mesh.collocate(reference, selected);
            expected = field(selectedPoints).';
            testCase.verifyEqual(space.eval(reference, coefficients, EI=selected), expected, ...
                'AbsTol', 3e-13);
            rows = space.localToGlobal(selected, []);
            testCase.verifyEqual(space.eval(reference, coefficients(rows, :), EI=selected), ...
                expected, 'AbsTol', 3e-13);
            derivative = repelem(mesh.Spacings{1}(selected).', numel(reference), 1)*[3, 0.5];
            testCase.verifyEqual(space.diff(reference, coefficients, 1, EI=selected), ...
                derivative, 'AbsTol', 3e-13);
            testCase.verifyEqual(space.feval([-2; 4], reference), ...
                repmat([-2, 4], numel(points), 1));
            testCase.verifyEqual(space.feval(field, 0, EI=selected), ...
                field(mesh.collocate(0, selected)).', 'AbsTol', 2e-14);
        end

        function discreteIntegrationByPartsIncludesFaceMeasure(testCase)
            for dimension = 1:2
                element = mdemat.common.fem.BH1OrthotopeElement.modal(dimension, 3);
                basis = element.Approximator.Basis;
                values = element.Volume.Values;
                weights = element.Volume.Weights;
                normals = element.Geometry.OutNormals;
                for axis = 1:dimension
                    order = double((1:dimension) == axis);
                    derivative = basis.diff(element.Volume.Nodes, order);
                    volume = (derivative.*weights)*values.';
                    boundary = zeros(element.NDofs);
                    for face = 1:element.NFluxes
                        testCase.verifyEqual(sum(element.Flux(face).Weights), 1, ...
                            'AbsTol', 3e-14);
                        trace = element.Flux(face).Values;
                        boundary = boundary+normals(axis, face)* ...
                            (trace.*element.Flux(face).Weights)*trace.';
                    end
                    testCase.verifyEqual(volume+volume.', boundary, 'AbsTol', 3e-13);
                end
            end
        end

        function sourceProjectionAndMultiplicationAreComplexLinear(testCase)
            geometry = mdemat.data.geometry.Orthotope([-0.5, 0.5]);
            volume = mdemat.common.quadrature.Integrator( ...
                mdemat.common.quadrature.GaussLegendreRule(1));
            volume.setPoints(geometry, 3);
            basis = mdemat.data.function.Linear([0; 1+1i], b=[1; 0]);
            approximation = mdemat.common.approximation.ModalApproximator(basis);
            approximation.setMass(volume.Nodes, volume.Weights);
            approximation.setDesign(volume.Nodes, volume.Weights);
            element = mdemat.common.fem.L2Element(geometry, approximation, volume);
            space = mdemat.common.fem.FiniteElementSpace( ...
                mdemat.data.mesh.Grid({[-1, -0.2, 1.5]}), element);
            field = @(x) [2+(1+1i)*x; -1+2i*x];
            source = mdemat.common.fem.SourceAssembly(space);
            testCase.verifyEqual(full(source.assembleVector(field)), space.fit(field), ...
                'AbsTol', 3e-13);
            multiplier = mdemat.common.fem.MultiplierAssembly(space);
            testCase.verifyEqual(full(multiplier.assembleMatrix(@(x) ones(size(x)))), ...
                eye(space.NGlobalDofs), 'AbsTol', 3e-13);
        end

        function physicalErrorNormsOnANonuniformMesh(testCase)
            mesh = mdemat.data.mesh.Grid({[-2, -1.7, -0.2, 1, 3]});
            element = mdemat.common.fem.BH1OrthotopeElement.modal(1, 2);
            space = mdemat.common.fem.FiniteElementSpace(mesh, element);
            plan = mdemat.common.fem.samplingPlan(space,4);
            rows = numel(plan.Weights);
            report = mdemat.tools.analyze.weightedError(zeros(2,rows), ...
                repmat([2;3],1,rows),plan.Weights);
            testCase.verifyEqual([report.ComponentL1,report.ComponentL2,report.ComponentLInf].', ...
                [5; sqrt(5); 1]*[2, 3], 'AbsTol', 2e-13);
            globalIndices = space.localToGlobal([], []);
            [elements, local] = space.globalToLocal(globalIndices);
            reconstructed = local+(elements-1)*space.NLocalDofs;
            testCase.verifyEqual(reconstructed, globalIndices);
        end

        function hermiteTemperatureUsesPhysicalVelocityNodes(testCase)
            for temperature = [0.4, 1, 2.3]
                element = mdemat.common.fem.L2OrthotopeElement.hermite(1, 4, T=temperature);
                distribution = mdemat.data.function.Maxwellian(1, temperature=temperature);
                density = distribution.eval(element.Volume.Nodes);
                probabilityWeights = element.Volume.Weights.*density.^2;
                testCase.verifyEqual(sum(probabilityWeights), 1, 'AbsTol', 3e-14);
                testCase.verifyEqual(probabilityWeights*element.Volume.Nodes(:).^2, ...
                    temperature, 'AbsTol', 2e-13);
                values = element.Volume.Values;
                gram = (values.*element.Volume.Weights)*values.';
                testCase.verifyEqual(gram, eye(element.NDofs), 'AbsTol', 3e-12);
            end
        end
    end
end
