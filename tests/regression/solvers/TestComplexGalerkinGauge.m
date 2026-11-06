classdef TestComplexGalerkinGauge < matlab.unittest.TestCase
    %TESTCOMPLEXGALERKINGAUGE Check basis covariance, not just complex storage.
    methods (Test)
        function physicalFieldsAndAllWeakOperatorsRespectAUnitaryPhase(testCase)
            for modal = [false, true]
                [realSpace, phaseSpace] = TestComplexGalerkinGauge.spaces(modal, false, false);
                field = @(x) x(1, :).^2+0.3*x(1, :)+0.7;
                original = realSpace.fit(field);
                changed = phaseSpace.fit(field);
                testCase.verifyEqual(changed, -1i*original, 'AbsTol', 2e-12);
                nodes = realSpace.Element.Volume.Nodes;
                testCase.verifyEqual(phaseSpace.eval(nodes, changed), ...
                    realSpace.eval(nodes, original), 'AbsTol', 2e-12);
                testCase.verifyEqual(phaseSpace.Element.Approximator.Mass, ...
                    realSpace.Element.Approximator.Mass, 'AbsTol', 2e-12);

                penalty = {'boundary', 'right'; 'left', ''};
                old = mdemat.common.fem.AdjointAssembly(realSpace, 'dirichlet', penalty);
                new = mdemat.common.fem.AdjointAssembly(phaseSpace, 'dirichlet', penalty);
                for callback = [false, true]
                    if callback, coefficient = @(x) 2+0*x(1, :); else, coefficient = 2; end
                    a = old.assembleMatrix(coefficient);
                    b = new.assembleMatrix(coefficient);
                    TestComplexGalerkinGauge.compareBlocks(testCase, b, a, 1);
                    a = old.assembleVector(coefficient, field);
                    b = new.assembleVector(coefficient, field);
                    TestComplexGalerkinGauge.compareBlocks(testCase, b, a, -1i);
                end
                weak = new.assembleMatrix(1);
                boundary = new.assembleVector(1, field);
                gradient = weak{2, 1}{1}*changed+boundary{2, 1}{1};
                testCase.verifyEqual(phaseSpace.eval(nodes, gradient), ...
                    phaseSpace.feval(@(x) 2*x(1, :)+0.3, nodes), 'AbsTol', 1e-11);

                old = mdemat.common.fem.UpwindAssembly(realSpace, 'dirichlet');
                new = mdemat.common.fem.UpwindAssembly(phaseSpace, 'dirichlet');
                testCase.verifyEqual(full(new.assembleMatrix(1)), full(old.assembleMatrix(1)), 'AbsTol', 1e-11);
                testCase.verifyEqual(full(new.assembleVector(1, 1, field)), ...
                    -1i*full(old.assembleVector(1, 1, field)), 'AbsTol', 1e-11);

                old = mdemat.common.fem.EllipticAssembly(realSpace, 'dirichlet', penalty);
                new = mdemat.common.fem.EllipticAssembly(phaseSpace, 'dirichlet', penalty);
                testCase.verifyEqual(full(new.assembleMatrix(1)), full(old.assembleMatrix(1)), 'AbsTol', 1e-11);
                testCase.verifyEqual(full(new.assembleVector(1, field)), ...
                    -1i*full(old.assembleVector(1, field)), 'AbsTol', 1e-11);

                old = mdemat.common.fem.WeightedAdjointAssembly(realSpace, 'dirichlet', penalty);
                new = mdemat.common.fem.WeightedAdjointAssembly(phaseSpace, 'dirichlet', penalty);
                weight = @(x) exp(0.2*x(1, :));
                derivative = @(x) 0.2*exp(0.2*x(1, :));
                a = old.assembleMatrix(1, weight, dWeight=derivative);
                b = new.assembleMatrix(1, weight, dWeight=derivative);
                TestComplexGalerkinGauge.compareBlocks(testCase, b, a, 1);
                TestComplexGalerkinGauge.compareBlocks(testCase, ...
                    new.assembleVector(1, weight, field), old.assembleVector(1, weight, field), -1i);
            end
        end

        function nonuniformPhysicalGaugeAndManufacturedSolutionAreCovariant(testCase)
            [~, space] = TestComplexGalerkinGauge.spaces(false, false, true);
            source = mdemat.common.fem.SourceAssembly(space);
            assembly = mdemat.common.fem.EllipticAssembly(space, 'periodic', ...
                {'boundary', 'right'; 'left', ''});
            matrix = -assembly.assembleMatrix(1);
            force = source.assembleVector(@(x) -2+12*x(1, :)-12*x(1, :).^2);
            [matrix, force] = assembly.forceCompletion(matrix, force);
            packed = matrix\force;
            values = space.eval(space.Element.Volume.Nodes, packed(1:end-1));
            exact = space.feval(@(x) x(1, :).^2.*(1-x(1, :)).^2-1/30, space.Element.Volume.Nodes);
            testCase.verifyEqual(values, exact, 'AbsTol', 2e-10);
            testCase.verifyEqual(space.integrationWeights().'*packed(1:end-1), 0, 'AbsTol', 1e-12);
        end

        function semiLagrangianTransferAndInflowRespectTheSamePhase(testCase)
            [realSpace, phaseSpace] = TestComplexGalerkinGauge.spaces(false, true, false);
            for boundary = ["periodic", "dirichlet"]
                old = mdemat.common.fem.SemiLagrangianAssembly(realSpace, char(boundary));
                new = mdemat.common.fem.SemiLagrangianAssembly(phaseSpace, char(boundary));
                testCase.verifyEqual(full(new.assembleMatrix(-1, 0.15, 0)), ...
                    full(old.assembleMatrix(-1, 0.15, 0)), 'AbsTol', 1e-11);
                callback = @(x, varargin) 1+x(1, :);
                testCase.verifyEqual(full(new.assembleVector(-1, 0.15, 0, callback)), ...
                    -1i*full(old.assembleVector(-1, 0.15, 0, callback)), 'AbsTol', 1e-11);
            end
        end
    end

    methods (Static, Access = private)
        function compareBlocks(testCase, actual, expected, phase)
            for row = 1:2
                for column = 1:2
                    if row == 2 && column == 2, continue; end
                    testCase.verifyEqual(full(actual{row, column}{1}), ...
                        phase*full(expected{row, column}{1}), 'AbsTol', 2e-11);
                end
            end
        end

        function [original, changed] = spaces(modal, semiLagrangian, nonuniform)
            if nonuniform
                mesh = mdemat.data.mesh.Grid({[0, 0.12, 0.6, 1]});
                count = 5;
            else
                mesh = mdemat.data.mesh.UniformGrid([0, 1], 4);
                count = 3;
            end
            if semiLagrangian
                base = mdemat.common.fem.BH1SLOrthotopeElement.nodal(1, count);
            elseif modal
                base = mdemat.common.fem.BH1OrthotopeElement.modal(1, count);
            else
                base = mdemat.common.fem.BH1OrthotopeElement.nodal(1, count);
            end
            original = mdemat.common.fem.FiniteElementSpace(mesh, base);
            phase = PhaseBasis(base.Approximator.Basis, 1i);
            volume = base.Volume.copy();
            flux = arrayfun(@(index) base.Flux(index).copy(), 1:base.NFluxes);
            if modal
                approximator = mdemat.common.approximation.ModalApproximator(phase);
                approximator.setMass(volume.Nodes, volume.Weights);
                approximator.setDesign(volume.Nodes, volume.Weights);
            else
                approximator = mdemat.common.approximation.NodalApproximator(phase);
                approximator.setMass(volume.Nodes, volume.Weights);
                approximator.setDesign(volume.Nodes);
            end
            if semiLagrangian
                element = mdemat.common.fem.BH1SLOrthotopeElement(base.Geometry, approximator, volume, flux);
            else
                element = mdemat.common.fem.BH1OrthotopeElement(base.Geometry, approximator, volume, flux);
            end
            changed = mdemat.common.fem.FiniteElementSpace(mesh, element);
        end
    end
end
