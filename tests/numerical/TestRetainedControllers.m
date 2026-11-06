classdef TestRetainedControllers < matlab.unittest.TestCase
    % Retained study controllers must execute the same physical equations.
    methods (Test)
        function steadyRunRecoversValueAndFlux(testCase)
            space = TestRetainedControllers.space(3, false);
            exact = @(x) 1+x(1, :)+x(1, :).^2;
            scheme = mdemat.solver.elliptic.LdgScheme(config=struct( ...
                'bc', exact, 'force', @(x) -2+0*x(1, :)));
            state = scheme.run(mdemat.common.fem.SpatialState(space));
            nodes = space.Element.Volume.Nodes;
            testCase.verifyEqual(space.eval(nodes, state.Dofs.U), ...
                space.feval(exact, nodes), 'AbsTol', 3e-11);
            testCase.verifyEqual(space.eval(nodes, state.Dofs.Q), ...
                space.feval(@(x) 1+2*x(1, :), nodes), 'AbsTol', 3e-11);
            timer = scheme.getObserver('timer');
            timing = timer.find('Integration');
            testCase.verifyGreaterThanOrEqual(timing.Duration, 0);
        end

        function steadyConvergenceAssessesExactAndNestedSolutions(testCase)
            exact = @(x) 1+x(1, :)+sin(pi*x(1, :));
            for hasReference = [true, false]
                config = struct('bc', exact, 'force', @(x) pi^2*sin(pi*x(1, :)), ...
                    'components', struct('U', 1));
                if hasReference
                    config.exact = struct('U', exact);
                end
                scheme = mdemat.solver.elliptic.PdgScheme(config=config);
                scheme.setObservation('analysis',struct('PointsPerDimension',7,'NLevels',2));
                [state, tables] = scheme.converge(mdemat.common.fem.SpatialState( ...
                    TestRetainedControllers.space(2, false)));
                testCase.verifyGreaterThan(state.XDisc.NMeshElements, 2);
                testCase.verifyGreaterThan(tables.U1.L2_Error(1), 0);
                testCase.verifyLessThan(tables.U1.L2_Error(2), ...
                    tables.U1.L2_Error(1)/4);
                timer = scheme.getObserver('timer');
                assessment = timer.find('Analysis');
                testCase.verifyNotEmpty(assessment);
                testCase.verifyFalse(assessment.IsActive);
            end
        end

        function transportRunsPreserveTheAffineCharacteristic(testCase)
            for semiLagrangian = [false, true]
                space = TestRetainedControllers.space(3, semiLagrangian);
                config = struct('advection', -1, 'ic', @(x) 2+x(1, :), ...
                    'bc', @(x, ~, t) 2+x(1, :)+t, 'tFinal', 0.03, ...
                    'dt', 0.01, 'tOdeIntName', 'ssprk3');
                if semiLagrangian
                    scheme = mdemat.solver.hyperbolic.SldgScheme(config=config);
                else
                    scheme = mdemat.solver.hyperbolic.UwdgScheme(config=config);
                end
                state = scheme.run(mdemat.common.fem.SpatialState(space));
                nodes = space.Element.Volume.Nodes;
                testCase.verifyEqual(space.eval(nodes, state.Dofs.U), ...
                    space.feval(@(x) 2+x(1, :)+0.03, nodes), 'AbsTol', 5e-11);
                testCase.verifyEqual(scheme.TDisc.Timeline.Now, 0.03, 'AbsTol', eps);
            end
        end

        function diffusionRunsWithTimeDependentDirichletData(testCase)
            space = TestRetainedControllers.space(8, false);
            initial = @(x) 1+x(1, :)+0.2*x(1, :).^2+0.1*sin(pi*x(1, :));
            exact = @(x, t) 1+x(1, :)+0.2*(x(1, :).^2+2*t) ...
                +0.1*exp(-pi^2*t)*sin(pi*x(1, :));
            config = struct('diffusion', 1, 'ic', initial, ...
                'bc', @(x, ~, t) exact(x, t), ...
                'tFinal', 0.002, 'dt', 0.0005, 'tOdeIntName', 'sdirk2');
            scheme = mdemat.solver.parabolic.LdgScheme(config=config);
            state = scheme.run(mdemat.common.fem.SpatialState(space));
            nodes = space.Element.Volume.Nodes;
            testCase.verifyEqual(space.eval(nodes, state.Dofs.U), ...
                space.feval(@(x) exact(x, 0.002), nodes), 'AbsTol', 2e-4);
        end

        function disabledConvergenceReturnsWithoutSolving(testCase)
            state = mdemat.common.fem.SpatialState(TestRetainedControllers.space(1, false));
            scheme = mdemat.solver.elliptic.PdgScheme();
            [unchanged, result] = scheme.converge(state);
            testCase.verifyTrue(unchanged == state);
            testCase.verifyEmpty(result);
            testCase.verifyEmpty(fieldnames(state.Dofs));
        end

        function emptySelectedComponentsDisableConvergence(testCase)
            state = mdemat.common.fem.SpatialState(TestRetainedControllers.space(1,false));
            scheme = mdemat.solver.elliptic.PdgScheme();
            scheme.setObservation('analysis',struct('PointsPerDimension',3, ...
                'Components',struct('U',[])));
            [unchanged,result] = scheme.converge(state);
            testCase.verifySameHandle(unchanged,state);
            testCase.verifyEmpty(result);
            testCase.verifyEmpty(fieldnames(state.Dofs));
            testCase.verifyEmpty(scheme.Observers.timer);
        end

        function modelAdaptedExactReferencesDoNotAddAnInterlevelSolve(testCase)
            spatial = TestRetainedControllers.space(3,false);
            velocity = mdemat.common.fem.SumSpace( ...
                mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2SphereElement.modal(1,1,reduction='topology')), ...
                mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2SphereElement.nodal(1,2,reduction='topology')));
            config = struct('nDims',1,'nu',1,'nv',2,'E',0.5, ...
                'vDimReduction','topology','ic',1, ...
                'exact',@(x,v,varargin) ones(1,size(x,2)), ...
                'dt',0.001,'cfl',[],'tFinal',0.002);
            scheme = mdemat.application.radiation.MmDgScheme(config=config);
            scheme.setObservation('analysis',struct('PointsPerDimension',3, ...
                'NLevels',1,'Components',struct('U',1)));
            [state,result] = scheme.converge( ...
                mdemat.application.radiation.MacroMicroState(spatial,velocity));
            testCase.verifyEqual(state.XDisc.NMeshElements,3);
            testCase.verifySize(result.U1.L2_Error,[1,1]);
            testCase.verifyLessThan(result.U1.L2_Error,2e-12);
        end
    end

    methods (Static, Access = private)
        function space = space(cells, semiLagrangian)
            mesh = mdemat.data.mesh.UniformGrid([0, 1], cells);
            if semiLagrangian
                element = mdemat.common.fem.BH1SLOrthotopeElement.nodal(1, 3);
            else
                element = mdemat.common.fem.BH1OrthotopeElement.nodal(1, 3);
            end
            space = mdemat.common.fem.FiniteElementSpace(mesh, element);
        end
    end
end
