classdef TestSolverContracts < matlab.unittest.TestCase
    %TESTSOLVERCONTRACTS Independent field/history/ownership regressions.
    methods (Test)
        function poissonPreservesEveryGradientComponent(testCase)
            space = TestSolverContracts.makeSpace(2);
            state = mdemat.common.fem.SpatialState(space);
            scheme = mdemat.solver.elliptic.LdgScheme(config=struct( ...
                'bc', @(x) zeros(1, size(x,2))));
            ng = space.NGlobalDofs;
            expected = [(1:ng).', 10+(1:ng).', -20-(1:ng).'];
            scheme.A = speye(numel(expected));
            scheme.b = expected(:);
            state = scheme.step(state);
            testCase.verifyEqual(state.Dofs.U, expected(:, 1));
            testCase.verifyEqual(state.Dofs.Q, expected(:, 2:end));
        end

        function poissonNeedsNoReferenceSolution(testCase)
            ldg = mdemat.solver.elliptic.LdgScheme();
            pdg = mdemat.solver.elliptic.PdgScheme();
            testCase.verifyTrue(all(structfun(@isempty, ldg.Observation.Analysis.Exacts)));
            testCase.verifyTrue(all(structfun(@isempty, pdg.Observation.Analysis.Exacts)));
            testCase.verifyEmpty(ldg.Observers.timer);
            testCase.verifyEmpty(ldg.Observers.visualizer);
            testCase.verifyEmpty(pdg.Observers.timer);
            testCase.verifyEmpty(pdg.Observers.visualizer);
            testCase.verifyTrue(ldg.Config.hasAuxiliaryVariable);
        end

        function poissonRecoversManufacturedQuadraticAndItsGradient(testCase)
            space = TestSolverContracts.makeSpace(2);
            solution = @(x) x(1, :).^2 + 2*x(2, :).^2 ...
                + 3*x(1, :) - 4*x(2, :) + 0.5;
            gradient = @(x) [2*x(1, :)+3; 4*x(2, :)-4];
            config = struct('bc', solution, ...
                'force', @(x) -6*ones(1, size(x, 2)));
            ldg = mdemat.solver.elliptic.LdgScheme(config=config);
            state = mdemat.common.fem.SpatialState(space);
            ldg.initialize(state);
            ldg.step(state);
            nodes = space.Element.Volume.Nodes;
            testCase.verifyEqual(space.eval(nodes, state.Dofs.U), ...
                space.feval(solution, nodes), 'AbsTol', 2e-10);
            testCase.verifySize(state.Dofs.Q, [space.NGlobalDofs, 2]);
            testCase.verifyEqual(space.eval(nodes, state.Dofs.Q), ...
                space.feval(gradient, nodes), 'AbsTol', 2e-10);

            pdg = mdemat.solver.elliptic.PdgScheme(config=config);
            primal = mdemat.common.fem.SpatialState(space);
            pdg.initialize(primal);
            pdg.step(primal);
            testCase.verifyEqual(primal.Dofs.U, state.Dofs.U, 'AbsTol', 2e-10);
        end

        function periodicPoissonUsesPhysicalMeanOnNonuniformMesh(testCase)
            mesh = mdemat.data.mesh.Grid({[0, 0.1, 0.35, 0.8, 1]});
            element = mdemat.common.fem.BH1OrthotopeElement.nodal(1, 5);
            space = mdemat.common.fem.FiniteElementSpace(mesh, element);
            % u and u_x match at both ends, and integral(u)=0 exactly.
            solution = @(x) x(1, :).^2.*(1-x(1, :)).^2-1/30;
            gradient = @(x) 2*x(1, :)-6*x(1, :).^2+4*x(1, :).^3;
            config = struct('force', @(x) -2+12*x(1, :)-12*x(1, :).^2);
            nodes = element.Volume.Nodes;
            state = mdemat.common.fem.SpatialState(space);
            ldg = mdemat.solver.elliptic.LdgScheme(config=config);
            ldg.initialize(state);
            ldg.step(state);
            values = space.eval(nodes, state.Dofs.U);
            testCase.verifyEqual(values, space.feval(solution, nodes), ...
                'AbsTol', 2e-10);
            testCase.verifyEqual(space.eval(nodes, state.Dofs.Q), ...
                space.feval(gradient, nodes), 'AbsTol', 2e-10);
            integrals = element.Volume.Weights ...
                * reshape(values, element.Volume.NPoints, []);
            determinants = mesh.computeElementJacobianDeterminants();
            testCase.verifyEqual(integrals*determinants(:), 0, ...
                'AbsTol', 1e-12);
            pdg = mdemat.solver.elliptic.PdgScheme(config=config);
            primal = mdemat.common.fem.SpatialState(space);
            pdg.initialize(primal);
            pdg.step(primal);
            testCase.verifyEqual(primal.Dofs.U, state.Dofs.U, 'AbsTol', 2e-10);
        end

        function observationHistoryUsesTimeIntervalAndPackedOffsets(testCase)
            space = TestSolverContracts.makeSpace(1);
            ng = space.NGlobalDofs;
            probe = MolHistoryProbe(config=struct('dt', 0.1));
            probe.TDisc = HistoryIntegratorProbe();
            probe.TDisc.Timeline.Now = 0.1;
            probe.TDisc.Timeline.Next = 0.2;
            probe.TDisc.Timeline.Count = 2;
            recent = repmat([11, 21, 22], ng, 1);
            older = repmat([12, 31, 32], ng, 1);
            probe.TDisc.U0 = {recent(:), older(:); ...
                -recent(:), -older(:)};
            probe.setObservation('plot', struct('Components', struct('B', 2)));
            dofs = struct('A', zeros(ng, 1), 'B', zeros(ng, 2));
            [values, count] = probe.history(space, dofs);
            testCase.verifyEqual(count, 2);
            testCase.verifyFalse(isfield(values{1}, 'A'));
            testCase.verifyEqual(values{1}.B, recent(:, 2:end));
            testCase.verifyEqual(values{2}.B, older(:, 2:end));

            probe.TDisc.Timeline.Next = probe.TDisc.Timeline.Now;
            [values, count] = probe.history(space, dofs);
            testCase.verifyEmpty(values);
            testCase.verifyEqual(count, 0);
        end

        function radiationCopyPreservesStateWithoutAliasingWeights(testCase)
            state = TestSolverContracts.makeMicroState();
            ng = state.NXDofs;
            state.setDecomposition('full');
            state.setScaling(0.3);
            state.setDof('G', reshape(1:(ng*state.NVMicroDofs), ng, []));
            state.setCoefficient('CS', 2);
            state.setHistory('time', [0, 0.1]);
            copied = state.copy();
            testCase.verifyEqual(copied.Dofs, state.Dofs);
            testCase.verifyEqual(copied.Coefs, state.Coefs);
            testCase.verifyEqual(copied.History, state.History);
            testCase.verifyEqual(copied.Decomposition, state.Decomposition);
            testCase.verifyEqual(copied.Levels, state.Levels);
            copied.VDisc.setRhsWeight(0.8 * ones(state.NVMicroDofs, 1));
            copied.Dofs.G(1) = -1;
            testCase.verifyEqual(state.VDisc.RhsWeight, ...
                0.3 * ones(state.NVMicroDofs, 1));
            testCase.verifyEqual(state.Dofs.G(1), 1);
            testCase.verifyEqual(double(copied.Status), 0);
        end

        function radiationMicroReconstructionAppliesScaling(testCase)
            state = TestSolverContracts.makeMicroState();
            ng = state.NXDofs;
            nv = state.NVMicroDofs;
            state.setDof('G', reshape(1:(ng*nv), ng, nv));
            state.VDisc.setRhsWeight((1:nv).'/10);
            x = state.XDisc.Element.Volume.Nodes;
            v = state.VDisc.Rhs.Element.Volume.Nodes;
            spatialValues = state.XDisc.eval(x, state.Dofs.G);
            basis = state.VDisc.Rhs.Element.Approximator.Basis.eval(v);
            expected = (spatialValues .* state.VDisc.RhsWeight.') * basis;
            testCase.verifyEqual(state.distribution(x, v), expected, ...
                'AbsTol', 1e-12);
            % A stale reconstructed cache must not change the physical density.
            state.Dofs.F = zeros(size(state.Dofs.G));
            weights = state.VDisc.Rhs.Element.Volume.Weights;
            testCase.verifyEqual(state.density(x), expected * weights(:), ...
                'AbsTol', 1e-12);
        end

        function radiationMicroFitAndReconstructionAreInverse(testCase)
            state = TestSolverContracts.makeMicroState();
            state.VDisc.setRhsWeight((1:state.NVMicroDofs).'/7);
            initial = @(x, v) 1+x(1, :)+v(1, :).^2;
            state.Dofs.G = state.microFit(initial);
            nodes = state.XDisc.Element.Volume.Nodes;
            velocity = state.VDisc.Rhs.Element.Volume.Nodes;
            points = state.XDisc.Mesh.collocate(nodes);
            expected = 1+points(1, :).'+velocity(1, :).^2;
            testCase.verifyEqual(state.distribution(nodes, velocity), ...
                expected, 'AbsTol', 1e-12);
        end

        function derivedRadiationHistoryIsReconstructedWithoutStateMutation(testCase)
            state = TestSolverContracts.makeMicroState();
            ng = state.NXDofs;
            nv = state.NVMicroDofs;
            state.Dofs.G = ones(ng, nv);
            state.kineticReconstruct();
            current = state.Dofs;
            probe = MolHistoryProbe(config=struct('dt', 0.1));
            probe.TDisc = HistoryIntegratorProbe();
            probe.TDisc.Timeline.Now = 0.1;
            probe.TDisc.Timeline.Next = 0.2;
            probe.TDisc.Timeline.Count = 1;
            historical = 3*ones(ng, nv);
            probe.TDisc.U0{1, 1} = historical(:);
            [history, count] = probe.history(state.XDisc, state.Dofs, ...
                reconstruct=@(packed) state.historyDofs(packed));
            testCase.verifyEqual(count, 1);
            testCase.verifyEqual(history{1}.G, historical);
            testCase.verifyEqual(history{1}.F, 3*current.F, 'AbsTol', 1e-12);
            testCase.verifyEqual(state.Dofs, current);
        end

        function stateCoefficientsAcceptEvaluationCallbacks(testCase)
            state = mdemat.data.variable.State();
            callback = @(x) x.^2;
            state.setCoefficient('source', callback);
            testCase.verifyEqual(state.Coefs.source(3), 9);
            testCase.verifyError(@() state.setCoefficient('bad', 'text'), ...
                'MATLAB:validators:mustBeA');
        end
    end

    methods (Static, Access = private)
        function space = makeSpace(dimensions)
            mesh = mdemat.data.mesh.UniformGrid(repmat([0, 1], 1, dimensions), ...
                2*ones(1, dimensions));
            element = mdemat.common.fem.BH1OrthotopeElement.nodal(dimensions, 3);
            space = mdemat.common.fem.FiniteElementSpace(mesh, element);
        end

        function state = makeMicroState()
            xSpace = TestSolverContracts.makeSpace(1);
            velocity = mdemat.common.fem.SpectralSpace(xSpace.Element);
            vSpace = mdemat.common.fem.SumSpace([], velocity);
            state = mdemat.application.radiation.MacroMicroState(xSpace, vSpace);
        end
    end
end
