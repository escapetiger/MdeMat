classdef TestIntegratorPersistence < matlab.unittest.TestCase
    methods (Test)
        function variableStepBdfRetainsMassAndSolutionHistory(testCase)
            integrator = mdemat.common.ode.Bdf3Integrator(0.5);
            integrator.Solver = mdemat.common.linalg.LinearSolver(nnzTh=20, ...
                condTh=1e8, maxIter=127, tol=2e-10, restart=11);
            steps = [0.04, 0.07, 0.03, 0.06, 0.08];
            initial = [1; 2i];
            value = initial;
            for dt = steps(1:3), value = advanceBdf(integrator, value, dt); end
            record = mdemat.common.ode.snapshotIntegrator(integrator);
            restored = mdemat.common.ode.restoreIntegrator(record);
            testCase.verifyEqual(restored.Startup, integrator.Startup);
            testCase.verifyEqual(restored.Solver.NnzTh, 20);
            testCase.verifyEqual(restored.Solver.CondTh, 1e8);
            testCase.verifyEqual(restored.Solver.KrylovOpts, integrator.Solver.KrylovOpts);
            testCase.verifyEmpty(restored.Solver.Lhs);
            testCase.verifyEqual(restored.U0, integrator.U0);
            testCase.verifyEqual(restored.MRhs, integrator.MRhs);
            testCase.verifyEqual(restored.Timeline.StepSizeQueue, integrator.Timeline.StepSizeQueue);
            testCase.verifyFalse(any(restored.Timeline.OldStepSizeQueue));
            continued = value;
            for dt = steps(4:end)
                value = advanceBdf(integrator, value, dt);
                continued = advanceBdf(restored, continued, dt);
            end
            testCase.verifyEqual(continued, value, 'AbsTol', 2e-13);
            testCase.verifyEqual(continued, initial/(1+sum(steps)), 'AbsTol', 2e-13);
        end

        function recursiveSplitRestartMatchesUninterruptedTrajectory(testCase)
            for kind = {'lie', 'strang'}
                factors = {mdemat.common.ode.BeIntegrator(0.5), ...
                    mdemat.common.ode.Exrk4Integrator(0.5)};
                if strcmp(kind{1}, 'lie')
                    original = mdemat.common.ode.LieTrotterIntegrator(factors=factors);
                else, original = mdemat.common.ode.StrangIntegrator(factors=factors); end
                value = [1; 2i];
                for index = 1:3, value = advanceSplit(original, value); end
                record = mdemat.common.ode.snapshotIntegrator(original);
                restored = mdemat.common.ode.restoreIntegrator(record);
                testCase.verifyClass(restored, class(original));
                testCase.verifyEqual(restored.Timeline.Now, original.Timeline.Now);
                for index = 1:original.NFactors
                    testCase.verifyEqual(restored.getFactor(index).U0, original.getFactor(index).U0);
                    testCase.verifyEqual(restored.getFactor(index).Timeline.Now, ...
                        original.getFactor(index).Timeline.Now);
                end
                continued = value;
                for index = 1:4
                    value = advanceSplit(original, value);
                    continued = advanceSplit(restored, continued);
                end
                testCase.verifyEqual(continued, value, 'AbsTol', 2e-13);
                testCase.verifyEqual(restored.Timeline.Now, original.Timeline.Now);
            end
        end

        function allRegisteredConcreteMethodsAndSettingsReconstruct(testCase)
            factories = {@mdemat.common.ode.BeIntegrator, @mdemat.common.ode.FeIntegrator, ...
                @mdemat.common.ode.Bdf2Integrator, @mdemat.common.ode.Bdf3Integrator, ...
                @mdemat.common.ode.EmpIntegrator, @mdemat.common.ode.HeunIntegrator, ...
                @mdemat.common.ode.Ssprk2Integrator, @mdemat.common.ode.Ssprk3Integrator, ...
                @mdemat.common.ode.Exrk4Integrator, @mdemat.common.ode.ImpIntegrator, ...
                @mdemat.common.ode.Sdirk2Integrator, @mdemat.common.ode.Sdirk3Integrator, ...
                @mdemat.common.ode.Sdirk4Integrator, @mdemat.common.ode.Esdirk3Integrator, ...
                @mdemat.common.ode.Ars111Integrator, @mdemat.common.ode.Ars222Integrator, ...
                @mdemat.common.ode.Ars443Integrator};
            for index = 1:numel(factories)
                original = factories{index}(1);
                original.NonlinearTolerance = 2e-9;
                original.NonlinearMaxIterations = 17;
                restored = mdemat.common.ode.restoreIntegrator( ...
                    mdemat.common.ode.snapshotIntegrator(original));
                testCase.verifyClass(restored, class(original));
                testCase.verifyEqual(restored.NonlinearTolerance, 2e-9);
                testCase.verifyEqual(restored.NonlinearMaxIterations, 17);
            end
        end

        function unknownMethodsPhasesAndCorruptedQueuesAreRejected(testCase)
            record = mdemat.common.ode.snapshotIntegrator(mdemat.common.ode.Bdf2Integrator(1));
            bad = record; bad.method = 'arbitrary.class';
            testCase.verifyError(@() mdemat.common.ode.restoreIntegrator(bad), ...
                'mdemat:ode:UnsupportedCheckpoint');
            bad = record; bad.phase = 'inside-stage';
            testCase.verifyError(@() mdemat.common.ode.restoreIntegrator(bad), ...
                'mdemat:ode:InvalidCheckpoint');
            bad = record; bad.timeline.StepSizeQueue = [-1, 0, 0];
            testCase.verifyError(@() mdemat.common.ode.restoreIntegrator(bad), ...
                'mdemat:ode:InvalidCheckpoint');
            bad = record; bad.coefficients.Solver = 1;
            testCase.verifyError(@() mdemat.common.ode.restoreIntegrator(bad), ...
                'mdemat:ode:InvalidCheckpoint');
        end
    end
end

function value = advanceBdf(integrator, value, dt)
integrator.setTimeStep(dt=dt);
integrator.update(value);
value = integrator.step(M=@(t) (1+t)*speye(2));
integrator.advance();
end

function value = advanceSplit(integrator, value)
integrator.setTimeStep(dt=0.05);
integrator.update(value);
value = integrator.step(L={[0, 1; -1, 0], [-0.3, 0.4; 0, -0.7]});
integrator.advance();
end
