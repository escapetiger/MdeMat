classdef TestSchemeCheckpoint < matlab.unittest.TestCase
    methods (Test)
        function resumeRunnerPreservesAcceptedValuesAcrossInitialization(testCase)
            config = struct('tFinal', 0.3, 'dt', 0.1);
            scheme = ResumeSchemeProbe(config);
            space = mdemat.common.fem.FiniteElementSpace( ...
                mdemat.data.mesh.UniformGrid([0, 1], 2), ...
                mdemat.common.fem.BH1OrthotopeElement.modal(1, 2));
            state = mdemat.common.fem.SpatialState(space);
            initial = (1:space.NGlobalDofs).';
            state.setDof('U', step(scheme.TDisc, initial));
            state.setHistory('time', 0.1);
            result = scheme.run(state, resume=true);
            testCase.verifySameHandle(result, state);
            testCase.verifyEqual(result.Dofs.U, initial/(1.1^3), 'AbsTol', 2e-13);
            testCase.verifyEqual(result.History.time, 0.1);
            testCase.verifyEqual(scheme.TDisc.Timeline.Now, 0.3);
            testCase.verifyEqual(scheme.TDisc.Timeline.Count, 4);
        end

        function completedStateAndOdeHistoryRestartTogether(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            file = fullfile(fixture.Folder, 'checkpoint.mat');
            config = struct('tFinal', 0.5, 'dt', 0.1);
            scheme = mdemat.solver.hyperbolic.UwdgScheme(config=config);
            space = mdemat.common.fem.FiniteElementSpace( ...
                mdemat.data.mesh.UniformGrid([0, 1], 2), ...
                mdemat.common.fem.BH1OrthotopeElement.modal(1, 2));
            state = mdemat.common.fem.SpatialState(space);
            state.setDof('U', (1:space.NGlobalDofs).');
            state.Dofs.U = step(scheme.TDisc, state.Dofs.U);
            state.setHistory('time', scheme.TDisc.Timeline.Now);
            scheme.saveCheckpoint(file, state);
            resumed = mdemat.solver.hyperbolic.UwdgScheme(config=config);
            restored = resumed.loadCkptState(file, mdemat.common.fem.SpatialState(space));
            testCase.verifyEqual(restored.Dofs, state.Dofs);
            testCase.verifyEqual(restored.History, state.History);
            testCase.verifyEqual(resumed.TDisc.U0, scheme.TDisc.U0);
            testCase.verifyEqual(resumed.TDisc.Timeline.Now, 0.1);
            uninterrupted = step(scheme.TDisc, state.Dofs.U);
            continued = step(resumed.TDisc, restored.Dofs.U);
            testCase.verifyEqual(continued, uninterrupted, 'AbsTol', 2e-13);
            testCase.verifyEqual(continued, (1:space.NGlobalDofs).'/(1.1^2), 'AbsTol', 2e-13);
            % Existing plotting workflows may load only the state from a checkpoint.
            checkpoint = mdemat.internal.execution.readCheckpoint(file);
            stateOnly = mdemat.common.fem.SpatialState.fromRecord(checkpoint.state);
            testCase.verifyEqual(stateOnly.Dofs, state.Dofs);
            testCase.verifyError(@() scheme.saveCheckpoint(file, state), 'mdemat:io:FileExists');
        end

        function failedRestoreDoesNotMutateExistingWorkspace(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            file = fullfile(fixture.Folder, 'checkpoint.json');
            config = struct('tFinal', 1, 'dt', 0.1);
            scheme = mdemat.solver.hyperbolic.UwdgScheme(config=config);
            space = mdemat.common.fem.FiniteElementSpace( ...
                mdemat.data.mesh.UniformGrid([0, 1], 2), ...
                mdemat.common.fem.BH1OrthotopeElement.modal(1, 2));
            state = mdemat.common.fem.SpatialState(space);
            state.setDof('U', ones(space.NGlobalDofs, 1));
            scheme.saveCheckpoint(file, state);
            original = scheme.TDisc;
            record = mdemat.tools.io.readRecord(file);
            record.state.dofs.U = zeros(3, 1);
            mdemat.tools.io.writeRecord(file, record, overwrite=true);
            testCase.verifyError(@() scheme.loadCkptState(file, state), ...
                'mdemat:persistence:InvalidState');
            testCase.verifySameHandle(scheme.TDisc, original);
            testCase.verifyEqual(state.Dofs.U, ones(space.NGlobalDofs, 1));
            record.state = state.snapshot();
            record.configuration.values.dt = 0.2;
            mdemat.tools.io.writeRecord(file, record, overwrite=true);
            testCase.verifyError(@() scheme.loadCkptState(file, state), ...
                'mdemat:persistence:InvalidCheckpoint');
            testCase.verifySameHandle(scheme.TDisc, original);
        end

        function modelDescriptionsExcludeObserversAndNormalizeOldFiles(testCase)
            config = struct('dt',0.1,'tFinal',1,'ic',@sin,'components',struct('U',1), ...
                'exact',struct('U',@cos),'titlePrefix','first');
            current = mdemat.internal.execution.checkpointConfiguration(config);
            changed = config;
            changed.components.U = [2,4];
            changed.exact.U = @tan;
            changed.titlePrefix = 'other';
            testCase.verifyEqual(mdemat.internal.execution.checkpointConfiguration(changed),current);
            old = struct('values',rmfield(config,'ic'),'callbackNames',{{'ic'}});
            % Actual old persisted files cannot contain the nested exact callback;
            % an ordinary numeric reference is a legal v1 value.
            old.values.exact = struct('U',0);
            normalized = mdemat.internal.execution.normalizeCheckpointConfiguration(old);
            testCase.verifyEqual(normalized,current);
            changed.dt = 0.2;
            testCase.verifyNotEqual(mdemat.internal.execution.checkpointConfiguration(changed),current);
            nested = struct('ic',struct('U',@sin,'Q',0),'dt',0.1);
            identity = mdemat.internal.execution.checkpointConfiguration(nested);
            testCase.verifyEqual(identity.callbackNames,{'ic.U'});
            testCase.verifyEqual(identity.values.ic.U,[]);
            testCase.verifyEqual(identity.values.ic.Q,0);
        end
    end
end

function value = step(integrator, value)
integrator.setTimeStep(dt=0.1);
integrator.update(value);
value = integrator.step(L=-speye(numel(value)));
integrator.advance();
end
