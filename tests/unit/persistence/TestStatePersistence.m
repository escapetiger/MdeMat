classdef TestStatePersistence < matlab.unittest.TestCase
    properties
        Folder
    end
    methods (TestMethodSetup)
        function temporaryDestination(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.Folder = fixture.Folder;
        end
    end
    methods (Test)
        function everyStateKindRoundtripsWithItsActualRepresentation(testCase)
            x = spatialSpace();
            v = mdemat.common.fem.SpectralSpace( ...
                mdemat.common.fem.L2OrthotopeElement.hermite(1, 4, T=0.7));
            macro = mdemat.common.fem.SpectralSpace(mdemat.common.fem.L2SphereElement.modal(1, 2));
            micro = mdemat.common.fem.SpectralSpace(mdemat.common.fem.L2SphereElement.modal(1, 3));
            velocity = mdemat.common.fem.SumSpace(macro, micro, ...
                lhsWeight=[1; 0.1], rhsWeight=[0.01; 0.02; 0.03]);
            states = {mdemat.data.variable.State(), mdemat.common.fem.SpatialState(x), ...
                mdemat.common.fem.KineticState(x, v), ...
                mdemat.application.plasma.HermiteState(x, v), ...
                mdemat.application.plasma.ApHermiteState(x, v), ...
                mdemat.application.radiation.MacroMicroState(x, velocity)};
            loaders = {@mdemat.data.variable.State.load, @mdemat.common.fem.SpatialState.load, ...
                @mdemat.common.fem.KineticState.load, @mdemat.application.plasma.HermiteState.load, ...
                @mdemat.application.plasma.ApHermiteState.load, ...
                @mdemat.application.radiation.MacroMicroState.load};
            for index = 1:numel(states)
                state = states{index};
                state.setDof('sample', reshape(1:2*x.NGlobalDofs, [], 2) + 1i/7);
                state.setCoefficient('constant', 0.7);
                state.setHistory('time', [0; 0.1]);
                state.setHistory('mass', [1; 1]);
                if isa(state, 'mdemat.application.radiation.MacroMicroState')
                    state.setDecomposition('partial');
                    state.Cache = struct('ignoredCallback', @sin);
                    state.Status = state.IsAllComputed;
                end
                for extension = {'.mat', '.json'}
                    file = fullfile(testCase.Folder, sprintf('state-%d%s', index, extension{1}));
                    state.save(file);
                    rebuilt = loaders{index}(file);
                    testCase.verifyClass(rebuilt, class(state));
                    testCase.verifyEqual(rebuilt.Dofs, state.Dofs);
                    testCase.verifyEqual(rebuilt.Coefs, state.Coefs);
                    testCase.verifyEqual(rebuilt.History, state.History);
                    testCase.verifyEqual(rebuilt.snapshot(), state.snapshot());
                    if isa(rebuilt, 'mdemat.common.fem.SpatialState')
                        testCase.verifyNotSameHandle(rebuilt.XDisc, state.XDisc);
                        sample = rebuilt.XDisc.eval(0.1, rebuilt.Dofs.sample);
                        testCase.verifyEqual(sample, state.XDisc.eval(0.1, state.Dofs.sample));
                    end
                    if isa(rebuilt, 'mdemat.application.radiation.MacroMicroState')
                        testCase.verifyEmpty(rebuilt.Cache);
                        testCase.verifyEqual(rebuilt.Status, uint8(0));
                        testCase.verifyEqual(rebuilt.Decomposition, state.Decomposition);
                        testCase.verifyEqual(rebuilt.VDisc.Weights, state.VDisc.Weights);
                    end
                end
            end
        end

        function callbacksRequireExplicitOmissionAndRebinding(testCase)
            state = mdemat.data.variable.State();
            state.setCoefficient('source', @(x) 3*x.^2);
            file = fullfile(testCase.Folder, 'callbacks.mat');
            testCase.verifyError(@() state.save(file), 'mdemat:persistence:UnsupportedCallback');
            testCase.verifyFalse(isfile(file));
            state.save(file, omitCallbacks=true);
            testCase.verifyError(@() mdemat.data.variable.State.load(file), ...
                'mdemat:persistence:CallbacksRequired');
            rebuilt = mdemat.data.variable.State.load(file, ...
                callbacks=struct('source', @(x) 3*x.^2));
            testCase.verifyEqual(rebuilt.Coefs.source([1, 2]), [3, 12]);
            testCase.verifyError(@() mdemat.data.variable.State.load(file, ...
                callbacks=struct('source', @sin, 'extra', @cos)), ...
                'mdemat:persistence:CallbacksRequired');
        end

        function overwriteAndInvalidSnapshotLeaveExistingDataIntact(testCase)
            state = mdemat.common.fem.SpatialState(spatialSpace());
            state.setDof('u', ones(state.XDisc.NGlobalDofs, 1));
            file = fullfile(testCase.Folder, 'state.json');
            state.save(file);
            original = fileread(file);
            testCase.verifyError(@() state.save(file), 'mdemat:io:FileExists');
            state.XDisc.Element.Approximator.Basis = containers.Map();
            testCase.verifyError(@() state.save(file, overwrite=true), ...
                'mdemat:persistence:UnsupportedRepresentation');
            testCase.verifyEqual(fileread(file), original);
            rebuilt = mdemat.common.fem.SpatialState.load(file);
            testCase.verifyEqual(rebuilt.Dofs.u, ones(rebuilt.XDisc.NGlobalDofs, 1));
        end

        function invalidDofPackingAndWrongStateKindAreRejected(testCase)
            state = mdemat.common.fem.SpatialState(spatialSpace());
            record = state.snapshot(); record.dofs.u = ones(2, 1);
            testCase.verifyError(@() mdemat.common.fem.SpatialState.fromRecord(record), ...
                'mdemat:persistence:InvalidState');
            testCase.verifyError(@() mdemat.data.variable.State.fromRecord(state.snapshot()), ...
                'mdemat:persistence:InvalidState');
            record = state.snapshot(); record.spaces.x.type = 'arbitrary.class';
            testCase.verifyError(@() mdemat.common.fem.SpatialState.fromRecord(record), ...
                'mdemat:persistence:UnsupportedTag');
        end

        function opaqueLegacyMatFilesAreNotLoaded(testCase)
            file = fullfile(testCase.Folder, 'opaque.mat');
            state = struct('xDisc', spatialSpace(), 'dofs', struct()); %#ok<NASGU>
            save(file, 'state');
            testCase.verifyError(@() mdemat.common.fem.SpatialState.load(file), ...
                'mdemat:io:InvalidRecord');
        end
    end
end

function space = spatialSpace()
space = mdemat.common.fem.FiniteElementSpace( ...
    mdemat.data.mesh.UniformGrid([0, 1], 3), ...
    mdemat.common.fem.BH1OrthotopeElement.modal(1, 3));
end
