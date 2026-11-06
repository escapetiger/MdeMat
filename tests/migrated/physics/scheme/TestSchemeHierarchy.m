classdef TestSchemeHierarchy < matlab.unittest.TestCase
    methods (Test)
        function schemesComposeValueConfigAndExposeIsolatedSnapshots(test)
            scheme = mdemat.internal.execution.Scheme();
            test.verifyClass(scheme.Configuration, 'mdemat.tools.environment.Config');
            test.verifyEqual(scheme.getConfig('verbose'), 0);
            values = scheme.Config;
            values.verbose = 9;
            test.verifyEqual(scheme.Config.verbose, 0);
            separate = scheme.Configuration.set('verbose', 1);
            test.verifyEqual(separate.get('verbose'), 1);
            test.verifyEqual(scheme.Config.verbose, 0);
        end

        function unknownKeysRequireExplicitConfigDefinitions(test)
            test.verifyError(@() mdemat.internal.execution.Scheme( ...
                config=struct('customParam', 42)), 'mdemat:tools:environment:Config:UnknownKey');
            config = mdemat.tools.environment.Config();
            config = config.define('customParam', 42, @isnumeric);
            scheme = mdemat.internal.execution.Scheme(config=config);
            test.verifyEqual(scheme.getConfig('customParam'), 42);
            bad = config.define('verbose', -1, @(value) true);
            test.verifyError(@() mdemat.internal.execution.Scheme(config=bad), ...
                'mdemat:tools:environment:Config:InvalidValue');
            reserved = config.define('bcType', 'forged');
            test.verifyError(@() mdemat.internal.execution.Scheme(config=reserved), ...
                'mdemat:execution:Scheme:ReservedKey');
        end

        function odeDefaultsAndConstructedMethodAgree(test)
            scheme = mdemat.internal.execution.OdeScheme();
            test.verifyEqual(scheme.Config.tOdeIntName, {'be'});
            test.verifyEqual(scheme.Config.tFinal, 0);
            test.verifyEmpty(scheme.Config.tSplitName);
            test.verifyClass(scheme.TDisc, 'mdemat.common.ode.BeIntegrator');
            config = struct('tFinal', 2, 'tOdeIntName', 'sdirk3');
            configured = mdemat.internal.execution.OdeScheme(config=config);
            test.verifyClass(configured.TDisc, 'mdemat.common.ode.Sdirk3Integrator');
            test.verifyEqual(configured.TDisc.Timeline.Final, 2);
            split = mdemat.internal.execution.OdeScheme(config=struct( ...
                'tFinal', 1, 'tOdeIntName', {{'fe','be'}}, 'tSplitName', 'strang'));
            test.verifyClass(split.TDisc, 'mdemat.common.ode.StrangIntegrator');
        end

        function frozenNumericsCannotDisagreeWithExistingIntegrator(test)
            scheme = mdemat.internal.execution.OdeScheme(config=struct('tFinal', 1));
            before = scheme.Config;
            test.verifyError(@() scheme.load(config=struct( ...
                'tFinal', 2, 'tOdeIntName', 'fe')), ...
                'mdemat:execution:Scheme:FrozenConfiguration');
            test.verifyError(@() scheme.setConfig('tFinal', 2), ...
                'mdemat:execution:Scheme:FrozenConfiguration');
            test.verifyEqual(scheme.Config, before);
            test.verifyClass(scheme.TDisc, 'mdemat.common.ode.BeIntegrator');
            test.verifyEqual(scheme.TDisc.Timeline.Final, 1);
            test.verifyError(@() scheme.loadFromStruct(struct('verbose', 2, 'tFinal', -1)), ...
                'mdemat:execution:Scheme:FrozenConfiguration');
            test.verifyEqual(scheme.Config, before);
        end

        function fileAndExplicitPatchAreMergedBeforeValidation(test)
            fixture = test.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            file = fullfile(fixture.Folder, 'configuration.txt');
            stream = fopen(file, 'w', 'n', 'UTF-8');
            fprintf(stream, 'OPTION=verbose=1\nOPTION=tFinal=-1\nOPTION=tOdeIntName=bdf2\n');
            fclose(stream);
            scheme = mdemat.internal.execution.OdeScheme(file=file, ...
                config=struct('verbose', 2, 'tFinal', 5));
            test.verifyEqual(scheme.getConfig('verbose'), 2);
            test.verifyEqual(scheme.TDisc.Timeline.Final, 5);
            test.verifyClass(scheme.TDisc, 'mdemat.common.ode.Bdf2Integrator');
        end

        function invalidValuesAndCrossKeyCombinationsFailBeforeConstruction(test)
            test.verifyError(@() mdemat.internal.execution.OdeScheme( ...
                config=struct('tFinal', -1)), 'mdemat:tools:environment:Config:InvalidValue');
            test.verifyError(@() mdemat.internal.execution.OdeScheme( ...
                config=struct('tOdeIntName', 'invalid')), ...
                'mdemat:tools:environment:Config:InvalidValue');
            test.verifyError(@() mdemat.internal.execution.OdeScheme( ...
                config=struct('tFinal', Inf)), 'mdemat:tools:environment:Config:InvalidValue');
            test.verifyError(@() mdemat.internal.execution.OdeScheme( ...
                config=struct('tSplitName', 'strang')), ...
                'mdemat:execution:OdeScheme:InvalidSplitting');
            test.verifyError(@() mdemat.internal.execution.OdeScheme( ...
                config=struct('tOdeIntName', {{'fe','be'}})), ...
                'mdemat:execution:OdeScheme:InvalidSplitting');
        end

        function predicatesDoNotCaptureTheOriginalScheme(test)
            scheme = mdemat.internal.execution.OdeScheme();
            [~, predicate] = scheme.Configuration.get('tOdeIntName');
            delete(scheme);
            test.verifyTrue(predicate('fe'));
            test.verifyFalse(predicate('invalid'));
        end

        function applicationPredicatesDoNotRetainTheScheme(test)
            schemes = {mdemat.application.plasma.DghsScheme(), ...
                mdemat.application.plasma.ApDghsScheme()};
            for index = 1:numel(schemes)
                for key = {'tOdeIntName','tSplitName'}
                    [~, predicate] = schemes{index}.Configuration.get(key{1});
                    test.verifyFalse(capturesScheme(functions(predicate)));
                end
            end
        end

        function deferredInputsAreCheckedByTheOperationThatUsesThem(test)
            location = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
            test.applyFixture(matlab.unittest.fixtures.PathFixture(location));
            partial = DeferredSchemeProbe();
            partial.validateConfig('verbose');
            test.verifyEmpty(partial.Config.initial);
            test.verifyError(@() partial.initialize(3), ...
                'mdemat:tools:environment:Config:InvalidValue');
            ready = DeferredSchemeProbe(config=struct('initial', @(x) x));
            test.verifyEqual(ready.initialize(3), 3);
        end

        function observationChangesDoNotChangeNumericalBindings(test)
            scheme = mdemat.internal.execution.OdeScheme(config=struct('tFinal', 1));
            before = scheme.Config;
            scheme.setObservation('plot', struct('Components', struct('U', [2,4]), ...
                'PointsPerDimension', 4));
            scheme.setObservation('plot', struct('Title', 'explicit'));
            test.verifyEqual(scheme.Observation.Plot.Components, struct('U', [2,4]));
            test.verifyEqual(scheme.Observation.Plot.PointsPerDimension, 4);
            test.verifyError(@() scheme.setObservation('plot', struct('Density', 4)), ...
                'mdemat:execution:Scheme:UnknownObservationKey');
            test.verifyEqual(scheme.Config, before);
            timer = mdemat.tools.analyze.Timer(verbose=false);
            scheme.setObserver('timer', timer);
            test.verifySameHandle(scheme.getObserver('timer'), timer);
            test.verifyFalse(isfield(scheme.Config, 'timer'));
            test.verifyError(@() scheme.setObserver('analyzer', []), ...
                'mdemat:execution:Scheme:UnknownObserver');
            test.verifyEqual(scheme.TDisc.Timeline.Final, 1);
        end

        function derivedWorkspaceCannotReplaceFrozenInputs(test)
            scheme = mdemat.internal.execution.Scheme();
            scheme.setDerived('xId', 'X2');
            test.verifyEqual(scheme.getConfig('xId'), 'X2');
            test.verifyTrue(scheme.hasConfig('xId'));
            test.verifyFalse(scheme.Configuration.has('xId'));
            test.verifyError(@() scheme.setDerived('verbose', 9), ...
                'mdemat:execution:Scheme:ReservedKey');
            test.verifyEqual(scheme.Config.verbose, 0);
        end

        function studyVocabularyIsExplicitAndDoesNotWeakenOwnerRules(test)
            location = fullfile(fileparts(mfilename('fullpath')), 'fixtures');
            test.applyFixture(matlab.unittest.fixtures.PathFixture(location));
            config = struct('nDims', 2, 'nx', [4,6], 'xBBox', [0,1,-1,1], ...
                'ampl', [0,1], 'omega', [pi,2*pi], 'tau', Inf, ...
                'vBasisType', 'hermite', 'useSplitting', 1);
            scheme = StudySchemeProbe(config=config);
            test.verifyEqual(scheme.Config.nx, [4,6]);
            test.verifyEqual(scheme.Config.L, 2);
            test.verifyEqual(scheme.Config.tau, Inf);
            test.verifyError(@() StudySchemeProbe(config=struct('L',3)), ...
                'mdemat:tools:environment:Config:InvalidValue');
            test.verifyError(@() StudySchemeProbe(config=struct('nx',[4,2.5])), ...
                'mdemat:tools:environment:Config:InvalidValue');
            test.verifyError(@() StudySchemeProbe(config=struct('xBBox',[1,0])), ...
                'mdemat:tools:environment:Config:InvalidValue');
            test.verifyError(@() StudySchemeProbe(config=struct('nDimz',2)), ...
                'mdemat:tools:environment:Config:UnknownKey');
        end
    end
end

function tf = capturesScheme(value)
if isa(value, 'mdemat.internal.execution.Scheme')
    tf = true;
elseif iscell(value)
    tf = any(cellfun(@capturesScheme, value), 'all');
elseif isstruct(value)
    tf = any(cellfun(@capturesScheme, struct2cell(value)), 'all');
else
    tf = false;
end
end
