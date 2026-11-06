classdef TestConfig < matlab.unittest.TestCase
    % Independent configuration ownership, validation and failure contracts.
    methods (Test)
        function entriesRemainPairedAndUpdatesHaveValueSemantics(test)
            original = mdemat.tools.environment.Config();
            predicate = @(value) isnumeric(value) && isscalar(value) && value > 0;
            original = original.define('step', 1, predicate);
            changed = original.set('step', 2);
            original.set('step', 3);
            [value, storedPredicate] = original.get('step');
            test.verifyEqual(value, 1);
            test.verifyEqual(storedPredicate, predicate);
            test.verifyEqual(changed.get('step'), 2);
            test.verifyFalse(isa(original, 'handle'));
            test.verifyEqual(original.keys(), {'step'});
            test.verifyEqual(original.toStruct(), struct('step', 1));
        end

        function defineRequiresExplicitOverride(test)
            config = mdemat.tools.environment.Config();
            config = config.define('value', 1);
            test.verifyError(@() config.define('value', 2), ...
                'mdemat:tools:environment:Config:DuplicateKey');
            config = config.define('value', 2, @isnumeric, override=true);
            test.verifyEqual(config.get('value'), 2);
            test.verifyTrue(config.has('value'));
            test.verifyFalse(config.has('missing'));
            test.verifyError(@() config.set('missing', 1), ...
                'mdemat:tools:environment:Config:UnknownKey');
            test.verifyError(@() config.define('bad name', 1), ...
                'mdemat:tools:environment:Config:InvalidKey');
            test.verifyError(@() config.define('badValidator', 1, 42), ...
                'mdemat:tools:environment:Config:InvalidValidator');
        end

        function failedPatchDoesNotChangeAnEarlierKey(test)
            config = mdemat.tools.environment.Config();
            config = config.define('first', 1);
            config = config.define('positive', 2, @(value) isscalar(value) && value > 0);
            test.verifyError(@() config.update(struct('first', 99, 'positive', -1)), ...
                'mdemat:tools:environment:Config:InvalidValue');
            test.verifyEqual(config.toStruct(), struct('first', 1, 'positive', 2));
            test.verifyError(@() config.update(struct('first', 99, 'typo', 1)), ...
                'mdemat:tools:environment:Config:UnknownKey');
            test.verifyEqual(config.get('first'), 1);
        end

        function deferredDefaultsAreCheckedByOperation(test)
            config = mdemat.tools.environment.Config();
            config = config.define('initial', [], @(value) isa(value, 'function_handle'));
            config = config.define('step', [], @(value) isempty(value) || ...
                (isnumeric(value) && isscalar(value) && value > 0));
            config.validate('step');
            config.validate({});
            test.verifyError(@() config.validate('initial'), ...
                'mdemat:tools:environment:Config:InvalidValue');
            resolved = config.set('initial', @(x) x);
            resolved.validate();
            test.verifyEmpty(config.get('initial'));
        end

        function predicatesAreCalledOnceAndHaveOneExplicitProtocol(test)
            calls = 0;
            config = mdemat.tools.environment.Config();
            config = config.define('value', 0, @predicate);
            test.verifyEqual(calls, 0);
            config = config.set('value', 1);
            test.verifyEqual(calls, 1);
            config.validate();
            test.verifyEqual(calls, 2);
            bad = config.define('bad', 1, @(value) double(value));
            test.verifyError(@() bad.set('bad', 2), ...
                'mdemat:tools:environment:Config:InvalidValidatorResult');
            function tf = predicate(value)
                calls = calls + 1;
                tf = isscalar(value) && value >= 0;
            end
        end

        function originalValidatorExceptionRemainsAvailable(test)
            config = mdemat.tools.environment.Config();
            config = config.define('value', 1, @failingPredicate);
            try
                config.set('value', 2);
                test.assertFail('Expected the predicate to fail.');
            catch exception
                test.verifyEqual(exception.identifier, ...
                    'mdemat:tools:environment:Config:InvalidValue');
                test.assertNumElements(exception.cause, 1);
                test.verifyEqual(exception.cause{1}.identifier, 'foundation:PredicateFailure');
            end
        end

        function containedHandlesRemainSharedWithoutSharingTheMap(test)
            timer = mdemat.tools.analyze.Timer(verbose=false);
            first = mdemat.tools.environment.Config();
            first = first.define('observer', timer);
            second = first;
            observer = second.get('observer');
            observer.Verbose = true;
            test.verifyTrue(timer.Verbose);
            second = second.set('observer', []);
            test.verifySameHandle(first.get('observer'), timer);
            test.verifyEmpty(second.get('observer'));
        end
    end
end

function result = failingPredicate(~) %#ok<STOUT>
error('foundation:PredicateFailure', 'Original validator context.');
end
