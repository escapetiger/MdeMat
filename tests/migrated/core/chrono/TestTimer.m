classdef TestTimer < matlab.unittest.TestCase
    methods (Test)
        function namesAreCompleteAndAnonymousNamesCannotCollide(test)
            timer = mdemat.tools.analyze.Timer(names={'a','c';'b','d'}, verbose=false);
            test.verifyNumElements(timer.Records, 4);
            test.verifyEqual([timer.Records.Name], ["a","b","c","d"]);
            test.verifyError(@() mdemat.tools.analyze.Timer(names={'same','same'}), ...
                'mdemat:tools:timing:Timer:DuplicateName');
            test.verifyError(@() mdemat.tools.analyze.Timer(names={''}), ...
                'mdemat:tools:timing:Timer:InvalidName');
            timer = mdemat.tools.analyze.Timer(verbose=false);
            timer.start(record="Record_2"); timer.stop();
            timer.start(); timer.stop();
            test.verifyNumElements(unique([timer.Records.Name]), 2);
            test.verifyFalse(any([timer.Records.IsActive]));
        end

        function exclusiveIntervalsAccumulateAndSnapshotsRemainValues(test)
            timer = mdemat.tools.analyze.Timer(verbose=false);
            [~, interval] = timer.start(record='solve');
            active = timer.find('solve');
            test.verifyTrue(active.IsActive);
            test.verifyEqual(active.Duration, 0);
            test.verifyError(@() timer.start(), ...
                'mdemat:tools:timing:Timer:recordAlreadyActive');
            timer.stop(interval=interval);
            first = timer.find('solve');
            timer.start(record='solve'); timer.stop(record='solve');
            second = timer.find('solve');
            test.verifyGreaterThanOrEqual(second.Duration, first.Duration);
            test.verifyTrue(active.IsActive);
            first.Duration = -1;
            test.verifyGreaterThanOrEqual(timer.find('solve').Duration, 0);
            test.verifyEmpty(timer.find('missing'));
            test.verifyError(@() timer.stop(), 'mdemat:tools:timing:Timer:NoActiveRecord');
        end

        function missingAndInactiveNamesHaveSpecificErrors(test)
            timer = mdemat.tools.analyze.Timer(names={'ready'}, verbose=false);
            test.verifyError(@() timer.stop(record='missing'), ...
                'mdemat:tools:timing:Timer:RecordNotFound');
            test.verifyError(@() timer.stop(record='ready'), ...
                'mdemat:tools:timing:Timer:RecordNotActive');
            [~, first] = timer.start(record='ready');
            timer.stop();
            [~, second] = timer.start(record='ready');
            test.verifyError(@() timer.stop(interval=first), ...
                'mdemat:tools:timing:Timer:WrongInterval');
            timer.stop(interval=second);
        end

        function cleanupCancelsOnlyItsOwnFailedInterval(test)
            timer = mdemat.tools.analyze.Timer(verbose=false);
            test.verifyError(@() failDuringTiming(timer), 'foundation:WorkFailed');
            test.verifyFalse(any([timer.Records.IsActive]));
            test.verifyEqual(timer.find('failed').Duration, 0);
            [~, stale] = timer.start(record='first');
            timer.stop();
            [~, current] = timer.start(record='second');
            timer.cancel(stale);
            test.verifyTrue(timer.find('second').IsActive);
            timer.stop(interval=current);
            timer.cancel(current);
            test.verifyFalse(any([timer.Records.IsActive]));
        end

        function resetRetainsNamesAndHasAnExplicitActivePolicy(test)
            timer = mdemat.tools.analyze.Timer(verbose=false);
            timer.start(record='solve'); timer.stop();
            timer.reset();
            test.verifyEqual(timer.find('solve').Duration, 0);
            timer.start(record='solve');
            test.verifyWarning(@() timer.reset(), ...
                'mdemat:tools:timing:Timer:ResetWhileActive');
            test.verifyFalse(timer.find('solve').IsActive);
            test.verifyEqual(timer.find('solve').Duration, 0);
        end

        function reportsUseWallClockUnitsAndExcludeActiveWork(test)
            timer = mdemat.tools.analyze.Timer(verbose=true);
            output = evalc('timer.start(record="solve", msg="begin"); timer.stop(msg="end");');
            test.verifySubstring(output, 'begin');
            test.verifySubstring(output, 'end');
            report = timer.snapshot();
            test.verifyEqual(report.Total, sum([report.Records.Duration]));
            test.verifyEqual(report.Clock, "wall");
            [~, interval] = timer.start(record='active');
            active = timer.snapshot();
            test.verifyTrue(active.IsActive);
            test.verifyEqual(active.Total, report.Total);
            timer.cancel(interval);
            output = evalc('timer.report(); timer.report(perc=false);');
            test.verifySubstring(output, 'Wall-clock');
            test.verifyFalse(contains(output, 'CPU'));
            empty = mdemat.tools.analyze.Timer(verbose=false);
            output = evalc('empty.report();');
            test.verifyFalse(contains(output, 'NaN'));
        end

        function reportFilesAndBooleanInputsAreValidated(test)
            fixture = test.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            file = fullfile(fixture.Folder, 'timing.txt');
            stream = fopen(file, 'w');
            guard = onCleanup(@() fclose(stream)); %#ok<NASGU>
            timer = mdemat.tools.analyze.Timer(verbose=0);
            timer.start(); timer.stop();
            timer.report(file=stream);
            test.verifyError(@() timer.report(file=-1), ...
                'mdemat:tools:timing:Timer:InvalidFile');
            test.verifyError(@() timer.report(perc=[true,false]), ...
                'mdemat:tools:timing:Timer:InvalidBoolean');
            test.verifyError(@() mdemat.tools.analyze.Timer(verbose=2), ...
                'mdemat:tools:timing:Timer:InvalidBoolean');
        end
    end
end

function failDuringTiming(timer)
[~, interval] = timer.start(record='failed');
guard = onCleanup(@() timer.cancel(interval)); %#ok<NASGU>
error('foundation:WorkFailed', 'Original failure survives cleanup.');
end
