classdef TestTimer < matlab.unittest.TestCase
    % TESTTIMER Unit tests for Timer and Record classes.
    %
    %   TestTimer provides test coverage for the Timer and Record classes
    %   including timing operations, error handling, and reporting.
    %
    % Examples:
    %   % Run all tests
    %   results = runtests(core.chrono.TestTimer);
    %
    % See Also:
    %   core.chrono.Timer, core.chrono.Record

    properties
        testRecordNames = {'Setup', 'Calculation', 'Cleanup'} % Names for test records
        timer % Timer object for testing
        pauseDuration = 0.01 % Duration for timing tests
    end

    methods (TestMethodSetup)
        function createTimer(testCase)
            % CREATETIMER Initialize timer for each test method.

            testCase.timer = core.chrono.Timer(testCase.testRecordNames);
            testCase.timer.verbose = false;
        end
    end

    methods (Test)
        function testConstructor(testCase)
            % TESTCONSTRUCTOR Test Timer constructor functionality.

            testCase.verifyNumElements(testCase.timer.records, numel(testCase.testRecordNames));

            for i = 1:numel(testCase.testRecordNames)
                testCase.verifyEqual(testCase.timer.records(i).name, testCase.testRecordNames{i});
            end

            testCase.verifyError(@() core.chrono.Timer(), ...
                'core:chrono:Timer:InvalidInput');
            testCase.verifyError(@() core.chrono.Timer('not_a_cell'), ...
                'core:chrono:Timer:InvalidInput');
        end

        function testStartStopByIndex(testCase)
            % TESTSTARTstopbyindex Test start/stop operations using record indices.

            testCase.timer.start(1);
            pause(testCase.pauseDuration);
            testCase.timer.stop(1);

            testCase.verifyGreaterThan(testCase.timer.records(1).duration, 0);

            testCase.verifyError(@() testCase.timer.start(0), ...
                'core:chrono:Timer:InvalidIndex');
            testCase.verifyError(@() testCase.timer.start(numel(testCase.testRecordNames)+1), ...
                'core:chrono:Timer:InvalidIndex');
            testCase.verifyError(@() testCase.timer.stop(0), ...
                'core:chrono:Timer:InvalidIndex');
            testCase.verifyError(@() testCase.timer.stop(numel(testCase.testRecordNames)+1), ...
                'core:chrono:Timer:InvalidIndex');
            testCase.verifyError(@() testCase.timer.stop(2), ...
                'core:chrono:Timer:RecordNotActive');
        end

        function testStartStopByName(testCase)
            % TESTSTARTOPBYNAME Test start/stop operations using record names.

            testCase.timer.startByName('Calculation');
            pause(testCase.pauseDuration);
            testCase.timer.stopByName('Calculation');

            testCase.verifyGreaterThan(testCase.timer.records(2).duration, 0);

            testCase.verifyError(@() testCase.timer.startByName('NonexistentRecord'), ...
                'core:chrono:Timer:RecordNotFound');
            testCase.verifyError(@() testCase.timer.stopByName('NonexistentRecord'), ...
                'core:chrono:Timer:RecordNotFound');

            testCase.verifyError(@() testCase.timer.stopByName('Cleanup'), ...
                'core:chrono:Timer:RecordNotActive');
        end

        function testMultipleTimings(testCase)
            % TESTMULTIPLETIMINGS Test accumulation of multiple timing sessions.

            testCase.timer.start(1);
            pause(testCase.pauseDuration);
            testCase.timer.stop(1);

            initialDuration = testCase.timer.records(1).duration;

            testCase.timer.start(1);
            pause(testCase.pauseDuration);
            testCase.timer.stop(1);

            finalDuration = testCase.timer.records(1).duration;

            testCase.verifyGreaterThan(finalDuration, initialDuration);
        end

        function testMutualExclusion(testCase)
            % TESTMUTUALEXCLUSION Test that only one record can be active at a time.

            testCase.timer.start(1);
            testCase.verifyError(@() testCase.timer.start(2), ...
                'core:chrono:Timer:recordAlreadyActive');
            testCase.timer.stop(1);

            testCase.timer.start(2);
            testCase.timer.stop(2);
        end

        function testReset(testCase)
            % TESTRESET Test reset functionality for all records.

            testCase.timer.start(1);
            pause(testCase.pauseDuration);
            testCase.timer.stop(1);

            testCase.verifyGreaterThan(testCase.timer.records(1).duration, 0);

            testCase.timer.reset();

            for i = 1:numel(testCase.timer.records)
                testCase.verifyEqual(testCase.timer.records(i).duration, 0);
            end
        end

        function testRecordStatus(testCase)
            % TESTRECORDSTATUS Test record active status detection.

            % Test no active records initially
            testCase.verifyFalse(testCase.timer.records(1).isActive());
            testCase.verifyFalse(testCase.timer.records(2).isActive());

            % Test record becomes active when started
            testCase.timer.start(2);
            testCase.verifyTrue(testCase.timer.records(2).isActive());

            % Test record becomes inactive when stopped
            testCase.timer.stop(2);
            testCase.verifyFalse(testCase.timer.records(2).isActive());
        end

        function testTimingData(testCase)
            % TESTTIMINGDATA Test timing data structure generation.

            testCase.timer.start(1);
            pause(testCase.pauseDuration);
            testCase.timer.stop(1);

            testCase.timer.start(2);
            pause(testCase.pauseDuration*2);
            testCase.timer.stop(2);

            % Test basic timing functionality
            testCase.verifyGreaterThan(testCase.timer.records(1).duration, 0);
            testCase.verifyGreaterThan(testCase.timer.records(2).duration, 0);
            testCase.verifyGreaterThan(testCase.timer.records(2).duration, testCase.timer.records(1).duration);
        end

        function testReport(testCase)
            % TESTREPORT Test report generation functionality.

            testCase.timer.start(1);
            pause(testCase.pauseDuration);
            testCase.timer.stop(1);

            testCase.timer.start(2);
            pause(testCase.pauseDuration*2);
            testCase.timer.stop(2);

            tempFile = [tempname, '.txt'];
            fileId = fopen(tempFile, 'w');

            testCase.timer.report(true, fileId);
            testCase.timer.report(false, fileId);

            fclose(fileId);

            testCase.verifyTrue(isfile(tempFile));

            if isfile(tempFile)
                delete(tempFile);
            end
        end

        function testVerboseOutput(testCase)
            % TESTVERBOSEOUTPUT Test verbose message output control.

            testCase.timer.verbose = true;

            output = evalc('testCase.timer.start(1, ''Test message'')');
            testCase.timer.stop(1);

            testCase.verifySubstring(output, 'Test message');

            testCase.timer.verbose = false;

            output = evalc('testCase.timer.start(1, ''Test message'')');
            testCase.timer.stop(1);

            testCase.verifyEmpty(output);
        end

        function testComplexTimingSequence(testCase)
            % TESTCOMPLEXTIMINGSEQUENCE Test complex timing operations.

            % First record timing
            testCase.timer.start(1);
            pause(testCase.pauseDuration);
            testCase.timer.stop(1);

            % Second record timing (longer)
            testCase.timer.start(2);
            pause(testCase.pauseDuration*2);
            testCase.timer.stop(2);

            % Third record timing
            testCase.timer.start(3);
            pause(testCase.pauseDuration);
            testCase.timer.stop(3);

            % First record again (accumulates)
            testCase.timer.start(1);
            pause(testCase.pauseDuration);
            testCase.timer.stop(1);

            % Verify timing relationships
            testCase.verifyGreaterThan(testCase.timer.records(1).duration, testCase.timer.records(3).duration);
            testCase.verifyGreaterThan(testCase.timer.records(2).duration, testCase.timer.records(1).duration);
        end

        function testZeroDurations(testCase)
            % TESTZERODURATIONS Test behavior with zero durations.

            % Verify all records start with zero duration
            for i = 1:numel(testCase.timer.records)
                testCase.verifyEqual(testCase.timer.records(i).duration, 0);
            end
        end

        function testMissingIndex(testCase)
            % TESTMISSINGINDEX Test error handling for missing index arguments.

            testCase.verifyError(@() testCase.timer.start(), ...
                'core:chrono:Timer:MissingIndex');
            testCase.verifyError(@() testCase.timer.stop(), ...
                'core:chrono:Timer:MissingIndex');
        end

        function testMissingName(testCase)
            % TESTMISSINGNAME Test error handling for missing name arguments.

            testCase.verifyError(@() testCase.timer.startByName(), ...
                'core:chrono:Timer:MissingName');
            testCase.verifyError(@() testCase.timer.stopByName(), ...
                'core:chrono:Timer:MissingName');
        end

        function testIndividualRecord(testCase)
            % TESTINDIVIDUALRECORD Test individual Record class functionality.

            record = core.chrono.Record('TestRecord');

            testCase.verifyEqual(record.name, 'TestRecord');
            testCase.verifyEqual(record.duration, 0);
            testCase.verifyFalse(record.isActive());

            record.startTime = tic;
            testCase.verifyTrue(record.isActive());

            pause(testCase.pauseDuration);

            record.duration = record.duration + toc(record.startTime);
            record.startTime = [];

            testCase.verifyGreaterThan(record.duration, 0);
            testCase.verifyFalse(record.isActive());

            record.reset();
            testCase.verifyEqual(record.duration, 0);
        end
    end
end