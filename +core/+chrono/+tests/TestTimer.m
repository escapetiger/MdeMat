classdef TestTimer < matlab.unittest.TestCase
    % TESTTIMER Comprehensive test suite for Timer and Record classes.
    %
    %   TestTimer provides complete test coverage for the Timer and Record
    %   classes including timing operations, error handling, state
    %   management, and report generation. Tests validate both individual
    %   functionality and integration between classes.
    %
    % See also:
    %   core.chrono.Timer, core.chrono.Record

    properties
        Timer % Timer object for testing
        PauseDuration = 0.01 % Duration for timing tests
    end

    methods (TestMethodSetup)
        function createTimer(testCase)
            % CREATETIMER Initialize fresh timer instance for each test.

            testCase.Timer = core.chrono.Timer(verbose=false);
        end
    end

    methods (Test)
        function testConstructor(testCase)
            % TESTCONSTRUCTOR Validate Timer constructor behavior and error handling.

            % Test empty constructor
            emptyTimer = core.chrono.Timer();
            testCase.verifyNumElements(emptyTimer.Records, 0);
            testCase.verifyTrue(emptyTimer.Verbose); % Default is true
            
            % Test constructor with names
            namedTimer = core.chrono.Timer(names={'Setup', 'Calculation'});
            testCase.verifyNumElements(namedTimer.Records, 2);
            testCase.verifyEqual(namedTimer.Records(1).Name, "Setup");
            testCase.verifyEqual(namedTimer.Records(2).Name, "Calculation");
            
            % Test constructor with verbose control
            quietTimer = core.chrono.Timer(verbose=false);
            testCase.verifyFalse(quietTimer.Verbose);
            
            % Test constructor with both options
            fullTimer = core.chrono.Timer(names={'test'}, verbose=false);
            testCase.verifyNumElements(fullTimer.Records, 1);
            testCase.verifyFalse(fullTimer.Verbose);
        end

        function testStartStop(testCase)
            % TESTSTARTSTOP Validate basic start/stop operations.

            testCase.Timer.start();
            pause(testCase.PauseDuration);
            testCase.Timer.stop();

            testCase.verifyNumElements(testCase.Timer.Records, 1);
            testCase.verifyGreaterThan(testCase.Timer.Records(1).Duration, 0);
            testCase.verifyEqual(testCase.Timer.Records(1).Name, "Record_1");
            
            % Test error when no active record to stop
            testCase.verifyError(@() testCase.Timer.stop(), ...
                'core:chrono:Timer:NoActiveRecord');
        end

        function testStartStopWithMessage(testCase)
            % TESTSTARTSTOPWITHMESSAGE Validate start/stop operations with messages.

            verboseTimer = core.chrono.Timer(verbose=true);
            output = evalc('verboseTimer.start(msg="Starting test"); pause(testCase.PauseDuration); verboseTimer.stop(msg="Stopping test");');
            
            testCase.verifySubstring(output, 'Starting test');
            testCase.verifySubstring(output, 'Stopping test');
            testCase.verifyGreaterThan(verboseTimer.Records(1).Duration, 0);
        end

        function testMultipleRecords(testCase)
            % TESTMULTIPLERECORDS Validate creation of multiple timing records.

            % Create first record
            testCase.Timer.start();
            pause(testCase.PauseDuration);
            testCase.Timer.stop();

            % Create second record
            testCase.Timer.start();
            pause(testCase.PauseDuration * 2);
            testCase.Timer.stop();

            testCase.verifyNumElements(testCase.Timer.Records, 2);
            testCase.verifyGreaterThan(testCase.Timer.Records(1).Duration, 0);
            testCase.verifyGreaterThan(testCase.Timer.Records(2).Duration, testCase.Timer.Records(1).Duration);
        end

        function testMutualExclusion(testCase)
            % TESTMUTUALEXCLUSION Validate exclusive record activation constraint.

            testCase.Timer.start();
            testCase.verifyError(@() testCase.Timer.start(), ...
                'core:chrono:Timer:recordAlreadyActive');
            testCase.Timer.stop();

            % Should be able to start again after stopping
            testCase.Timer.start();
            testCase.Timer.stop();
        end

        function testReset(testCase)
            % TESTRESET Validate reset functionality for all timing records.

            testCase.Timer.start();
            pause(testCase.PauseDuration);
            testCase.Timer.stop();

            testCase.verifyGreaterThan(testCase.Timer.Records(1).Duration, 0);

            testCase.Timer.reset();

            for i = 1:numel(testCase.Timer.Records)
                testCase.verifyEqual(testCase.Timer.Records(i).Duration, 0);
            end
        end

        function testRecordStatus(testCase)
            % TESTRECORDSTATUS Validate record active status detection.

            % Test no active records initially (empty timer)
            testCase.verifyNumElements(testCase.Timer.Records, 0);

            % Test record becomes active when started
            testCase.Timer.start();
            testCase.verifyTrue(testCase.Timer.Records(1).IsActive);

            % Test record becomes inactive when stopped
            testCase.Timer.stop();
            testCase.verifyFalse(testCase.Timer.Records(1).IsActive);
        end

        function testTimingData(testCase)
            % TESTTIMINGDATA Validate timing data accuracy and relationships.

            testCase.Timer.start();
            pause(testCase.PauseDuration);
            testCase.Timer.stop();

            testCase.Timer.start();
            pause(testCase.PauseDuration * 2);
            testCase.Timer.stop();

            % Test basic timing functionality
            testCase.verifyGreaterThan(testCase.Timer.Records(1).Duration, 0);
            testCase.verifyGreaterThan(testCase.Timer.Records(2).Duration, 0);
            testCase.verifyGreaterThan(testCase.Timer.Records(2).Duration, testCase.Timer.Records(1).Duration);
        end

        function testReport(testCase)
            % TESTREPORT Validate report generation and file output.

            testCase.Timer.start();
            pause(testCase.PauseDuration);
            testCase.Timer.stop();

            testCase.Timer.start();
            pause(testCase.PauseDuration * 2);
            testCase.Timer.stop();

            tempFile = [tempname, '.txt'];
            fileId = fopen(tempFile, 'w');

            testCase.Timer.report(perc=true, file=fileId);
            testCase.Timer.report(perc=false, file=fileId);

            fclose(fileId);

            testCase.verifyTrue(isfile(tempFile));

            if isfile(tempFile)
                delete(tempFile);
            end
        end

        function testVerboseOutput(testCase)
            % TESTVERBOSEOUTPUT Validate verbose message output control.

            verboseTimer = core.chrono.Timer(verbose=true);

            output = evalc('verboseTimer.start(msg="Test message")');
            verboseTimer.stop();

            testCase.verifySubstring(output, 'Test message');
        end

        function testComplexTimingSequence(testCase)
            % TESTCOMPLEXTIMINGSEQUENCE Validate complex multi-record timing scenarios.

            % First record timing
            testCase.Timer.start();
            pause(testCase.PauseDuration);
            testCase.Timer.stop();

            % Second record timing (longer)
            testCase.Timer.start();
            pause(testCase.PauseDuration * 2);
            testCase.Timer.stop();

            % Third record timing
            testCase.Timer.start();
            pause(testCase.PauseDuration);
            testCase.Timer.stop();

            % Verify we have 3 records
            testCase.verifyNumElements(testCase.Timer.Records, 3);
            testCase.verifyGreaterThan(testCase.Timer.Records(2).Duration, testCase.Timer.Records(1).Duration);
            testCase.verifyGreaterThan(testCase.Timer.Records(2).Duration, testCase.Timer.Records(3).Duration);
        end

        function testZeroDurations(testCase)
            % TESTZERODURATIONS Validate initial state and zero duration handling.

            % Verify timer starts empty
            testCase.verifyNumElements(testCase.Timer.Records, 0);
            
            % Create a record and verify it starts with zero duration
            testCase.Timer.start();
            testCase.Timer.stop();
            testCase.verifyEqual(testCase.Timer.Records(1).Duration, 0, 'AbsTol', testCase.PauseDuration);
        end

        function testIndividualRecord(testCase)
            % TESTINDIVIDUALRECORD Validate individual Record class functionality.

            record = core.chrono.Record('TestRecord');

            testCase.verifyEqual(record.Name, "TestRecord");
            testCase.verifyEqual(record.Duration, 0);
            testCase.verifyFalse(record.IsActive);

            record.StartTime = tic;
            testCase.verifyTrue(record.IsActive);

            pause(testCase.PauseDuration);

            record.Duration = record.Duration + toc(record.StartTime);
            record.StartTime = [];

            testCase.verifyGreaterThan(record.Duration, 0);
            testCase.verifyFalse(record.IsActive);

            record.reset();
            testCase.verifyEqual(record.Duration, 0);
        end
        
        function testEmptyTimerOperations(testCase)
            % TESTEMPTYTIMEROPERATIONS Validate operations on empty timer.
            
            emptyTimer = core.chrono.Timer(verbose=false);
            testCase.verifyNumElements(emptyTimer.Records, 0);
            
            % Should be able to start timing (creates first record)
            emptyTimer.start();
            testCase.verifyNumElements(emptyTimer.Records, 1);
            testCase.verifyEqual(emptyTimer.Records(1).Name, "Record_1");
            
            emptyTimer.stop();
            testCase.verifyFalse(emptyTimer.Records(1).IsActive);
        end

        function testNamedRecords(testCase)
            % TESTNAMEDRECORDS Validate named record start/stop functionality.

            % Start timing for specific named record
            testCase.Timer.start(record="computation");
            pause(testCase.PauseDuration);
            testCase.Timer.stop(record="computation");

            testCase.verifyNumElements(testCase.Timer.Records, 1);
            testCase.verifyEqual(testCase.Timer.Records(1).Name, "computation");
            testCase.verifyGreaterThan(testCase.Timer.Records(1).Duration, 0);
            testCase.verifyFalse(testCase.Timer.Records(1).IsActive);
        end

        function testMultipleNamedRecords(testCase)
            % TESTMULTIPLENAMEDRECORDS Validate multiple named records functionality.

            % Create and time different named records
            testCase.Timer.start(record="setup");
            pause(testCase.PauseDuration);
            testCase.Timer.stop(record="setup");

            testCase.Timer.start(record="computation");
            pause(testCase.PauseDuration * 2);
            testCase.Timer.stop(record="computation");

            testCase.Timer.start(record="cleanup");
            pause(testCase.PauseDuration);
            testCase.Timer.stop(record="cleanup");

            testCase.verifyNumElements(testCase.Timer.Records, 3);
            testCase.verifyEqual(testCase.Timer.Records(1).Name, "setup");
            testCase.verifyEqual(testCase.Timer.Records(2).Name, "computation");
            testCase.verifyEqual(testCase.Timer.Records(3).Name, "cleanup");

            % Verify computation took longer than others
            testCase.verifyGreaterThan(testCase.Timer.Records(2).Duration, ...
                testCase.Timer.Records(1).Duration);
            testCase.verifyGreaterThan(testCase.Timer.Records(2).Duration, ...
                testCase.Timer.Records(3).Duration);
        end

        function testRecordAccumulation(testCase)
            % TESTRECORDACCUMULATION Validate timing accumulation for same record.

            % Time the same named record multiple times
            testCase.Timer.start(record="iterative");
            pause(testCase.PauseDuration);
            testCase.Timer.stop(record="iterative");

            firstDuration = testCase.Timer.Records(1).Duration;

            testCase.Timer.start(record="iterative");
            pause(testCase.PauseDuration);
            testCase.Timer.stop(record="iterative");

            % Should still have only one record, but with accumulated duration
            testCase.verifyNumElements(testCase.Timer.Records, 1);
            testCase.verifyEqual(testCase.Timer.Records(1).Name, "iterative");
            testCase.verifyGreaterThan(testCase.Timer.Records(1).Duration, firstDuration);
        end

        function testRecordNotFound(testCase)
            % TESTRECORDNOTFOUND Validate error handling for non-existent records.

            % Try to stop a record that doesn't exist
            testCase.verifyError(@() testCase.Timer.stop(record="nonexistent"), ...
                'core:chrono:Timer:RecordNotFound');
        end

        function testNamedRecordWithMessages(testCase)
            % TESTNAMEDRECORDWITHMESSAGES Validate named records with verbose output.

            verboseTimer = core.chrono.Timer(verbose=true);

            output = evalc(['verboseTimer.start(msg="Starting computation", record="compute"); ', ...
                'pause(testCase.PauseDuration); ', ...
                'verboseTimer.stop(msg="Computation complete", record="compute");']);

            testCase.verifySubstring(output, 'Starting computation');
            testCase.verifySubstring(output, 'Computation complete');
            testCase.verifyEqual(verboseTimer.Records(1).Name, "compute");
            testCase.verifyGreaterThan(verboseTimer.Records(1).Duration, 0);
        end

        function testMixedRecordTypes(testCase)
            % TESTMIXEDRECORDTYPES Validate mixing anonymous and named records.

            % Start with anonymous record
            testCase.Timer.start();
            pause(testCase.PauseDuration);
            testCase.Timer.stop();

            % Add named record
            testCase.Timer.start(record="named_task");
            pause(testCase.PauseDuration);
            testCase.Timer.stop(record="named_task");

            % Add another anonymous record
            testCase.Timer.start();
            pause(testCase.PauseDuration);
            testCase.Timer.stop();

            testCase.verifyNumElements(testCase.Timer.Records, 3);
            testCase.verifyEqual(testCase.Timer.Records(1).Name, "Record_1");
            testCase.verifyEqual(testCase.Timer.Records(2).Name, "named_task");
            testCase.verifyEqual(testCase.Timer.Records(3).Name, "Record_3");
        end
    end
end