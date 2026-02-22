classdef Timer < handle
    % TIMER High-precision timing framework for performance measurement.
    %
    %   Timer manages a collection of named timing records with mutually
    %   exclusive operation (only one record can be active at any time).
    %   It provides comprehensive timing functionality including start/stop
    %   operations, duration accumulation, and formatted reporting.
    %
    %   The Timer class maintains an array of Record objects, each tracking
    %   a specific operation or phase. When a record is started, all other
    %   records must be stopped to ensure accurate measurements. Timing
    %   data accumulates across multiple start/stop cycles for the same
    %   record, allowing measurement of iterative or multi-phase
    %   operations.
    %
    % See also:
    %   core.chrono.Record

    properties
        Records % Array of Record objects for timing measurement
        Verbose % Controls display of status messages
    end

    methods
        function obj = Timer(options)
            % TIMER Construct Timer instance with optional timing records.
            %
            %   obj = Timer() creates an empty Timer object with no records
            %   and enables verbose output.
            %
            %   obj = Timer(names={"rec1", "rec2"}) creates a Timer object
            %   with the specified record names and enables verbose output.
            %
            %   obj = Timer(verbose=false) creates a Timer object with no
            %   records and disables message output.
            %
            %   obj = Timer(names={"rec1"}, verbose=false) creates a Timer
            %   with specified records and controls message output.

            arguments
                options.names {mustBeA(options.names, 'cell')} = {}
                options.verbose {mustBeNumericOrLogical} = true
            end

            if ~isempty(options.names)
                n = length(options.names);
                obj.Records = arrayfun(@(k) core.chrono.Record(options.names{k}), 1:n);
            else
                obj.Records = core.chrono.Record.empty(0, 1);
            end

            obj.Verbose = options.verbose;
        end

        function record = find(obj, name)
            % FIND Find the record with the specified name.
            %
            %   record = find(obj, name) returns record by name.

            record = [];
            for i = 1:length(obj.Records)
                if strcmp(obj.Records(i).Name, name)
                    record = obj.Records(i);
                    break;
                end
            end
        end

        function obj = start(obj, options)
            % START Begin timing measurement by creating or reusing a
            % record.
            %
            %   obj = start(obj) creates a new anonymous record and starts
            %   timing.
            %
            %   obj = start(obj, msg="message") creates a new anonymous
            %   record, starts timing, and displays the message if verbose
            %   output is enabled.
            %
            %   obj = start(obj, record="recordName") starts timing for the
            %   specified named record. Creates the record if it doesn't
            %   exist.
            %
            %   obj = start(obj, msg="message", record="recordName") starts
            %   timing for the specified record with message display.

            arguments
                obj core.chrono.Timer
                options.msg {mustBeTextScalar} = ""
                options.record {mustBeTextScalar} = ""
            end

            [isRunning, activeRecord] = obj.isAnyRecordActive();
            if isRunning
                core.except.assert(0, 'recordAlreadyActive', ...
                    ['Record %s is already running.', ...
                    ' Please stop it before starting a new record.'], ...
                    activeRecord.Name);
            end

            if strlength(options.record) > 0
                idx = obj.findOrCreateRecord(options.record);
            else
                idx = obj.createAnonymousRecord();
            end

            if strlength(options.msg) > 0 && obj.Verbose
                disp(options.msg);
            end

            obj.Records(idx).StartTime = tic;
        end


        function obj = stop(obj, options)
            % STOP Complete timing measurement and accumulate duration.
            %
            %   obj = stop(obj) stops timing for the currently active record.
            %
            %   obj = stop(obj, msg="message") stops timing for the currently
            %   active record and displays the message if verbose output is
            %   enabled.
            %
            %   obj = stop(obj, record="recordName") stops timing for the
            %   specified named record.
            %
            %   obj = stop(obj, msg="message", record="recordName") stops
            %   timing for the specified record with message display.

            arguments
                obj core.chrono.Timer
                options.msg {mustBeTextScalar} = ""
                options.record {mustBeTextScalar} = ""
            end

            if strlength(options.record) > 0
                idx = obj.findRecordIndex(options.record);
                core.except.assert(idx > 0, 'RecordNotFound', ...
                    'Record "%s" not found.', options.record);
            else
                idx = obj.getActiveRecordIndex();
            end

            if strlength(options.msg) > 0 && obj.Verbose
                disp(options.msg);
            end

            core.except.assert(obj.Records(idx).IsActive, ...
                'RecordNotActive', ...
                'Record "%s" was not started.', obj.Records(idx).Name)

            elapsed = toc(obj.Records(idx).StartTime);
            obj.Records(idx).Duration = obj.Records(idx).Duration + elapsed;
            obj.Records(idx).StartTime = [];
        end


        function obj = reset(obj)
            % RESET Clear all accumulated timing data.
            %
            %   obj = reset(obj) resets the accumulated duration of all
            %   records to zero and clears any active timing operations.

            arguments
                obj core.chrono.Timer
            end

            [isRunning, activeRecord] = obj.isAnyRecordActive();

            if isRunning && obj.Verbose
                core.except.verify(0, 'ResetWhileActive', ...
                    'Resetting timer while record "%s" is active.', ...
                    activeRecord.Name);
            end

            for i = 1:length(obj.Records), obj.Records(i).reset(); end
        end

        function obj = report(obj, options)
            % REPORT Generate formatted timing performance summary.
            %
            %   obj = report(obj) generates a formatted report of all
            %   timing records with durations and percentages.
            %
            %   obj = report(obj, perc=false) controls percentage display
            %   with the perc option.
            %
            %   obj = report(obj, file=fid) specifies output destination
            %   with file identifier.

            arguments
                obj core.chrono.Timer
                options.perc {mustBeNumericOrLogical} = true
                options.file {mustBeInteger} = 1
            end

            totalTime = sum([obj.Records.Duration]);

            fprintf(options.file, '[R] CPU-times Report.\n');

            if options.perc
                formatStr = '[R] %-18s : %12.2fs %5.0f%%\n';
                for i = 1:length(obj.Records)
                    rec = obj.Records(i);
                    if totalTime > 0
                        pct = (rec.Duration / totalTime) * 100;
                    else
                        pct = 0;
                    end
                    fprintf(options.file, formatStr, rec.Name, rec.Duration, pct);
                end
                fprintf(options.file, '[R] %-18s : %12.2fs %5.0f%%\n', ...
                    'Total', totalTime, 100);
            else
                formatStr = ' %-18s : %12.2fs\n';
                for i = 1:length(obj.Records)
                    rec = obj.Records(i);
                    fprintf(options.file, formatStr, rec.Name, rec.Duration);
                end
                fprintf(options.file, '[R] %-18s : %12.2fs\n', ...
                    'Total', totalTime);
            end
        end
    end

    methods (Access = private)
        function [isRunning, activeRecord] = isAnyRecordActive(obj)
            % ISANYRECORDACTIVE Determine if any timing operation is in
            % progress.

            isRunning = false;
            activeRecord = [];

            for i = 1:length(obj.Records)
                if obj.Records(i).IsActive
                    isRunning = true;
                    activeRecord = obj.Records(i);
                    return;
                end
            end
        end

        function idx = createAnonymousRecord(obj)
            % CREATEANONYMOUSRECORD Create new anonymous record and return
            % its index.

            recordName = sprintf('Record_%d', length(obj.Records)+1);
            newRecord = core.chrono.Record(recordName);
            obj.Records(end+1) = newRecord;
            idx = length(obj.Records);
        end

        function idx = getActiveRecordIndex(obj)
            % GETACTIVERECORDINDEX Find index of currently active record.

            [isRunning, activeRecord] = obj.isAnyRecordActive();
            core.except.assert(isRunning, 'NoActiveRecord', ...
                'No record is currently active.');

            names = [obj.Records.Name];
            idx = find(strcmp(names, activeRecord.Name), 1);
        end

        function idx = findRecordIndex(obj, recordName)
            % FINDRECORDINDEX Find index of record by name.

            if isempty(obj.Records)
                idx = 0;
                return;
            end

            names = [obj.Records.Name];
            idx = find(strcmp(names, recordName), 1);
            if isempty(idx)
                idx = 0;
            end
        end

        function idx = findOrCreateRecord(obj, recordName)
            % FINDORCREATERECORD Find existing record or create new one.

            idx = obj.findRecordIndex(recordName);
            if idx == 0
                newRecord = core.chrono.Record(recordName);
                obj.Records(end+1) = newRecord;
                idx = length(obj.Records);
            end
        end
    end
end