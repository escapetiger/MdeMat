classdef Timer < handle
    % TIMER Measure elapsed time for multiple operations.
    %
    %   Timer provides a framework for measuring and reporting execution
    %   time across multiple operations. Manages a collection of Record
    %   objects with mutually exclusive timing (only one record active
    %   at a time).
    %
    % Examples:
    %   % Basic usage cycle
    %   timer = core.chrono.Timer({'Setup', 'Computation', 'IO'});
    %
    %   timer.startByName('Setup');
    %   pause(0.1);  % Simulate setup operations
    %   timer.stopByName('Setup');
    %
    %   timer.startByName('Computation', 'Starting main calculation...');
    %   pause(0.2);  % Simulate computation
    %   timer.stopByName('Computation');
    %
    %   % Generate report
    %   timer.report();  % Display timing summary with percentages
    %
    % See Also:
    %   core.chrono.Record

    properties
        records % Array of Record objects for timing measurement
        verbose % Controls display of status messages
    end

    methods
        function obj = Timer(names, verbose)
            % TIMER Constructor for Timer class.
            %
            %   obj = Timer(names) creates a Timer object with the
            %   specified record names and sets @a verbose = 1.
            %
            %   obj = Timer(names, verbose) creates a Timer object with the
            %   specified record names and verbose flag.
            %
            % Inputs:
            %   names - Record names
            %   verbose - Verbose flag (optional, default: 1) 
            %
            % Outputs:
            %   obj - The constructed Timer object

            if nargin < 2, verbose = true; end

            core.except.assert(nargin >= 1, 'InvalidInput', ...
                'A cell array of record names is required.');

            core.except.assert(iscell(names), 'InvalidInput', ...
                'Record names must be provided as a cell array.');

            n = length(names);
            obj.records = arrayfun(@(k) core.chrono.Record(names{k}), 1:n);
            obj.verbose = verbose;
        end

        function obj = start(obj, idx, msg)
            % START Start timing a record by index.
            %
            %   obj = start(obj, idx) starts timing the record at the specified
            %   @a idx.
            %
            %   obj = start(obj, idx, msg) starts timing the record at the
            %   specified @a idx and displays @a msg if @a verbose > 0.
            %
            % Inputs:
            %   obj - The Timer object
            %   idx - Record index
            %   msg - Message (optional, default: '') 
            %
            % Outputs:
            %   obj - The Timer object

            core.except.assert(nargin >= 2, 'MissingIndex', ...
                'Record index must be provided.');

            if nargin < 3, msg = ''; end

            n = numel(obj.records);
            core.except.assert(idx >= 1 && idx <= n, ...
                'InvalidIndex', 'Invalid record index: %d', idx);

            [isRunning, activeRecord] = obj.isAnyRecordActive();
            if isRunning
                core.except.assert(0, 'recordAlreadyActive', ...
                    ['Record %s is already running.', ...
                    ' Please stop it before starting a new record.'], ...
                    activeRecord.name);
            end

            if ~isempty(msg) && obj.verbose, disp(msg); end

            obj.records(idx).startTime = tic;
        end

        function obj = startByName(obj, name, msg)
            % STARTBYNAME Start timing a record by name.
            %
            %   obj = startByName(obj, name) starts a record identified by
            %   @a name.
            %
            %   obj = startByName(obj, name, msg) starts a record
            %   identified by @a name and displays @a msg if @a verbose >
            %   0.
            %
            % Inputs:
            %   obj - The Timer object
            %   name - Record name
            %   msg - Message (optional, default: '')
            %
            % Outputs:
            %   obj - The Timer object

            core.except.assert(nargin >= 2, 'MissingName', ...
                'Record name must be provided.');

            idx = obj.getRecordIndex(name);

            if nargin < 3
                obj.start(idx);
            else
                obj.start(idx, msg);
            end
        end

        function obj = stop(obj, idx, msg)
            % STOP Stop timing a record by index and accumulate elapsed
            % time.
            %
            %   obj = stop(obj, idx) stops the record at @a idx and adds
            %   elapsed time to the record's duration.
            %
            %   obj = stop(obj, idx, msg) stops the record at @a idx, adds
            %   elapsed time to the record's duration and displays @a msg
            %   if @a verbose > 0.
            %
            % Inputs:
            %   obj - The Timer object
            %   idx - Record index
            %   msg - Message (optional, default: '')
            %
            % Outputs:
            %   obj - The Timer object

            core.except.assert(nargin >= 2, 'MissingIndex', ...
                'Record index must be provided.');

            if nargin < 3, msg = ''; end

            n = numel(obj.records);
            core.except.assert(idx >= 1 && idx <= n, ...
                'InvalidIndex', 'Invalid record index: %d', idx);

            if ~isempty(msg) && obj.verbose, disp(msg); end

            core.except.assert(obj.records(idx).isActive(), ...
                'RecordNotActive', ...
                'Record "%s" was not started.', obj.records(idx).name)

            elapsed = toc(obj.records(idx).startTime);
            obj.records(idx).duration = obj.records(idx).duration + elapsed;
            obj.records(idx).startTime = [];
        end

        function obj = stopByName(obj, name, msg)
            % STOPBYNAME Stop timing a record by name.
            %
            %   obj = stopByName(obj, name) stops a record identified by @a
            %   name.
            %
            %   obj = stopByName(obj, name, msg) stops a record identified
            %   by @a name and displays @a msg.
            %
            % Inputs:
            %   obj - The Timer object
            %   name - Record name
            %   msg - Message (optional, default: '')
            %
            % Outputs:
            %   obj - The Timer object

            core.except.assert(nargin >= 2, ...
                'MissingName', 'Record name must be provided.');

            idx = obj.getRecordIndex(name);

            if nargin < 3
                obj.stop(idx);
            else
                obj.stop(idx, msg);
            end
        end

        function obj = reset(obj)
            % RESET Reset all record durations to zero.
            %
            %   obj = reset(obj) resets the accumulated duration of all
            %   records to zero. Issues warning if any record is currently
            %   active.
            %
            % Inputs:
            %   obj - The Timer object
            %
            % Outputs:
            %   obj - The Timer object

            [isRunning, activeRecord] = obj.isAnyRecordActive();
            if isRunning && obj.verbose
                core.except.verify(0, 'ResetWhileActive', ...
                    'Resetting timer while record "%s" is active.', ...
                    activeRecord.name);
            end

            for i = 1:length(obj.records), obj.records(i).reset(); end
        end

        function obj = report(obj, perc, file)
            % REPORT Print summary of all records and their durations.
            %
            %   obj = report(obj) generates formatted report of all timing
            %   records.
            %
            %   obj = report(obj, perc, file) controls percentage display
            %   and output destination.
            %
            % Inputs:
            %   obj - The Timer object
            %   perc - True for show percentages (optional, default: true)
            %   file - File identifier (optional, default: 1) 
            %
            % Outputs:
            %   obj - The Timer object

            if nargin < 2, perc = true; end

            if nargin < 3, file = 1; end

            totalTime = sum([obj.records.duration]);

            fprintf(file, 'CPU-times Summary\n');
            fprintf(file, '===============\n');

            if perc
                formatStr = ' %-18s : %12.2fs %5.0f%%\n';
                for i = 1:length(obj.records)
                    rec = obj.records(i);
                    if totalTime > 0
                        pct = (rec.duration / totalTime) * 100;
                    else
                        pct = 0;
                    end
                    fprintf(file, formatStr, rec.name, rec.duration, pct);
                end
                fprintf(file, '---------------\n');
                fprintf(file, ' %-18s : %12.2fs %5.0f%%\n', 'Total', totalTime, 100);
            else
                formatStr = ' %-18s : %12.2fs\n';
                for i = 1:length(obj.records)
                    rec = obj.records(i);
                    fprintf(file, formatStr, rec.name, rec.duration);
                end
                fprintf(file, '---------------\n');
                fprintf(file, ' %-18s : %12.2fs\n', 'Total', totalTime);
            end
        end
    end

    methods (Access = private)
        function idx = getRecordIndex(obj, name)
            % GETRECORDINDEX Get the index of a record by its name.

            names = {obj.records.name};
            idx = find(strcmp(names, name), 1);

            core.except.assert(~isempty(idx), 'RecordNotFound', ...
                'Record "%s" not found.', name);
        end

        function [isRunning, activeRecord] = isAnyRecordActive(obj)
            % ISANYRECORDACTIVE Check if any record is currently active.

            isRunning = false;
            activeRecord = [];

            for i = 1:length(obj.records)
                if obj.records(i).isActive()
                    isRunning = true;
                    activeRecord = obj.records(i);
                    return;
                end
            end
        end
    end
end