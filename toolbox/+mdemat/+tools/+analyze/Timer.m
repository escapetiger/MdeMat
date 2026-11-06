classdef Timer < handle
    %TIMER Accumulate mutually exclusive named wall-clock intervals.
    %   start/stop accumulate elapsed seconds for a named record. Anonymous
    %   starts create unique names. Records and find return value snapshots:
    %   changing a snapshot cannot change this Timer. Duration includes only
    %   stopped intervals, never a currently running interval.
    %
    %   [~, interval] = timer.start(record="solve") returns an interval ID.
    %   onCleanup(@() timer.cancel(interval)) safely cancels that interval on
    %   failure without stopping a newer interval or counting failed work.
    %   A successful timer.stop(interval=interval) commits elapsed time.

    properties
        Verbose = true
    end
    properties (Dependent, SetAccess = private)
        Records
    end
    properties (Access = private)
        Items = struct('Name', {}, 'Duration', {})
        ActiveIndex = 0
        ActiveToken = []
        Generation = uint64(0)
    end

    methods
        function obj = Timer(options)
            arguments
                options.names = {}
                options.verbose = true
            end
            names = recordNames(options.names);
            if numel(unique(names)) ~= numel(names)
                error('mdemat:tools:timing:Timer:DuplicateName', ...
                    'Timing record names must be unique.');
            end
            obj.Verbose = options.verbose;
            for index = 1:numel(names)
                obj.Items(end+1) = struct('Name', names(index), 'Duration', 0);
            end
        end

        function set.Verbose(obj, value)
            mustBeBooleanScalar(value);
            obj.Verbose = logical(value);
        end

        function records = get.Records(obj)
            records = obj.Items;
            for index = 1:numel(records)
                records(index).IsActive = index == obj.ActiveIndex;
            end
            if isempty(records)
                records = struct('Name', {}, 'Duration', {}, 'IsActive', {});
            end
        end

        function record = find(obj, name)
            %FIND Return a record snapshot, or [] when its name is absent.
            arguments
                obj (1, 1) mdemat.tools.analyze.Timer
                name {mustBeTextScalar}
            end
            name = recordName(name);
            index = obj.findIndex(name);
            record = [];
            if index > 0
                records = obj.Records;
                record = records(index);
            end
        end

        function [obj, interval] = start(obj, options)
            %START Begin one interval; message output precedes measurement.
            arguments
                obj (1, 1) mdemat.tools.analyze.Timer
                options.record {mustBeTextScalar} = ""
                options.msg {mustBeTextScalar} = ""
            end
            validateMessages(options.record, options.msg);
            if obj.ActiveIndex > 0
                error('mdemat:tools:timing:Timer:recordAlreadyActive', ...
                    'Record "%s" is already running.', obj.Items(obj.ActiveIndex).Name);
            end
            if strlength(options.record) == 0
                number = numel(obj.Items) + 1;
                name = "Record_" + number;
                while obj.findIndex(name) > 0
                    number = number + 1;
                    name = "Record_" + number;
                end
            else
                name = recordName(options.record);
            end
            index = obj.findIndex(name);
            if index == 0
                obj.Items(end+1) = struct('Name', name, 'Duration', 0);
                index = numel(obj.Items);
            end
            if strlength(options.msg) > 0 && obj.Verbose, disp(options.msg); end
            obj.Generation = obj.Generation + uint64(1);
            interval = obj.Generation;
            obj.ActiveIndex = index;
            obj.ActiveToken = tic;
        end

        function obj = stop(obj, options)
            %STOP Commit the selected active interval before displaying text.
            arguments
                obj (1, 1) mdemat.tools.analyze.Timer
                options.record {mustBeTextScalar} = ""
                options.msg {mustBeTextScalar} = ""
                options.interval = []
            end
            validateMessages(options.record, options.msg);
            if strlength(options.record) > 0
                index = obj.findIndex(recordName(options.record));
                if index == 0
                    error('mdemat:tools:timing:Timer:RecordNotFound', ...
                        'Timing record "%s" was not found.', options.record);
                end
                if index ~= obj.ActiveIndex
                    error('mdemat:tools:timing:Timer:RecordNotActive', ...
                        'Timing record "%s" is not active.', options.record);
                end
            else
                index = obj.ActiveIndex;
                if index == 0
                    error('mdemat:tools:timing:Timer:NoActiveRecord', ...
                        'No timing record is active.');
                end
            end
            if ~isempty(options.interval)
                validateInterval(options.interval);
                if options.interval ~= obj.Generation
                    error('mdemat:tools:timing:Timer:WrongInterval', ...
                        'The requested timing interval is no longer active.');
                end
            end
            elapsed = toc(obj.ActiveToken);
            obj.Items(index).Duration = obj.Items(index).Duration + elapsed;
            obj.ActiveIndex = 0;
            obj.ActiveToken = [];
            if strlength(options.msg) > 0 && obj.Verbose, disp(options.msg); end
        end

        function cancel(obj, interval)
            %CANCEL Discard only this still-active interval; stale IDs are safe.
            arguments
                obj (1, 1) mdemat.tools.analyze.Timer
                interval
            end
            validateInterval(interval);
            if obj.ActiveIndex > 0 && obj.Generation == interval
                obj.ActiveIndex = 0;
                obj.ActiveToken = [];
            end
        end

        function obj = reset(obj)
            %RESET Clear durations, retain names and cancel any active work.
            if obj.ActiveIndex > 0
                warning('mdemat:tools:timing:Timer:ResetWhileActive', ...
                    'Reset cancels the active timing interval.');
            end
            obj.ActiveIndex = 0;
            obj.ActiveToken = [];
            for index = 1:numel(obj.Items), obj.Items(index).Duration = 0; end
        end

        function result = snapshot(obj)
            %SNAPSHOT Return stopped seconds and activity without side effects.
            result = struct('Records', obj.Records, ...
                'Total', sum([obj.Items.Duration]), 'IsActive', obj.ActiveIndex > 0, ...
                'ActiveName', "", 'Unit', "s", 'Clock', "wall");
            if result.IsActive, result.ActiveName = obj.Items(obj.ActiveIndex).Name; end
        end

        function obj = report(obj, options)
            %REPORT Print a snapshot to stdout or an already open writable file.
            arguments
                obj (1, 1) mdemat.tools.analyze.Timer
                options.perc = true
                options.file = 1
            end
            mustBeBooleanScalar(options.perc);
            validateFile(options.file);
            result = obj.snapshot();
            fprintf(options.file, '[R] Wall-clock report (stopped intervals).\n');
            for index = 1:numel(result.Records)
                record = result.Records(index);
                if options.perc
                    percent = 0;
                    if result.Total > 0, percent = 100 * record.Duration / result.Total; end
                    fprintf(options.file, '[R] %-18s : %12.6gs %5.1f%%\n', ...
                        record.Name, record.Duration, percent);
                else
                    fprintf(options.file, '[R] %-18s : %12.6gs\n', record.Name, record.Duration);
                end
            end
            fprintf(options.file, '[R] %-18s : %12.6gs\n', 'Total', result.Total);
        end
    end

    methods (Access = private)
        function index = findIndex(obj, name)
            if isempty(obj.Items), index = 0; return; end
            index = find([obj.Items.Name] == name, 1);
            if isempty(index), index = 0; end
        end
    end
end

function name = recordName(name)
if ~((ischar(name) && isrow(name)) || (isstring(name) && isscalar(name))) || ...
        ismissing(string(name)) || strlength(string(name)) == 0
    error('mdemat:tools:timing:Timer:InvalidName', ...
        'Record names must be nonempty text scalars.');
end
name = string(name);
end

function names = recordNames(names)
if isempty(names), names = strings(1, 0); return; end
if ischar(names)
    names = {names};
elseif isstring(names)
    names = num2cell(names);
end
if ~iscell(names)
    error('mdemat:tools:timing:Timer:InvalidName', 'Record names must be a text list.');
end
names = cellfun(@recordName, names(:), 'UniformOutput', false);
names = reshape([names{:}], 1, []);
end

function mustBeBooleanScalar(value)
if ~(islogical(value) || isnumeric(value)) || ~isscalar(value) || ...
        ~isreal(value) || ~isfinite(value) || ~(value == 0 || value == 1)
    error('mdemat:tools:timing:Timer:InvalidBoolean', ...
        'Boolean options must be logical scalars or numeric 0/1.');
end
end

function validateInterval(interval)
if ~isa(interval, 'uint64') || ~isscalar(interval) || interval == 0
    error('mdemat:tools:timing:Timer:InvalidInterval', ...
        'An interval ID must be the uint64 scalar returned by start.');
end
end

function validateMessages(record, message)
if ismissing(string(record)) || ismissing(string(message))
    error('mdemat:tools:timing:Timer:InvalidName', ...
        'Record names and status messages cannot be missing.');
end
end

function validateFile(file)
if ~isnumeric(file) || ~isscalar(file) || ~isreal(file) || ...
        ~isfinite(file) || fix(file) ~= file || file < 1
    error('mdemat:tools:timing:Timer:InvalidFile', 'Output requires a writable file ID.');
end
if file == 1 || file == 2, return; end
try
    [name, permission] = fopen(file);
catch
    name = '';
    permission = '';
end
if isempty(name) || ~(contains(permission, 'w') || contains(permission, 'a') || contains(permission, '+'))
    error('mdemat:tools:timing:Timer:InvalidFile', 'Output requires a writable file ID.');
end
end
