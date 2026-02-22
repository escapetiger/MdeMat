classdef Record < handle
    % RECORD Individual timing measurement with state management.
    %
    %   Record represents a single named timing operation that tracks
    %   active/inactive state and accumulates elapsed time across multiple
    %   measurement sessions. Each record maintains a name identifier,
    %   cumulative duration, and timing token for active measurements.
    %
    %   The Record class provides the fundamental timing primitive used by
    %   the Timer class. It maintains state information to distinguish
    %   between active timing (when a measurement is in progress) and
    %   inactive timing (when no measurement is running). Duration values
    %   accumulate across multiple start/stop cycles.
    %
    % See also:
    %   core.chrono.Timer

    properties
        Name % Record name
        Duration % Accumulated elapsed time in seconds
        StartTime % Timing token when active, empty when inactive
    end

    properties (Dependent)
        IsActive % Logical flag indicating if timing is active
    end

    methods
        function obj = Record(name)
            % RECORD Construct Record instance with name identifier.
            %
            %   obj = Record(name) creates a new Record object with the
            %   specified @a name and initializes timing state to zero
            %   duration and inactive status.

            arguments
                name string
            end

            obj.Name = name;
            obj.reset();
        end

        function obj = reset(obj)
            % RESET Clear accumulated timing data and state.
            %
            %   obj = reset(obj) resets the accumulated duration to zero
            %   and clears any active timing operation.

            arguments
                obj core.chrono.Record
            end

            obj.Duration = 0;
            obj.StartTime = [];
        end

        function tf = get.IsActive(obj)
            % GET.ISACTIVE Determine if timing measurement is in progress.

            tf = ~isempty(obj.StartTime);
        end
    end
end