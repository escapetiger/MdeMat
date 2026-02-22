classdef Record < handle
    % RECORD Single timing record with state tracking.
    %
    %   Record represents a named timing measurement that maintains active
    %   or inactive state and accumulates elapsed time.
    %
    % Examples:
    %   % Basic usage
    %   rec = core.utilities.Record('Computation');
    %   rec.startTime = tic;
    %   pause(0.1);  % Simulate computation
    %   rec.duration = rec.duration + toc(rec.startTime);
    %   rec.startTime = [];  % Mark as inactive
    %
    %   % Check if active
    %   if ~rec.isActive()
    %     rec.startTime = tic;  % Start timing again
    %   end
    %
    % See Also:
    %   core.chrono.Timer

    properties
        name % Record name
        duration % Accumulated elapsed time in seconds
        startTime % Timing token when active, empty when inactive
    end

    methods
        function obj = Record(name)
            % RECORD Constructor for Record class.
            %
            %   obj = Record(name) creates a new Record object with the
            %   specified @a name and initializes duration to zero.
            %
            % Inputs:
            %   name - Name identifier
            %
            % Outputs:
            %   obj - The constructed Record object

            core.except.assert(nargin >= 1, 'InvalidInput', ...
                'Record name must be provided.');
            obj.name = name;
            obj.reset();
        end

        function tf = isActive(obj)
            % ISACTIVE Return true if timing is currently active.
            %
            %   tf = isActive(obj) determines whether the Record is
            %   currently in an active timing state.
            %
            % Inputs:
            %   obj - The Record object
            %
            % Outputs:
            %   tf - True if the record is active, false otherwise

            tf = ~isempty(obj.startTime);
        end

        function obj = reset(obj)
            % RESET Reset record duration to zero.
            %
            %   obj = reset(obj) resets the accumulated duration to zero,
            %   clearing the timing history for this Record.
            %
            % Inputs:
            %   obj - The Record object
            %
            % Outputs:
            %   obj - The Record object

            obj.duration = 0;
            obj.startTime = [];
        end
    end
end