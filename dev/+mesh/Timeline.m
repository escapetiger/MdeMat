classdef Timeline < handle
    % TIMELINE Abstract base class for temporal discretization.
    %
    %   Timeline provides a foundation for representing sequences of time
    %   points and managing temporal progression in numerical simulations.
    %   This abstract class defines the common interface and basic
    %   functionality that must be implemented by concrete timeline
    %   classes.
    %
    % See also:
    %   approx.mesh.StaticTimeline, approx.mesh.DynamicTimeline
    
    properties (Access = public)
        final % Final time (positive scalar)
        now   % Current time (scalar)
        count % Event counter (positive integer)
    end

    methods
        function obj = Timeline(final)
            % TIMELINE Constructor for Timeline.
            %
            %   obj = Timeline(final) creates a timeline object with
            %   the specified final time. This constructor is called
            %   by concrete subclasses during their initialization.
            %
            % Inputs:
            %   final - Final time (positive scalar)
            %
            % Outputs:
            %   obj - Constructed Timeline object
            
            core.except.assert(isscalar(final) && final >= 0, ...
                'InvalidInput', 'Final time must be a nonnegative scalar.');
            
            obj.final = final;
            obj.now = 0;
            obj.count = 1;
        end

        function obj = reset(obj)
            % RESET Reset the timeline to initial state.
            %
            %   obj = reset(obj) resets the timeline to time zero and
            %   reinitializes the event counter.
            %
            % Inputs:
            %   obj - The Timeline object
            %
            % Outputs:
            %   obj - The Timeline object
            
            obj.now = 0;
            obj.count = 1;
        end
    end

    methods (Abstract)
        % ADVANCE Advance the timeline by one time step.
        obj = advance(obj)
    end
end