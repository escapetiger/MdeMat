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
        Final % Final time (positive scalar)
        Now % Current time (scalar)
        Count % Event counter (positive integer)
    end

    methods
        function obj = Timeline(final)
            % TIMELINE Constructor for Timeline.
            %
            %   obj = Timeline(final) creates a timeline object with
            %   the specified @a final time. This constructor is called
            %   by concrete subclasses during their initialization.

            arguments
                final{mustBeNonnegative}
            end

            obj.Final = final;
            obj.Now = 0;
            obj.Count = 1;
        end

        function obj = reset(obj)
            % RESET Reset the timeline to initial state.
            %
            %   obj = reset(obj) resets the timeline to time zero and
            %   reinitializes the event counter.

            obj.Now = 0;
            obj.Count = 1;
        end

        function tf = isBefore(obj, t)
            % ISBEFORE Check the current time is before the specified time.
            %
            %   tf = isBefore(obj, t) returns true if the current time does
            %   not exceed the specified time @a t.

            arguments
                obj approx.mesh.Timeline
                t(1, 1) {mustBeReal}
            end

            tol = 1e-14;
            tf = obj.Now - tol <= t;
        end
    end

    methods (Abstract)
        % ADVANCE Advance the timeline by one time step.
        obj = advance(obj)
    end
end