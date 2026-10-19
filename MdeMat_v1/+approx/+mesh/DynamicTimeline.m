classdef DynamicTimeline < approx.mesh.Timeline
    % DYNAMICTIMELINE Adaptive time stepping timeline.
    %
    %   DynamicTimeline represents a temporal discretization with
    %   dynamically computed time step sizes. This class is particularly
    %   useful for simulations where the time step needs to be adjusted
    %   based on spatial discretization constraints (e.g., CFL conditions)
    %   or solution behavior.
    %
    % Examples:
    %   % Create dynamic timeline with final time 10
    %   timeline = DynamicTimeline(10);
    %   
    %   % Set time step based on CFL condition
    %   timeline.setTimeStep(0.1, 0.5);  % h=0.1, CFL=0.5
    %   
    %   % Advance through time
    %   while ~timeline.isFinished()
    %       timeline.advance();
    %       % Update time step if needed
    %       timeline.setTimeStep(newH, cflNumber);
    %   end
    %
    % Notes:
    %   The time step is computed as dt = C * h^p where h is the spatial
    %   mesh parameter, C is a convergence constant (e.g., CFL number),
    %   and p is an exponent (default 1).
    %
    % See also:
    %   approx.mesh.Timeline, approx.mesh.StaticTimeline
    
    properties (Access = public)
        h0     % Previous time step sizes (from new to old)
        h      % Current time step sizes (from new to old)
        next   % Next time point (scalar)
        nSteps % Number of time steps to track (positive integer)
    end
    
    methods
        function obj = DynamicTimeline(final, nSteps)
            % DYNAMICTIMELINE Constructor for DynamicTimeline.
            %
            %   obj = DynamicTimeline(final) creates a dynamic timeline
            %   with the specified final time and default step tracking.
            %
            %   obj = DynamicTimeline(final, nSteps) creates a timeline
            %   that tracks the specified number of previous time steps.
            %
            % Inputs:
            %   final - Final time (positive scalar)
            %   nSteps - Number of previous steps to track (optional, default: 1)
            %
            % Outputs:
            %   obj - Constructed DynamicTimeline object
            
            if nargin < 2, nSteps = 1; end
            
            obj@approx.mesh.Timeline(final);
            obj.nSteps = nSteps;
            obj.reset();
        end

        function reset(obj)
            % RESET Reset the dynamic timeline to initial state.
            %
            %   reset(obj) resets the timeline to time zero and clears
            %   all time step history.
            %
            % Inputs:
            %   obj - The DynamicTimeline object
            
            reset@approx.mesh.Timeline(obj);
            obj.next = obj.now;            
            obj.h = zeros(1, obj.nSteps+1);
            obj.h0 = zeros(1, obj.nSteps+1);
        end

        function setTimeStep(obj, h, C)
            % SETTIMESTEP Set the current time step size.
            %
            %   setTimeStep(obj, h, C) computes the time step as dt = C*h
            %   and ensures it doesn't exceed the remaining time.
            %
            % Inputs:
            %   obj - The DynamicTimeline object
            %   h - Spatial mesh parameter (positive scalar)
            %   C - Convergence constant, e.g., CFL number (positive scalar)
     
            core.except.assert(C > 0 && h > 0, ...
                'InvalidInput', 'C and h must be positive.');

            p = max(obj.nSteps / obj.count, 1);

            dt = C * h^p;
            
            obj.h(1) = min(dt, obj.final - obj.now);
            obj.next = obj.now + obj.h(1);
        end
        
        function advance(obj)
            % ADVANCE Advance the timeline by one time step.
            %
            %   advance(obj) advances the current time by the current
            %   time step and updates the time step history.
            %
            % Inputs:
            %   obj - The DynamicTimeline object
            
            core.except.assert(~obj.isFinished(), ...
                'TimeExceeded', 'Cannot advance beyond the final time.');
            
            obj.now = obj.now + obj.h(1);
            obj.count = obj.count + 1;
            obj.h0 = obj.h;
            if numel(obj.h) > 1
                obj.h = circshift(obj.h, 1);
            end
            obj.next = obj.now + obj.h(1);
        end
        
        function TF = isFinished(obj)
            % ISFINISHED Check if the final time has been reached.
            %
            %   TF = isFinished(obj) returns true if the current time
            %   has reached or exceeded the final time.
            %
            % Inputs:
            %   obj - The DynamicTimeline object
            %
            % Outputs:
            %   TF - True if timeline is complete (logical scalar)
            
            TF = obj.now >= obj.final;
        end
        
        function TF = hasStepSizeChanged(obj)
            % HASSTEPSIZECHANGED Check if time step size has changed.
            %
            %   TF = hasStepSizeChanged(obj) returns true if the current
            %   time step differs from the previous one, or if this is
            %   the first step.
            %
            % Inputs:
            %   obj - The DynamicTimeline object
            %
            % Outputs:
            %   TF - True if step size changed (logical scalar)
            
            TF = ~isequal(obj.h, obj.h0) || (obj.count == 1);
        end
        
        function dt = getTimeStep(obj)
            % GETTIMESTEP Get the current time step size.
            %
            %   dt = getTimeStep(obj) returns the current time step
            %   size that will be used for the next advance.
            %
            % Inputs:
            %   obj - The DynamicTimeline object
            %
            % Outputs:
            %   dt - Current time step size (positive scalar)
            
            dt = obj.h(1);
        end
    end
end