classdef DynamicTimeline < approx.mesh.Timeline
    % DYNAMICTIMELINE Adaptive time stepping timeline.
    %
    %   DynamicTimeline represents a temporal discretization with
    %   dynamically computed time step sizes. This class is particularly
    %   useful for simulations where the time step needs to be adjusted
    %   based on spatial discretization constraints (e.g., CFL conditions)
    %   or solution behavior.
    %
    %   The time step is computed as \f$\Delta t = C * h^p\f$ where \f$h\f$
    %   is the spatial mesh parameter, \f$C\f$ is a convergence constant
    %   (e.g., CFL number), and \f$p\f$ is an exponent (default 1).
    %
    % See also:
    %   approx.mesh.Timeline, approx.mesh.StaticTimeline
    
    properties (Access = public)
        OldStepSizeQueue % Previous time step sizes (from new to old)
        StepSizeQueue % Current time step sizes (from new to old)
        Next % Next time point (scalar)
        NSteps % Number of time steps to track (positive integer)
    end
    
    properties (Dependent)
        StepSize % Current time step
        IsExhausted % Whether the final time has been reached
        HasStepSizeChanged % Whether the step size changed
    end
    
    methods
        function obj = DynamicTimeline(final, options)
            % DYNAMICTIMELINE Constructor for DynamicTimeline.
            %
            %   obj = DynamicTimeline(final) creates a dynamic timeline
            %   with the specified final time and default step tracking.
            %
            %   obj = DynamicTimeline(final, nSteps=nSteps) creates a
            %   timeline that tracks the specified number of previous time
            %   steps.
            
            arguments
                final {mustBeNonnegative}
                options.nSteps {mustBePositive, mustBeInteger} = 1
            end
            
            obj@approx.mesh.Timeline(final);
            obj.NSteps = options.nSteps;
            obj.reset();
        end
        
        function obj = reset(obj)
            % RESET Reset the dynamic timeline to initial state.
            %
            %   obj = reset(obj) resets the timeline to time zero and
            %   clears all time step history.
            
            arguments
                obj approx.mesh.DynamicTimeline
            end
            
            reset@approx.mesh.Timeline(obj);
            obj.Next = obj.Now;
            obj.StepSizeQueue = zeros(1, obj.NSteps+1);
            obj.OldStepSizeQueue = zeros(1, obj.NSteps+1);
        end
        
        function obj = setTimeStep(obj, options)
            % SETTIMESTEP Set the current time step size.
            %
            %   obj = setTimeStep(obj, h=h, C=C) sets time step as dt = C*h
            %   for dynamic timelines.
            %
            %   obj = setTimeStep(obj, h=h, C=C, p=p) sets time step as dt
            %   = C*h^p for dynamic timelines.
            %
            %   obj = setTimeStep(obj, dt=dt) directly sets the time step
            %   for manual control (use with caution).
            
            arguments
                obj approx.mesh.DynamicTimeline
                options.dt {mustBePositive} = []
                options.h {mustBeNumeric} = []
                options.C {mustBePositive, mustBeNonempty} = 1
                options.p {mustBePositive, mustBeNonempty} = 1
            end
            
            if isempty(options.dt) && isempty(options.h)
                core.except.assert(0, 'InvalidInput', ...
                    'Either dt or h must be provided to set the time step.');
            end
            
            if ~isempty(options.dt)
                dt = options.dt;
            else
                h = options.h;
                C = options.C;
                p = max(obj.NSteps/obj.Count, options.p);
                dt = C*h^p;
            end
            
            obj.StepSizeQueue(1) = min(dt, obj.Final - obj.Now);
            obj.Next = obj.Now + obj.StepSizeQueue(1);
        end
        
        function obj = advance(obj)
            % ADVANCE Advance the timeline by one time step.
            %
            %   obj = advance(obj) advances the current time by the current
            %   time step and updates the time step history.
            
            core.except.assert(~obj.IsExhausted, ...
                'TimeExceeded', 'Cannot advance beyond the final time.');
            
            obj.Now = obj.Now + obj.StepSizeQueue(1);
            obj.Count = obj.Count + 1;
            obj.OldStepSizeQueue = obj.StepSizeQueue;
            if numel(obj.StepSizeQueue) > 1
                obj.StepSizeQueue = circshift(obj.StepSizeQueue, 1);
            end
            obj.Next = obj.Now + obj.StepSizeQueue(1);
        end
        
        function TF = get.IsExhausted(obj)
            % GET.ISEXHAUSTED Check if the final time has been reached.
            
            TF = obj.Now >= obj.Final;
        end
        
        function TF = get.HasStepSizeChanged(obj)
            % GET.HASSTEPSIZECHANGED Check if time step size has changed.
            
            TF = ~isequal(obj.StepSizeQueue, obj.OldStepSizeQueue) || (obj.Count == 1);
        end
        
        function dt = get.StepSize(obj)
            % GET.DT Get the current time step size.
            
            dt = obj.StepSizeQueue(1);
        end
    end
end