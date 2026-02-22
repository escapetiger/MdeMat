classdef StaticTimeline < approx.mesh.Timeline
    % STATICTIMELINE Predefined time node timeline.
    %
    %   StaticTimeline represents a temporal discretization with predefined
    %   time points. It advances through these points sequentially and
    %   provides information about the current state and properties of
    %   the timeline, including uniformity checks and time step queries.
    %
    %   Time nodes must be provided in strictly ascending order. The
    %   class supports both uniform and non-uniform time stepping and
    %   can detect uniformity automatically.
    %
    % See also:
    %   approx.mesh.Timeline, approx.mesh.DynamicTimeline
    
    properties (Access = public)
        Nodes % Prescribed time nodes (row vector)
    end
    
    properties (Dependent)
        StepSize % Current time step
        IsExhausted % Whether the final time has been reached
        NTotalTimeSteps % Total number of time steps
    end
    
    methods
        function obj = StaticTimeline(varargin)
            % STATICTIMELINE Constructor for StaticTimeline.
            %
            %   obj = StaticTimeline(nodes) creates a timeline with the
            %   specified time nodes.
            %
            %   obj = StaticTimeline(nSteps, final) creates a uniform
            %   timeline with nSteps from 0 to final time.
            
            core.except.assert(nargin >= 1 && nargin <= 2, ...
                'InvalidInput', ...
                'StaticTimeline requires 1 or 2 input arguments.');
            
            if nargin == 1
                nodes = varargin{1};
                
                core.except.assert(isvector(nodes) && isnumeric(nodes), ...
                    'InvalidInput', ...
                    'nodes must be a numeric vector.');
                
                core.except.assert(all(diff(nodes) > 0), ...
                    'InvalidInput', ...
                    'Time nodes must be in strictly ascending order.');
                
                finalTime = nodes(end);
                
            else
                nSteps = varargin{1};
                finalTime = varargin{2};
                
                core.except.assert(isscalar(nSteps) && nSteps > 0 && mod(nSteps, 1) == 0, ...
                    'InvalidInput', ...
                    'nSteps must be a positive integer.');
                
                core.except.assert(isscalar(finalTime) && finalTime >= 0, ...
                    'InvalidInput', ...
                    'final must be a nonnegative scalar.');
                
                nodes = linspace(0, finalTime, nSteps + 1);
            end
            
            obj@approx.mesh.Timeline(finalTime);
            obj.Nodes = [nodes(:).', inf];  % Add sentinel value
            obj.Now = obj.Nodes(1);
        end
        
        function obj = advance(obj)
            % ADVANCE Advance the timeline to the next time node.
            %
            %   advance(obj) advances the timeline to the next predefined
            %   time node and increments the event counter.
            
            core.except.assert(~obj.IsExhausted, 'TimeExceeded', ...
                'Cannot advance beyond the final time.');
            
            obj.Count = obj.Count + 1;
            obj.Now = obj.Nodes(obj.Count);
        end
        
        function TF = get.IsExhausted(obj)
            % GET.ISEXHAUSTED Check if timeline has reached final time.
            
            TF = obj.Count >= length(obj.Nodes);
        end
        
        function dt = get.StepSize(obj)
            % GET.DT Get the current time step size.
            
            if obj.IsExhausted
                dt = NaN;
            else
                dt = obj.Nodes(obj.Count + 1) - obj.Nodes(obj.Count);
            end
        end
        
        function n = get.NTotalTimeSteps(obj)
            % GET.NTOTALTIMESTEPS Get the total number of time steps.
            
            n = length(obj.Nodes) - 1;
        end
    end
end