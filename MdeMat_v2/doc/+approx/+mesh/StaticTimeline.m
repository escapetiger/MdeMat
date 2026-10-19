classdef StaticTimeline < approx.mesh.Timeline
    % STATICTIMELINE Predefined time node timeline.
    %
    %   StaticTimeline represents a temporal discretization with predefined
    %   time points. It advances through these points sequentially and
    %   provides information about the current state and properties of
    %   the timeline, including uniformity checks and time step queries.
    %
    % Examples:
    %   % Create timeline with uniform time steps
    %   timeline = StaticTimeline(10, 5.0);  % 10 steps to time 5.0
    %   
    %   % Create timeline with custom time nodes
    %   customNodes = [0, 0.1, 0.5, 1.0, 2.0, 5.0];
    %   timeline = StaticTimeline(customNodes);
    %   
    %   % Advance through timeline
    %   while ~timeline.isFinished()
    %       currentTime = timeline.now;
    %       dt = timeline.getTimeStep();
    %       timeline.advance();
    %   end
    %
    % Notes:
    %   Time nodes must be provided in strictly ascending order. The
    %   class supports both uniform and non-uniform time stepping and
    %   can detect uniformity automatically.
    %
    % See also:
    %   approx.mesh.Timeline, approx.mesh.DynamicTimeline
    
    properties (Access = public)
        nodes % Prescribed time nodes (row vector)
    end

    properties (Dependent)
        dt % Current time step
        isFinished % Whether the final time has been reached
        nTotalTimeSteps % Total number of time steps
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
            %
            % Inputs:
            %   varargin - Inputs arguments
            %<   nodes - Vector of time nodes in ascending order
            %<   nSteps - Number of time steps (positive integer)
            %<   final - Final time (positive scalar)
            %
            % Outputs:
            %   obj - Constructed StaticTimeline object
            
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
            obj.nodes = [nodes(:).', inf];  % Add sentinel value
            obj.now = obj.nodes(1);
        end
        
        function obj = advance(obj)
            % ADVANCE Advance the timeline to the next time node.
            %
            %   advance(obj) advances the timeline to the next predefined
            %   time node and increments the event counter.
            %
            % Inputs:
            %   obj - The StaticTimeline object
            %
            % Outputs:
            %   obj - The StaticTimeline object
            
            core.except.assert(~obj.isFinished(), 'TimeExceeded', ...
                'Cannot advance beyond the final time.');
            
            obj.count = obj.count + 1;
            obj.now = obj.nodes(obj.count);
        end
        
        function TF = get.isFinished(obj)
            % GET.ISFINISHED Check if timeline has reached final time.
            
            TF = obj.count >= length(obj.nodes);
        end
        
        function dt = get.dt(obj)
            % GET.DT Get the current time step size.
            
            if obj.isFinished
                dt = NaN;
            else
                dt = obj.nodes(obj.count + 1) - obj.nodes(obj.count);
            end
        end

        function n = get.nTotalTimeSteps(obj)
            % GET.NTOTALTIMESTEPS Get the total number of time steps.
            
            n = length(obj.nodes) - 1;
        end
    end
end