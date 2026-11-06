classdef SteadyScheme < mdemat.internal.execution.Scheme
    %STEADYSCHEME Solve once, with optional explicitly sampled observations.
    methods
        function state = run(obj,state)
            arguments
                obj mdemat.internal.execution.SteadyScheme
                state mdemat.common.fem.SpatialState
            end
            obj.prepareObservers(false);
            timer=obj.getObserver('timer'); timer.reset();
            % Elliptic initialization owns its Assembly interval. Timer has
            % one active interval, so the controller must not nest another.
            state=obj.initialize(state);
            state=timed(timer,'Integration',@()obj.step(state));
            timed(timer,'Visualization',@()obj.plot(state));
            state=obj.finalize(state);
        end
        function [state,tables] = converge(obj,state)
            [state,tables]=mdemat.internal.execution.runConvergence(obj,state);
        end
        function fields = samplePlotFields(obj,state,plan)
            fields=mdemat.common.fem.plotFieldSamples( ...
                state.XDisc,state.Dofs,obj.Observation.Plot.Components,plan);
        end
        function fields = samplePlotReferences(obj,state,plan,time)
            fields=mdemat.internal.execution.plotReferenceFields(state.XDisc, ...
                obj.Observation.Plot.Exacts,obj.Observation.Plot.Components,plan,time);
        end
    end
    methods (Access = protected)
        function obj = setup(obj)
            setup@mdemat.internal.execution.Scheme(obj);
            mdemat.internal.execution.defineStudyInputs(obj,'spatial');
        end
        function obj = plot(obj,state)
            options=obj.Observation.Plot;
            if isempty(options.Modes)||isempty(options.PointsPerDimension)|| ...
                    any(options.PointsPerDimension<=0),return;end
            datasets=mdemat.internal.execution.sampleObservation(obj,state,[]);
            mdemat.internal.execution.configurePlotLayouts(obj,datasets);
            session=obj.getObserver('visualizer');
            if ~isempty(session)&&session.IsEnabled,session.render(datasets,title=options.Title);end
        end
    end
end
function output = timed(timer,name,operation)
[~,interval]=timer.start(record=name);
guard=onCleanup(@()timer.cancel(interval));
output=operation();timer.stop(interval=interval);
end
