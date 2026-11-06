classdef ObservationSchemeProbe < mdemat.solver.hyperbolic.UwdgScheme
    %OBSERVATIONSCHEMEPROBE Accepted samples u(t)=t at deliberately unequal steps.
    properties
        FailStep = false
    end
    properties (SetAccess = private)
        AcceptedTimes = []
    end
    methods
        function obj = ObservationSchemeProbe()
            obj@mdemat.solver.hyperbolic.UwdgScheme(config=struct( ...
                'advection',0,'ic',@(x) zeros(1,size(x,2)), ...
                'dt',0.1,'cfl',[],'tFinal',0.3));
        end
        function state = initialize(obj,state)
            obj.AcceptedTimes = [];
            state.setDof('U',state.XDisc.fit(0));
        end
        function state = preStep(obj,state)
            steps = [0.1,0.2];
            obj.TDisc.setTimeStep(dt=steps(min(obj.TDisc.Timeline.Count,2)));
        end
        function state = step(obj,state)
            if obj.FailStep
                error('mdemat:test:ObservationStepFailure','Deliberate integration failure.');
            end
            time = obj.TDisc.Timeline.Next;
            state.setDof('U',state.XDisc.fit(time));
            obj.AcceptedTimes(end+1) = time;
        end
        function state = finalize(obj,state) %#ok<INUSD>
        end
    end
end
