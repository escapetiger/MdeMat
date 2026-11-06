classdef ApHermiteForceProbe < mdemat.application.plasma.ApDghsScheme
    %APHERMITEFORCEPROBE Expose the equilibrium electric substep.
    methods
        function obj = ApHermiteForceProbe(options)
            arguments
                options.config (1,1) struct = struct()
            end
            obj@mdemat.application.plasma.ApDghsScheme(config=options.config);
        end
        function [constant, meanWeights] = gaugeWeights(obj, state)
            [constant, meanWeights] = obj.poissonGauge(state);
        end

        function state = forceStep(obj, state, dt, method)
            switch method
                case 'be'
                    state = obj.stepNonlinearBe(state, dt);
                case 'sdirk2'
                    state = obj.stepNonlinearSdirk2(state, dt);
            end
        end
    end
end
