classdef HermiteForceProbe < mdemat.application.plasma.DghsScheme
    %HERMITFORCEPROBE Expose the protected electric substep for regression.
    methods
        function obj = HermiteForceProbe(options)
            arguments
                options.config (1,1) struct = struct()
            end
            obj@mdemat.application.plasma.DghsScheme(config=options.config);
        end
        function state = forceStep(obj, state, dt, method)
            switch method
                case 'be'
                    state = obj.stepNonlinearBe(state, dt);
                case 'midpoint'
                    state = obj.stepNonlinearImp(state, dt);
                case 'sdirk2'
                    state = obj.stepNonlinearSdirk2(state, dt);
            end
        end
    end
end
