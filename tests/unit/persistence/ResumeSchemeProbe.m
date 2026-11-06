classdef ResumeSchemeProbe < mdemat.solver.hyperbolic.UwdgScheme
    % Deliberately destructive initialization exposes an accidental restart reset.
    methods
        function obj = ResumeSchemeProbe(config)
            obj@mdemat.solver.hyperbolic.UwdgScheme(config=config);
        end
        function state = initialize(obj, state)
            state.setDof('U', 999*ones(state.XDisc.NGlobalDofs, 1));
            state.setHistory('time', -1);
            obj.L = -speye(state.XDisc.NGlobalDofs);
        end
        function state = finalize(obj, state) %#ok<INUSD>
        end
    end
end
