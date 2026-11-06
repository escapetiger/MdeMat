classdef HistoryIntegratorProbe < mdemat.common.ode.OdeIntegrator
    %HISTORYINTEGRATORPROBE Concrete two-stage/two-history observation fixture.
    methods
        function obj = HistoryIntegratorProbe()
            obj@mdemat.common.ode.OdeIntegrator(2, 2, 1);
        end
        function value = step(obj, varargin) %#ok<INUSD>
            value = obj.U0{1, 1};
        end
    end
end
