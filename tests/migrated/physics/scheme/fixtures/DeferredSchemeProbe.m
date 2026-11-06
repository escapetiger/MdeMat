classdef DeferredSchemeProbe < mdemat.internal.execution.Scheme
    methods
        function obj = DeferredSchemeProbe(varargin)
            obj@mdemat.internal.execution.Scheme(varargin{:});
        end
        function state = initialize(obj, state)
            obj.validateConfig('initial');
        end
    end
    methods (Access = protected)
        function obj = setup(obj)
            setup@mdemat.internal.execution.Scheme(obj);
            obj.addConfig('initial', default=[], ...
                validator=@(value) isa(value, 'function_handle'));
        end
    end
end
