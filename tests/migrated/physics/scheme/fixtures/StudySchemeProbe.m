classdef StudySchemeProbe < mdemat.internal.execution.Scheme
    methods
        function obj = StudySchemeProbe(varargin)
            obj@mdemat.internal.execution.Scheme(varargin{:});
        end
    end
    methods (Access = protected)
        function obj = setup(obj)
            setup@mdemat.internal.execution.Scheme(obj);
            obj.addConfig('L', default=2, ...
                validator=@(value) isnumeric(value) && isscalar(value) && value == 2);
            mdemat.internal.execution.defineStudyInputs(obj, 'plasma');
        end
    end
end
