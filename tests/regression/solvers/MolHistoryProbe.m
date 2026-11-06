classdef MolHistoryProbe < mdemat.internal.execution.MolScheme
    %MOLHISTORYPROBE Expose the protected temporal-view boundary for regression.
    methods
        function [values, count] = history(obj, space, dofs, varargin)
            [values, count] = obj.precompute(space, dofs, varargin{:});
        end
    end
end
