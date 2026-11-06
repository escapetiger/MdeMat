classdef PhaseBasis < mdemat.data.function.Function
    %PHASEBASIS Uniform unitary change of coordinates for Galerkin identities.
    properties
        Original
        Phase
    end
    properties (Dependent)
        NFactorCodims
    end
    methods
        function obj = PhaseBasis(original, phase)
            obj@mdemat.data.function.Function(nDims=original.NDims, nCodims=original.NCodims);
            obj.Original = original;
            obj.Phase = phase;
        end
        function counts = get.NFactorCodims(obj)
            counts = obj.Original.NFactorCodims;
        end
    end
    methods (Access = protected)
        function values = evalImpl(obj, points)
            values = obj.Phase*obj.Original.eval(points);
        end
        function values = diffImpl(obj, points, order)
            values = obj.Phase*obj.Original.diff(points, order);
        end
    end
end
