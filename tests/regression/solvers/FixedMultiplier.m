classdef FixedMultiplier
    %FIXEDMULTIPLIER Supply an independently chosen spatial force matrix.
    properties
        Matrix
    end
    methods
        function obj = FixedMultiplier(matrix)
            obj.Matrix = matrix;
        end
        function matrix = assembleMatrix(obj, ~)
            matrix = obj.Matrix;
        end
    end
end
