classdef Constant < core.function.Function
    % CONSTANTFUNCTION Function that returns constant values.
    %
    %   Constant represents a mathematical mapping that returns the same
    %   constant value for all input points:
    %
    %   \f[
    %     f(x) = c
    %   \f]
    %
    %   where \f$c\f$ is a constant vector in \f$R^m\f$.
    %
    % See also:
    %   core.function.Function, core.function.Linear

    properties
        c % The constant value returned by the function
    end

    methods
        function obj = Constant(nDims, constant)
            % CONSTANT Constructor for the Constant class.
            %
            %   obj = Constant(nDims, constant) creates a Constant that
            %   maps from \\f$R^n\\f$ to \\f$R^m\\f$, where m is determined
            %   by the size of @a constant.

            arguments
                nDims{mustBePositive, mustBeInteger}
                constant{mustBeNumeric}
            end

            nCodims = length(constant);
            obj@core.function.Function(nDims = nDims, nCodims = nCodims);
            obj.c = constant(:);
        end
    end

    methods (Access = protected)
        function Y = evalImpl(obj, X)
            % EVALIMPL Implementation of constant function evaluation.

            Y = repmat(obj.c, 1, size(X, 2));
        end

        function dY = diffImpl(obj, X, ~)
            % DIFFIMPL Implementation of constant function derivative.

            dY = zeros(obj.NCodims, size(X, 2));
        end
    end
end