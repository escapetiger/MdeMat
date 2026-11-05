classdef ConstantFunction < core.function.Function
    % CONSTANTFUNCTION Function that returns constant values.
    %
    %   ConstantFunction represents a mathematical mapping that returns
    %   the same constant value for all input points:
    %   
    %   \f[
    %     f(x) = c
    %   \f]
    %
    %   where \f$c\f$ is a constant vector in \f$R^m\f$.
    %
    % Examples:
    %   % Scalar constant function
    %   f = ConstantFunction(3, 5);  % R^3 -> R^1, returns 5
    %
    %   % Vector constant function
    %   f = ConstantFunction(2, [1; 2; 3]);  % R^2 -> R^3, returns [1;2;3]
    %
    %   % Evaluate
    %   X = rand(3, 100);
    %   Y = f.evaluate(X);  % Y will be 5 for all points
    %
    % See Also:
    %   core.function.Function, core.function.LinearFunction

    properties
        constant % The constant value returned by the function
    end

    methods
        function obj = ConstantFunction(nDims, constant)
            % CONSTANTFUNCTION Constructor for the ConstantFunction class.
            %
            %   obj = ConstantFunction(nDims, constant) creates a
            %   ConstantFunction that maps from R^nDims to R^m, where m
            %   is determined by the size of constant.
            %
            % Inputs:
            %   nDims - Number of input dimensions
            %   constant - Constant value (scalar or column vector)
            %
            % Outputs:
            %   obj - The created ConstantFunction object

            constant = constant(:);
            nCodims = length(constant);

            obj@core.function.Function(nDims, nCodims);
            obj.constant = constant;
        end
    end

    methods (Access = protected)
        function Y = evaluateImpl(obj, X)
            % EVALUATEIMPL Implementation of constant function evaluation.
            %
            %   Y = evaluateImpl(obj, X) returns the constant value
            %   replicated for each input point.
            %
            % Inputs:
            %   obj - The ConstantFunction object
            %   X - Coordinates (nDims x nPoints)
            %
            % Outputs:
            %   Y - Constant values (nCodims x nPoints)

            nPoints = size(X, 2);
            Y = repmat(obj.constant, 1, nPoints);
        end

        function dY = derivativeImpl(obj, X, r)
            % DERIVATIVEIMPL Implementation of constant function
            % derivative.
            %
            %   dY = derivativeImpl(obj, X, r) returns zeros since the
            %   derivative of a constant function is always zero.
            %
            % Inputs:
            %   obj - The ConstantFunction object
            %   X - Coordinates (nDims x nPoints)
            %   r - Derivative order vector (1 x nDims)
            %
            % Outputs:
            %   dY - Zero derivatives (nCodims x nPoints)

            nPoints = size(X, 2);
            dY = zeros(obj.nCodims, nPoints);
        end
    end
end