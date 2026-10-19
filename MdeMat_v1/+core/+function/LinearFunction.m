classdef LinearFunction < core.function.Function
    % LINEARFUNCTION Function representing affine transformations.
    %
    %   LinearFunction represents affine transformations of the form:
    %
    %   \f[
    %     f(x) = A*x + b
    %   \f]
    %
    %   where \f$A\f$ is an \f$m\times n\f$ matrix and \f$b\f$ is an
    %   \f$m\times 1\f$ vector. When \f$b\f$ is not provided or is zero,
    %   this reduces to a pure linear transformation.
    %
    % Examples:
    %   % Pure linear transformation f(x) = A*x
    %   A = [1 2; 3 4];
    %   f = LinearFunction(A);  % R^2 -> R^2
    %
    %   % Affine transformation f(x) = A*x + b
    %   b = [1; -1];
    %   f = LinearFunction(A, b);  % R^2 -> R^2
    %
    %   % Evaluate
    %   X = rand(2, 100);
    %   Y = f.evaluate(X);
    %
    % See Also:
    %   core.function.Function, core.function.ConstantFunction

    properties
        A % Transformation matrix (nCodims x nDims)
        b % Translation vector (nCodims x 1)
    end

    properties (Dependent)
        isAffine % True if the function includes a translation (b ~= 0)
        isPureLinear % True if the function is purely linear (b == 0)
    end

    methods
        function obj = LinearFunction(A, b)
            % LINEARFUNCTION Constructor for the LinearFunction class.
            %
            %   obj = LinearFunction(A) creates a pure linear
            %   transformation f(x) = A*x.
            %
            %   obj = LinearFunction(A, b) creates an affine transformation
            %   f(x) = A*x + b.
            %
            % Inputs:
            %   A - Transformation matrix (nCodims x nDims)
            %   b - Translation vector (nCodims x 1, optional)
            %
            % Outputs:
            %   obj - The created LinearFunction object

            core.except.assert(isnumeric(A) && ismatrix(A), ...
                'InvalidInput', ...
                'A must be a numeric matrix.');

            [nCodims, nDims] = size(A);

            core.except.assert(nCodims > 0 && nDims > 0, ...
                'InvalidInput', ...
                'A must be non-empty.');

            if nargin >= 2 && ~isempty(b)
                core.except.assert(isnumeric(b), ...
                    'InvalidInput', ...
                    'b must be numeric.');

                b = b(:);
                core.except.assert(length(b) == nCodims, ...
                    'InvalidInput', ...
                    'b must have the same number of rows as A.');
            else
                b = zeros(nCodims, 1);
            end

            obj@core.function.Function(nDims, nCodims);
            obj.A = A;
            obj.b = b;
        end

        function TF = get.isAffine(obj)
            % GET.ISAFFINE Check if the function is affine.

            TF = any(obj.b ~= 0);
        end

        function TF = get.isPureLinear(obj)
            % GET.ISPURELINEAR Check if the function is purely linear.

            TF = all(obj.b == 0);
        end
    end

    methods (Access = protected)
        function Y = evaluateImpl(obj, X)
            % EVALUATEIMPL Implementation of linear function evaluation.
            %
            %   Y = evaluateImpl(obj, X) computes the affine transformation
            %   Y = A*X + b for each column of X.
            %
            % Inputs:
            %   obj - The LinearFunction object
            %   X - Coordinates (nDims x nPoints)
            %
            % Outputs:
            %   Y - Transformed values (nCodims x nPoints)

            nPoints = size(X, 2);
            Y = obj.A * X;

            if obj.isAffine
                Y = Y + repmat(obj.b, 1, nPoints);
            end
        end

        function dY = derivativeImpl(obj, X, r)
            % DERIVATIVEIMPL Implementation of linear function derivative.
            %
            %   dY = derivativeImpl(obj, X, r) computes derivatives of the
            %   linear function. For first-order derivatives, this returns
            %   the appropriate columns/rows of matrix A. Higher-order
            %   derivatives are zero since linear functions have constant
            %   first derivatives.
            %
            % Inputs:
            %   obj - The LinearFunction object
            %   X - Coordinates (nDims x nPoints)
            %   r - Derivative order vector (1 x nDims)
            %
            % Outputs:
            %   dY - Derivative values (nCodims x nPoints)

            nPoints = size(X, 2);

            if sum(r) == 1
                dY = repmat(obj.A(:, find(r > 0, 1)), 1, nPoints);
            else
                dY = zeros(obj.nCodims, nPoints);
            end
        end
    end
end