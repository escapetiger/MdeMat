classdef Function < handle
    % FUNCTION Base class for all functions from Euclidean space to
    % Euclidean space.
    %
    %   Function represents a mathematical mapping from R^n to R^m, which
    %   can be evaluated at specified points. This class defines the
    %   interface for all function types in the system.
    %
    %   The implementation uses a two-step evaluation process:
    %   1. The public evaluate method handles input/output reshaping
    %   2. The protected evaluateImpl method handles the actual computation

    properties
        nDims % Number of dimensions of input space
        nCodims % Number of dimensions of output space
        parameters % Parameter registry
    end

    properties (Dependent)
        isWellDefined % True if both input and output dimensions are positive
        hasParameters % True if the function has registered parameters
    end

    methods
        function obj = Function(nDims, nCodims)
            % FUNCTION Constructor for the abstract Function class.
            %
            %   obj = Function() creates a Function with default dimensions
            %   (both nDims and nCodims are 0).
            %
            %   obj = Function(nDims) creates a Function with specified
            %   input dimension and default output dimension (0).
            %
            %   obj = Function(nDims, nCodims) creates a Function with
            %   specified input and output dimensions.
            %
            % Inputs:
            %   nDims - Number of input dimensions (optional, default:0)
            %   nCodims - Number of output dimensions (optional, default: 0)
            %
            % Outputs:
            %   obj - The created Function object

            if nargin >= 1 && ~isempty(nDims)
                obj.nDims = nDims;
            else
                obj.nDims = 0;
            end

            if nargin >= 2 && ~isempty(nCodims)
                obj.nCodims = nCodims;
            else
                obj.nCodims = 0;
            end
            obj.parameters = struct();
        end

        function Y = evaluate(obj, X)
            % EVALUATE Evaluate the function at specified points.
            %
            %   Y = evaluate(obj, X) evaluates the function at the points
            %   specified in X. The input X should have nDims rows, and
            %   the output Y will have nCodims rows.
            %
            % Inputs:
            %   obj - The Function object
            %   X - Coordinates (nDims x ...)
            %
            % Outputs:
            %   Y - Function values (nCodims x ...)
            %
            % Examples:
            %   % Evaluate at single point
            %   x = [1; 2; 3];
            %   y = evaluate(obj, x);
            %
            %   % Evaluate at multiple points
            %   X = rand(3, 100);
            %   Y = evaluate(obj, X);

            core.except.assert(obj.isWellDefined, ...
                'NotWellDefined', 'Function is not well-defined.');

            n = obj.nDims;
            core.except.assert(size(X, 1) == n, ...
                'InvalidInput', 'Input dimension mismatch.');

            s = size(X);
            X = reshape(X, n, []);
            Y = obj.evaluateImpl(X);

            m = obj.nCodims;
            core.except.assert(size(Y, 1) == m, ...
                'InvalidOutput', 'Output dimension mismatch.');

            if length(s) > 2, Y = reshape(Y, [m, s(2:end)]); end
        end

        function dY = derivative(obj, X, r)
            % DERIVATIVE Compute the derivative of specified order at given
            % points.
            %
            %   dY = derivative(obj, X, r) computes the derivative of order
            %   specified by the vector r at the points X. The derivative
            %   order vector r specifies the partial derivative order for
            %   each input dimension.
            %
            % Inputs:
            %   obj - The DifferentiableFunction object
            %   X - Points where to evaluate (nDims x ...)
            %   r - Derivative order vector (1 x nDims)
            %
            % Outputs:
            %   dY - Derivative values (nCodims x ...)

            core.except.assert(obj.isWellDefined, ...
                'NotWellDefined', 'Function is not well-defined.');

            n = obj.nDims;
            core.except.assert(isvector(r) && isnumeric(r) ...
                && length(r) == n && all(r >= 0) && sum(r) > 0, ...
                'InvalidInput', 'Bad derivative order.');

            core.except.assert(size(X, 1) == n, ...
                'InvalidInput', 'Input dimension mismatch.');

            s = size(X);
            X = reshape(X, n, []);
            dY = obj.derivativeImpl(X, r);

            m = obj.nCodims;
            core.except.assert(size(dY, 1) == m, ...
                'InvalidOutput', 'Output dimension mismatch.');

            if length(s) > 2
                dY = reshape(dY, [m, s(2:end)]);
            end
        end

        function J = jacobian(obj, X)
            % JACOBIAN Compute the Jacobian matrix at specified points.
            %
            %   J = jacobian(obj, X) computes the Jacobian matrix of the
            %   function at the specified points. For a function from
            %   \f$R^n\f$ to \f$R^m\f$, the Jacobian is an \f$m\times n\f$
            %   matrix of first partial derivatives.
            %
            % Inputs:
            %   obj - The DifferentiableFunction object
            %   X - Coordinates (nDims x ...)
            %
            % Outputs:
            %   J - Jacobian tensor (nCodims x nDims x ...)

            core.except.assert(obj.isWellDefined, ...
                'NotWellDefined', 'Function is not well-defined.');

            n = obj.nDims;
            core.except.assert(size(X, 1) == n, ...
                'InvalidInput', 'Input dimension mismatch.');

            s = size(X);
            X = reshape(X, n, []);

            p = size(X, 2);
            m = obj.nCodims;
            J = zeros(m, n, p);

            for i = 1:m
                for j = 1:n
                    r = zeros(1, n);
                    r(j) = 1;
                    if p > 0
                        dYj = obj.derivativeImpl(X, r);
                        if size(dYj, 1) == 1
                            J(i, j, :) = dYj;
                        else
                            J(i, j, :) = dYj(i, :);
                        end
                    end
                end
            end

            if length(s) > 2
                J = reshape(J, [m, n, s(2:end)]);
            end
        end

        function G = gradient(obj, X)
            % GRADIENT Compute the gradient at specified points.
            %
            %   G = gradient(obj, X) computes the gradient vector of a
            %   scalar-valued function at the specified points.
            %
            % Notes:
            %   The gradient is only defined for scalar-valued functions.
            %
            % Inputs:
            %   obj - The DifferentiableFunction object (must be scalar-valued)
            %   X - Points where to evaluate (nDims x ...)
            %
            % Outputs:
            %   G - Gradient vector (nDims x ...)

            n = obj.nDims;
            m = obj.nCodims;
            core.except.assert(m == 1, 'NotScalarValued', ...
                'Gradient is only defined for scalar-valued functions.');

            J = obj.jacobian(X);

            s = size(X);
            if length(s) > 2
                G = reshape(squeeze(J(1, :, :)), [n, s(2:end)]);
            else
                G = squeeze(J(1, :, :));
            end
        end

        function TF = get.isWellDefined(obj)
            % GET.ISWELLDEFINED Check if a function is well-defined.

            TF = (obj.nDims > 0) && (obj.nCodims > 0);
        end

        function TF = get.hasParameters(obj)
            % GET.HASPARAMETERS Check if a function has parameters.

            TF = ~isempty(fieldnames(obj.parameters));
        end
    end

    methods (Abstract, Access = protected)
        % EVALUATEIMPL Implement the actual function evaluation.
        Y = evaluateImpl(obj, X)

        % DERIVATIVEIMPL Implement the actual derivative evaluation.
        dY = derivativeImpl(obj, X, r)
    end
end