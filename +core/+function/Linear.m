classdef Linear < core.function.Function
    % LINEAR Function representing linear/affine transformations.
    %
    %   Linear represents linear/affine transformations of the form:
    %
    %   \f[
    %     f(x) = A*x + b
    %   \f]
    %
    %   where \f$A\f$ is an \f$m\times n\f$ matrix and \f$b\f$ is an
    %   \f$m\times 1\f$ vector. When \f$b\f$ is not provided or is zero,
    %   this reduces to a pure linear transformation.
    %
    % See also:
    %   core.function.Function, core.function.Constant

    properties
        A % Transformation matrix (nCodims x nDims)
        b % Translation vector (nCodims x 1)
    end

    methods
        function obj = Linear(A, options)
            % LINEAR Constructor for the Linear class.
            %
            %   obj = Linear(A) creates a pure linear transformation
            %   \f$f(x) = Ax\f$.
            %
            %   obj = Linear(A, b=b) creates an affine transformation
            %   \f$f(x) = Ax + b\f$.

            arguments
                A {mustBeNumeric, mustBeNonempty}
                options.b {mustBeNumeric} = []
            end

            core.except.assert(ismatrix(A), 'InvalidInput', ...
                'A must be a numeric matrix.');

            [nCodims, nDims] = size(A);

            core.except.assert(isempty(options.b) || length(options.b) == nCodims, ...
                'InvalidInput', 'b must have the same number of rows as A.');

            obj@core.function.Function(nDims=nDims, nCodims=nCodims);
            obj.A = A;
            obj.b = options.b;
        end
    end

    methods (Access = protected)
        function Y = evalImpl(obj, X)
            % EVALIMPL Implementation of linear function evaluation.

            np = size(X, 2);
            Y = obj.A * X;

            if ~isempty(obj.b)
                Y = Y + repmat(obj.b(:), 1, np);
            end
        end

        function dY = diffImpl(obj, X, order)
            % DIFFIMPL Implementation of linear function derivative.

            np = size(X, 2);

            if sum(order) == 1
                dY = repmat(obj.A(:, find(order > 0, 1)), 1, np);
            else
                dY = zeros(obj.NCodims, np);
            end
        end
    end
end