classdef ProductFunction < core.function.Function
    % PRODUCTFUNCTION Function representing Hadamard product of two functions.
    %
    %   ProductFunction represents the element-wise (Hadamard) product of
    %   two functions with broadcasting support:
    %     (f * g)(x) = f(x) .* g(x)
    %
    %   The class supports broadcasting following MATLAB's element-wise
    %   operation rules. Input dimensions must be identical, while output
    %   dimensions follow broadcasting compatibility rules. The class
    %   implements the product rule for derivatives: (f * g)' = f' * g + f * g'
    %
    % Examples:
    %   % Compatible functions with broadcasting
    %   f = SomeFunction(3, 5);  % R^3 -> R^5
    %   g = SomeFunction(3, 1);  % R^3 -> R^1 (broadcasts to R^5)
    %
    %   % Create product function
    %   prod_func = ProductFunction(f, g);  % R^3 -> R^5
    %
    %   % Evaluate
    %   X = rand(3, 100);
    %   Y = prod_func.evaluate(X);  % Result: R^3 -> R^5
    %
    % See Also:
    %   core.function.Function

    properties
        f % First function
        g % Second function
    end

    methods
        function obj = ProductFunction(f, g)
            % PRODUCTFUNCTION Constructor for the ProductFunction class.
            %
            %   obj = ProductFunction(f, g) creates a ProductFunction
            %   representing the Hadamard product f .* g of two functions.
            %   Input dimensions must be identical, while output dimensions
            %   must be broadcasting-compatible.
            %
            % Inputs:
            %   f - First function (core.function.Function)
            %   g - Second function (core.function.Function)
            %
            % Outputs:
            %   obj - The created ProductFunction object

            core.except.assert(isa(f, 'core.function.Function'), ...
                'InvalidInput', ...
                'First argument must be a Function object.');

            core.except.assert(isa(g, 'core.function.Function'), ...
                'InvalidInput', ...
                'Second argument must be a Function object.');

            core.except.assert(f.nDims == g.nDims, ...
                'InputDimensionMismatch', ...
                'Input dimensions must be identical. Got f: %d, g: %d.', ...
                f.nDims, g.nDims);

            compatible = (f.nCodims == g.nCodims) || (f.nCodims == 1) || (g.nCodims == 1);

            core.except.assert(compatible, ...
                'BroadcastIncompatible', ...
                'Output dimensions not broadcasting-compatible: %d and %d', ...
                f.nCodims, g.nCodims);

            obj@core.function.Function(f.nDims, max(f.nCodims, g.nCodims));
            obj.f = f;
            obj.g = g;
        end
    end

    methods (Access = protected)
        function Y = evaluateImpl(obj, X)
            % EVALUATEIMPL Implementation of product function evaluation.
            %
            %   Y = evaluateImpl(obj, X) evaluates the Hadamard product
            %   f(X) .* g(X) with broadcasting support for compatible
            %   output dimensions.
            %
            % Inputs:
            %   obj - The ProductFunction object
            %   X - Coordinates (nDims x nPoints)
            %
            % Outputs:
            %   Y - Function values (nCodims x nPoints)

            Y1 = obj.f.evaluate(X);
            Y2 = obj.g.evaluate(X);
            Y = Y1 .* Y2;
        end

        function dY = derivativeImpl(obj, X, r)
            % DERIVATIVEIMPL Implementation of product function derivative.
            %
            %   dY = derivativeImpl(obj, X, r) computes the derivative of
            %   the product function using the product rule:
            %   (f * g)' = f' * g + f * g'
            %
            %   For higher-order derivatives, this uses the generalized
            %   Leibniz rule (binomial expansion).
            %
            % Inputs:
            %   obj - The ProductFunction object
            %   X - Coordinates (nDims x nPoints)
            %   r - Derivative order vector (1 x nDims)
            %
            % Outputs:
            %   dY - Derivative values (nCodims x nPoints)

            derivativeOrder = sum(r);

            if derivativeOrder == 1
                F = obj.f.evaluate(X);
                G = obj.g.evaluate(X);
                dF = obj.f.derivative(X, r);
                dG = obj.g.derivative(X, r);
                dY = dF .* G + F .* dG;
            else
                n = length(r);
                dY = zeros(obj.nCodims, size(X, 2));

                indices = cell(1, n);
                for i = 1:n
                    indices{i} = 0:r(i);
                end

                grids = cell(1, n);
                [grids{:}] = ndgrid(indices{:});

                m = prod(cellfun(@numel, indices));
                R1 = zeros(m, n);

                for i = 1:n
                    R1(:, i) = grids{i}(:);
                end

                R2 = repmat(r, m, 1) - R1;

                for i = 1:m
                    r1 = R1(i, :);
                    r2 = R2(i, :);

                    coeff = 1;
                    for j = 1:length(r)
                        if r(j) > 0
                            coeff = coeff * nchoosek(r(j), r1(j));
                        end
                    end

                    if all(r1 == 0)
                        df = obj.f.evaluate(X);
                    else
                        df = obj.f.derivative(X, r1);
                    end

                    if all(r2 == 0)
                        dg = obj.g.evaluate(X);
                    else
                        dg = obj.g.derivative(X, r2);
                    end

                    dY = dY + coeff * (df .* dg);
                end
            end
        end
    end
end