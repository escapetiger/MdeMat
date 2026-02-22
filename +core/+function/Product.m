classdef Product < core.function.Function
    % PRODUCT Function representing Hadamard product of two functions.
    %
    %   Product represents the element-wise (Hadamard) product of
    %   two functions with broadcasting support:
    %
    %   \[
    %       (f * g)(x) = f(x) .* g(x)
    %   \]
    %
    %   The class supports broadcasting following MATLAB's element-wise
    %   operation rules. Input dimensions must be identical, while output
    %   dimensions follow broadcasting compatibility rules. The class
    %   implements the product rule for derivatives: 
    % 
    %   \[
    %       (f * g)' = f' * g + f * g'
    %   \]
    %
    % See Also:
    %   core.function.Function

    properties
        Lhs % LHS function
        Rhs % RHSfunction
    end

    methods
        function obj = Product(lhs, rhs)
            % PRODUCT Constructor for the Product class.
            %
            %   obj = Product(f, g) creates a Product representing the
            %   Hadamard product @a f .* @a g of two functions. Input
            %   dimensions must be identical, while output dimensions must
            %   be broadcasting-compatible.

            arguments
                lhs {mustBeA(lhs, "core.function.Function")}
                rhs {mustBeA(rhs, "core.function.Function")}
            end

            core.except.assert(lhs.NDims == rhs.NDims, ...
                'InputDimensionMismatch', ...
                'Input dimensions must be identical. Got f: %d, g: %d.', ...
                lhs.NDims, rhs.NDims);

            compatible = (lhs.NCodims == rhs.NCodims) || (lhs.NCodims == 1) || (rhs.NCodims == 1);

            core.except.assert(compatible, 'BroadcastIncompatible', ...
                'Output dimensions not broadcasting-compatible: %d and %d', ...
                lhs.NCodims, rhs.NCodims);

            nDims = lhs.NDims;
            nCodims = max(lhs.NCodims, rhs.NCodims);

            obj@core.function.Function(nDims=nDims, nCodims=nCodims);
            obj.Lhs = lhs;
            obj.Rhs = rhs;
        end
    end

    methods (Access = protected)
        function Y = evalImpl(obj, X)
            % EVALIMPL Implementation of product function evaluation.
            %
            %   Y = evalImpl(obj, X) evaluates the Hadamard product
            %   f(X) .* g(X) with broadcasting support for compatible
            %   output dimensions.

            Y1 = obj.Lhs.eval(X);
            Y2 = obj.Rhs.eval(X);
            Y = Y1 .* Y2;
        end

        function dY = diffImpl(obj, X, order)
            % DIFFIMPL Implementation of product function
            % derivative.
            %
            %   dY = diffImpl(obj, X, order) computes the
            %   derivative of the product function using the product rule:
            %   (f * g)' = f' * g + f * g'
            %
            %   For higher-order derivatives, this uses the generalized
            %   Leibniz rule (binomial expansion).

            if sum(order) == 1
                F = obj.Lhs.eval(X);
                G = obj.Rhs.eval(X);
                dF = obj.Lhs.diff(X, order);
                dG = obj.Rhs.diff(X, order);
                dY = dF .* G + F .* dG;
            else
                nd = length(order);
                np = size(X, 2);
                dY = zeros(obj.NCodims, np);

                indices = cell(1, nd);
                for i = 1:nd
                    indices{i} = 0:order(i);
                end

                grids = cell(1, nd);
                [grids{:}] = ndgrid(indices{:});

                m = prod(cellfun(@numel, indices));
                R1 = zeros(m, nd);

                for i = 1:nd
                    R1(:, i) = grids{i}(:);
                end

                R2 = repmat(order, m, 1) - R1;

                for i = 1:m
                    r1 = R1(i, :);
                    r2 = R2(i, :);

                    coeff = 1;
                    for j = 1:length(order)
                        if order(j) > 0
                            coeff = coeff * nchoosek(order(j), r1(j));
                        end
                    end

                    if all(r1 == 0)
                        dF = obj.Lhs.eval(X);
                    else
                        dF = obj.Lhs.diff(X, r1);
                    end

                    if all(r2 == 0)
                        dG = obj.Rhs.eval(X);
                    else
                        dG = obj.Rhs.diff(X, r2);
                    end

                    dY = dY + coeff * (dF .* dG);
                end
            end
        end
    end
end