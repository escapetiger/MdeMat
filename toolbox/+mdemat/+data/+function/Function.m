classdef Function < handle
    % FUNCTION Base class for all functions from Euclidean space to
    % Euclidean space.
    %
    %   Function represents a mathematical mapping from \f$R^n\f$ to
    %   \f$R^m\f$, which can be evaluated at specified points. This class
    %   defines the interface for all function types in the system.

    properties
        NDims % Number of dimensions of input space
        NCodims % Number of dimensions of output space
    end

    properties (Dependent)
        IsWellDefined % True if both NDims and NCodims are positive
    end

    methods
        function obj = Function(options)
            % FUNCTION Constructor for the abstract Function class.
            %
            %   obj = Function() creates a Function with default dimensions
            %   (both nDims and nCodims are 0).
            %
            %   obj = Function(nDims=nDims, nCodims=nCodims) creates a
            %   Function with specified input dimension @a nDims and output
            %   dimension @a nCodims.

            arguments
                options.nDims (1,1) {mustBeInteger} = 0
                options.nCodims (1,1) {mustBeInteger} = 0
            end

            obj.NDims = options.nDims;
            obj.NCodims = options.nCodims;
        end

        function Y = eval(obj, X)
            % EVAL Evaluate the function at specified points.
            %
            %   Y = eval(obj, X) evaluates the function at the points
            %   specified in @a X. The input @a X should have NDims rows, and
            %   the output @a Y will have NCodims rows.

            arguments
                obj mdemat.data.function.Function
                X {mustBeNumeric}
            end

            assert(all(obj.IsWellDefined, 'all'), 'mdemat:data:function:Function:NotWellDefined', ...
                '%s', 'Function is not well-defined.');

            nd = obj.NDims;
            assert(all(size(X, 1) == nd, 'all'), 'mdemat:data:function:Function:InvalidInput', ...
                '%s', 'Input dimension mismatch.');

            nx = size(X);
            X = reshape(X, nd, []);
            Y = obj.evalImpl(X);

            nc = obj.NCodims;
            assert(all(size(Y, 1) == nc, 'all'), 'mdemat:data:function:Function:InvalidOutput', ...
                '%s', 'Output dimension mismatch.');

            if length(nx) > 2, Y = reshape(Y, [nc, nx(2:end)]); end
        end

        function dY = diff(obj, X, order)
            % DIFF Compute the derivative of specified order at
            % given points.
            %
            %   dY = diff(obj, X, order) computes the derivative
            %   of order specified by the vector @a order at the points @a
            %   X. The derivative order vector @a order specifies the
            %   partial derivative order for each input dimension.

            arguments
                obj mdemat.data.function.Function
                X {mustBeNumeric}
                order {mustBeVector, mustBeNonnegative}
            end

            assert(all(obj.IsWellDefined, 'all'), ...
                'mdemat:data:function:Function:NotWellDefined', '%s', 'Function is not well-defined.');

            nd = obj.NDims;
            isValidOrder = length(order) == nd && sum(order) > 0;
            assert(all(isValidOrder, 'all'), 'mdemat:data:function:Function:InvalidInput', ...
                '%s', 'Bad derivative order.');

            assert(all(size(X, 1) == nd, 'all'), 'mdemat:data:function:Function:InvalidInput', ...
                '%s', 'Input dimension mismatch.');

            nx = size(X);
            X = reshape(X, nd, []);
            dY = obj.diffImpl(X, order);

            nc = obj.NCodims;
            assert(all(size(dY, 1) == nc, 'all'), 'mdemat:data:function:Function:InvalidOutput', ...
                '%s', 'Output dimension mismatch.');

            if length(nx) > 2
                dY = reshape(dY, [nc, nx(2:end)]);
            end
        end

        function J = jacobian(obj, X)
            % JACOBIAN Compute the Jacobian matrix at specified points.
            %
            %   J = jacobian(obj, X) computes the Jacobian matrix of the
            %   function at the specified points. For a function from
            %   \f$R^n\f$ to \f$R^m\f$, the Jacobian is an \f$m\times n\f$
            %   matrix of first partial derivatives.

            arguments
                obj mdemat.data.function.Function
                X {mustBeNumeric}
            end

            assert(all(obj.IsWellDefined, 'all'), 'mdemat:data:function:Function:NotWellDefined', ...
                '%s', 'Function is not well-defined.');

            nd = obj.NDims;
            assert(all(size(X, 1) == nd, 'all'), 'mdemat:data:function:Function:InvalidInput', ...
                '%s', 'Input dimension mismatch.');

            nx = size(X);
            X = reshape(X, nd, []);

            np = size(X, 2);
            nc = obj.NCodims;
            J = zeros(nc, nd, np);

            for i = 1:nc
                for j = 1:nd
                    order = zeros(1, nd);
                    order(j) = 1;
                    if np > 0
                        dYj = obj.diffImpl(X, order);
                        if size(dYj, 1) == 1
                            J(i, j, :) = dYj;
                        else
                            J(i, j, :) = dYj(i, :);
                        end
                    end
                end
            end

            if length(nx) > 2
                J = reshape(J, [nc, nd, nx(2:end)]);
            end
        end

        function G = gradient(obj, X)
            % GRADIENT Compute the gradient at specified points.
            %
            %   G = gradient(obj, X) computes the gradient vector of a
            %   scalar-valued function at the specified points.
            %
            %   The gradient is only defined for scalar-valued functions.

            arguments
                obj mdemat.data.function.Function
                X {mustBeNumeric}
            end

            nd = obj.NDims;
            nc = obj.NCodims;
            assert(all(nc == 1, 'all'), 'mdemat:data:function:Function:NotScalarValued', ...
                '%s', 'Gradient is only defined for scalar-valued functions.');

            J = obj.jacobian(X);

            nx = size(X);
            if length(nx) > 2
                G = reshape(squeeze(J(1, :, :)), [nd, nx(2:end)]);
            else
                G = squeeze(J(1, :, :));
            end
        end

        function TF = get.IsWellDefined(obj)
            % GET.ISWELLDEFINED Check if a function is well-defined.

            TF = (obj.NDims > 0) && (obj.NCodims > 0);
        end
    end

    methods (Access = protected)
        function Y = evalImpl(obj, X)
            % EVALIMPL Implementation of function evaluation.

            Y = X;
            assert(all(0, 'all'), 'mdemat:data:function:Function:NotImplemented', ...
                'Method evalImpl is not implemented by class %s.', class(obj));
        end

        function dY = diffImpl(obj, X, ~)
            % DIFFIMPL Implementation of function derivative
            % evaluation.

            dY = X;
            assert(all(0, 'all'), 'mdemat:data:function:Function:NotImplemented', ...
                'Method diffImpl is not implemented by class %s.', class(obj));
        end
    end
end