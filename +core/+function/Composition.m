classdef Composition < core.function.Function
    % COMPOSITION Function representing the composition of two functions.
    %
    %   Composition represents the mathematical composition of two
    %   functions:
    %
    %   \f[
    %       (f \circ g)(x) = f(g(x))
    %   \f]
    %
    %   where the output dimension of the inner function @a rhs must match
    %   the input dimension of the outer function @a lhs. The composition
    %   maps from the input space of @a rhs to the output space of @a lhs.
    %
    %   For derivatives, the chain rule is applied:
    %   \f[
    %       \frac{d}{dx}(f \circ g)(x) = f'(g(x)) \cdot g'(x)
    %   \f]
    %
    %   Higher-order derivatives use Faà di Bruno's formula for the
    %   composition of functions.
    %
    % See also:
    %   core.function.Function, core.function.ProductFunction

    properties
        Lhs % Outer function (core.function.Function)
        Rhs % Inner function (core.function.Function)
    end

    methods
        function obj = Composition(lhs, rhs)
            % COMPOSITION Constructor for the Composition class.
            %
            %   obj = Composition(f, g) creates a Composition representing
            %   the composition \f$f(g(x))\f$ where @a lhs is the outer
            %   function \f$f\f$ and @a rhs is the inner function \f$g\f$.
            %   The output dimension of @a rhs must match the input
            %   dimension of @a lhs.

            arguments
                lhs{mustBeA(lhs, "core.function.Function")}
                rhs{mustBeA(rhs, "core.function.Function")}
            end

            core.except.assert(rhs.NCodims == lhs.NDims, ...
                'DimensionMismatch', ...
                ['Output dimension of inner function (%d) must ' ...
                'match input dimension of outer function (%d).'], ...
                rhs.NCodims, lhs.NDims);

            nDims = rhs.NDims;
            nCodims = lhs.NCodims;

            obj@core.function.Function(nDims = nDims, nCodims = nCodims);
            obj.Lhs = lhs;
            obj.Rhs = rhs;
        end
    
        function J = jacobian(obj, X)
            % JACOBIAN Compute the Jacobian matrix at specified points.
            %
            %   J = jacobian(obj, X) computes the Jacobian matrix of the
            %   function at the specified points. For a function from
            %   \f$R^n\f$ to \f$R^m\f$, the Jacobian is an \f$m\times n\f$
            %   matrix of first partial derivatives.

            arguments
                obj core.function.Function
                X {mustBeNumeric}
            end

            core.except.assert(obj.IsWellDefined, 'NotWellDefined', ...
                'Function is not well-defined.');

            nd = obj.NDims;
            core.except.assert(size(X, 1) == nd, 'InvalidInput', ...
                'Input dimension mismatch.');

            nx = size(X);
            X = reshape(X, nd, []);
            
            nc = obj.NCodims;

            G = obj.Rhs.eval(X);
            JF = obj.Lhs.jacobian(G);
            JG = obj.Rhs.jacobian(X);
            J = pagemtimes(JF, JG);

            if length(nx) > 2
                J = reshape(J, [nc, nd, nx(2:end)]);
            end
        end

    end

    methods (Access = protected)
        function Y = evalImpl(obj, X)
            % EVALIMPL Implementation of composition function evaluation.
            %
            %   Y = evalImpl(obj, X) evaluates the composition
            %   f(g(X)) by first evaluating the inner function g(X) and
            %   then evaluating the outer function f at those points.

            G = obj.Rhs.eval(X);
            Y = obj.Lhs.eval(G);
        end

        function dY = diffImpl(obj, X, order)
            % DIFFIMPL Implementation of composition function derivative.
            %
            %   dY = diffImpl(obj, X, order) computes the derivative of
            %   the composition using the chain rule for first-order
            %   derivatives and Faà di Bruno's formula for higher orders.

            core.except.assert(sum(order) == 1, 'NotImplemented', ...
                'Higher-order derivatives not yet implemented for Composition.');

            G = obj.Rhs.eval(X);
            dF = obj.Lhs.diff(G, order);
            dG = obj.Rhs.diff(X, order);
            dG = reshape(dG, obj.Rhs.NCodims, 1, []);
            dY = pagemtimes(dF, dG);
        end
    end
end