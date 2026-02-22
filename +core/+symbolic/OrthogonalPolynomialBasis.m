classdef (Abstract) OrthogonalPolynomialBasis < core.symbolic.SymbolicFunction
    % ORTHOGONALPOLYNOMIALBASIS Abstract base class for orthogonal
    % polynomial basis functions.
    %
    %   OrthogonalPolynomialBasis provides a unified framework for
    %   generating orthogonal polynomial basis functions using three-term
    %   recurrence relations. Concrete subclasses implement specific
    %   polynomial families such as Legendre, Chebyshev, and Hermite
    %   polynomials.
    %
    %   The three-term recurrence relation has the general form:
    %
    %   \f[
    %       P_{n+1}(x) = (A_n \cdot x + B_n) \cdot P_n(x) + C_n \cdot
    %       P_{n-1}(x)
    %   \f]
    %
    % See also:
    %   core.symbolic.LegendreBasis, core.symbolic.ChebyshevBasis,
    %   core.symbolic.HermiteBasis

    properties
        MaxDegree % Maximum polynomial degree
        IsNormalized % Flag for orthonormal basis functions
        IsMonic % Flag for monic polynomials
    end

    methods
        function obj = OrthogonalPolynomialBasis(variables, maxDegree, options)
            % ORTHOGONALPOLYNOMIALBASIS Construct orthogonal polynomial
            % basis functions.
            %
            %   obj = OrthogonalPolynomialBasis(variables, maxDegree)
            %   creates orthogonal polynomial basis functions using the
            %   specified @a variables up to @a maxDegree.

            arguments
                variables {mustBeA(variables, 'sym')}
                maxDegree {mustBeNonnegative, mustBeInteger}
                options struct = struct()
            end

            obj = obj@core.symbolic.SymbolicFunction(variables);
            obj.MaxDegree = maxDegree;
            obj.initialize(options);
            obj.computeExpressions();
        end
    end

    methods (Access = protected)
        function obj = initialize(obj, options)
            % INITIALIZE Initialize optional properties.

            obj.IsNormalized = options.isNormalized;
            obj.IsMonic = options.isMonic;
        end

        function computeExpressions(obj)
            % COMPUTEEXPRESSIONS Generate polynomial basis function
            % expressions using recursion.

            obj.Expressions = repmat(sym([]), obj.MaxDegree+1, 1);
            x = obj.Variables(1);
            x = obj.transform(x);

            for n = 0:obj.MaxDegree
                p = obj.recursion(x, n);

                if obj.IsMonic
                    c = sym2poly(p);
                    p = p / c(1);
                end

                if obj.IsNormalized
                    normalizer = obj.computeNormalizer(n);
                    expr = simplify(normalizer*p);
                else
                    expr = simplify(p);
                end

                expr = horner(simplify(expr));
                obj.Expressions(n+1) = expr;
            end
        end

        function p = recursion(obj, x, n)
            % RECURSION Compute orthogonal polynomial using three-term
            % recurrence relation.
            %
            %   p = recursion(x, n) computes the nth orthogonal polynomial
            %   using the recurrence coefficients defined by subclasses.

            if n == 0
                p = sym(1);
            elseif n == 1
                [a0, b0, ~] = obj.computeCoefficients(0);
                p = a0 * x + b0;
            else
                [an1, bn1, cn1] = obj.computeCoefficients(n-1);
                p1 = obj.recursion(x, n-1);
                p2 = obj.recursion(x, n-2);
                p = (an1 * x + bn1) * p1 + cn1 * p2;
            end
        end

    end

    methods (Abstract, Access = protected)
        % COMPUTECOEFFICIENTS Compute coefficients for three-term
        % recurrence relation.
        %
        %   [a, b, c] = computeCoefficients(n) returns the coefficients
        %   @a a, @a b, and @a c for the three-term recurrence relation
        %   at degree @a n. The recurrence relation is: P_{n+1}(x) =
        %   (a*x + b)*P_n(x) + c*P_{n-1}(x)
        [a, b, c] = computeCoefficients(obj, n)

        % COMPUTENORMALIZER Compute normalizer for orthonormal basis.
        %
        %   factor = computeNormalizer(n) returns the normalization @a
        %   factor for the nth polynomial to create orthonormal basis
        %   functions.
        normalizer = computeNormalizer(obj, n)

        % TRANSFORM Transform variable to polynomial's canonical domain.
        %
        %   x = transform(x) transforms @a x to the canonical domain
        %   for the specific polynomial family.
        x = transform(obj, x)
    end
end