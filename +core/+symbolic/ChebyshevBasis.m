classdef ChebyshevBasis < core.symbolic.OrthogonalPolynomialBasis
    % CHEBYSHEVBASIS Chebyshev orthogonal basis functions.
    %
    %   ChebyshevBasis generates Chebyshev polynomial basis functions of
    %   the first kind up to a specified maximum degree. The basis
    %   functions are orthogonal on the interval [-1, 1] with weight
    %   function \f$ 1/\sqrt{1-x^2} \f$.
    %
    %   The Chebyshev polynomials of the first kind satisfy the recurrence:
    %
    %   \f[
    %       T_{n+1}(x) = 2x \cdot T_n(x) - T_{n-1}(x)
    %   \f]
    %
    % See also:
    %   core.symbolic.OrthogonalPolynomialBasis,
    %   core.symbolic.LegendreBasis

    properties
        Lower % Lower bound
        Upper % Upper bound
    end

    methods
        function obj = ChebyshevBasis(variables, maxDegree, options)
            % CHEBYSHEVBASIS Construct Chebyshev basis functions.
            %
            %   obj = ChebyshevBasis(variables, maxDegree) creates
            %   Chebyshev basis functions using the specified @a variables
            %   up to @a maxDegree on the domain \f$[-1, 1]\f$.
            %
            %   obj = ChebyshevBasis(variables, maxDegree,
            %   isNormalized=true) creates orthonormal Chebyshev basis
            %   functions.
            %
            %   obj = ChebyshevBasis(variables, maxDegree, lower=a,
            %   upper=b) creates basis functions on the domain \f$[a,
            %   b]\f$.

            arguments
                variables {mustBeA(variables, 'sym')}
                maxDegree {mustBeNonnegative, mustBeInteger}
                options.isNormalized {mustBeNumericOrLogical} = false
                options.isMonic {mustBeNumericOrLogical} = false
                options.lower {mustBeNumeric, mustBeFinite} = -1
                options.upper {mustBeNumeric, mustBeFinite} = 1
            end

            obj = obj@core.symbolic.OrthogonalPolynomialBasis( ...
                variables, maxDegree, options);
        end
    end

    methods (Access = protected)
        function obj = initialize(obj, options)
            % INITIALIZE Initialize optional properties.

            initialize@core.symbolic.OrthogonalPolynomialBasis(obj, options);
            obj.Lower = options.lower;
            obj.Upper = options.upper;
        end

        function [a, b, c] = computeCoefficients(~, n)
            % COMPUTECOEFFICIENTS Compute Chebyshev recurrence
            % coefficients.

            if n == 0
                a = 1;
                b = 0;
                c = 0;
            else
                a = 2;
                b = 0;
                c = -1;
            end
        end

        function normalizer = computeNormalizer(obj, n)
            % COMPUTENORMALIZER Compute Chebyshev normalizer.

            l = obj.Upper - obj.Lower;
            if n == 0
                normalizer = 1 / sqrt(sym(pi)) / sqrt(l);
            else
                normalizer = sqrt(2/sym(pi)) / sqrt(l);
            end
        end

        function x = transform(obj, x)
            % TRANSFORM Transform variable to polynomial's canonical domain.
            %
            %   x = transform(x) transforms @a x to the canonical domain
            %   for the specific polynomial family.

            l = obj.Upper - obj.Lower;
            c = (obj.Upper + obj.Lower) / 2;
            x = (2 * x - c) / l;
        end
    end
end