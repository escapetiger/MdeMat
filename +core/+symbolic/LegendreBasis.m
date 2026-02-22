classdef LegendreBasis < core.symbolic.OrthogonalPolynomialBasis
    % LEGENDREBASIS Legendre orthogonal basis functions.
    %
    %   LegendreBasis generates Legendre polynomial basis functions up to
    %   a specified maximum degree. The basis functions are orthogonal on
    %   the interval [-1, 1] with uniform weight function.
    %
    %   The Legendre polynomials satisfy the three-term recurrence
    %   relation:
    %
    %   \f[
    %       P_{n+1}(x) = \frac{(2n+1)x \cdot P_n(x) - n \cdot
    %       P_{n-1}(x)}{n+1}
    %   \f]
    %
    % See also:
    %   core.symbolic.OrthogonalPolynomialBasis,
    %   core.symbolic.ChebyshevBasis

    properties
        Lower % Lower bound
        Upper % Upper bound
    end

    methods
        function obj = LegendreBasis(variables, maxDegree, options)
            % LEGENDREBASIS Construct Legendre basis functions.
            %
            %   obj = LegendreBasis(variables, maxDegree) creates Legendre
            %   basis functions using the specified @a variables up to @a
            %   maxDegree on the domain \f$[-1, 1]\f$.
            %
            %   obj = LegendreBasis(variables, maxDegree,
            %   isNormalized=true) creates orthonormal Legendre basis
            %   functions.
            %
            %   obj = LegendreBasis(variables, maxDegree, lower=a, upper=b)
            %   creates basis functions on the domain \f$[a, b]\f$.

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
            % COMPUTECOEFFICIENTS Compute Legendre recurrence coefficients.

            a = (2 * n + 1) / (n + 1);
            b = 0;
            c = -n / (n + 1);
        end

        function normalizer = computeNormalizer(obj, n)
            % COMPUTENORMALIZER Compute Legendre normalizer.

            normalizer = sqrt(2*n+1) / sqrt(obj.Upper-obj.Lower);
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