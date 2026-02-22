classdef HermiteBasis < core.symbolic.OrthogonalPolynomialBasis
    % HERMITEBASIS Hermite orthogonal basis functions.
    %
    %   HermiteBasis generates probabilist's Hermite polynomial basis
    %   functions up to a specified maximum degree. The basis functions are
    %   orthogonal on the interval (-∞, ∞) with weight function
    %   \f$ e^{-x^2/2}/\sqrt{2\pi} \f$.
    %
    %   The probabilist's Hermite polynomials satisfy the recurrence:
    %
    %   \f[
    %       H_{n+1}(x) = x \cdot H_n(x) - n \cdot H_{n-1}(x)
    %   \f]
    %
    % See also:
    %   core.symbolic.OrthogonalPolynomialBasis,
    %   core.symbolic.LegendreBasis

    methods
        function obj = HermiteBasis(variables, maxDegree, options)
            % HERMITEBASIS Construct Hermite basis functions.
            %
            %   obj = HermiteBasis(variables, maxDegree) creates Hermite
            %   basis functions using the specified @a variables up to @a
            %   maxDegree on the unbounded domain \f$(-\infty, \infty)\f$.
            %
            %   obj = HermiteBasis(variables, maxDegree, isNormalized=true)
            %   creates orthonormal Hermite basis functions.

            arguments
                variables {mustBeA(variables, 'sym')}
                maxDegree {mustBeNonnegative, mustBeInteger}
                options.isNormalized {mustBeNumericOrLogical} = false
                options.isMonic {mustBeNumericOrLogical} = false
            end

            obj = obj@core.symbolic.OrthogonalPolynomialBasis( ...
                variables, maxDegree, options);
        end
    end

    methods (Access = protected)
        function [a, b, c] = computeCoefficients(~, n)
            % COMPUTECOEFFICIENTS Compute Hermite recurrence coefficients.

            a = 1;
            b = 0;
            c = -n;
        end

        function normalizer = computeNormalizer(~, n)
            % COMPUTENORMALIZER Compute Hermite normalizer.

            normalizer = 1 / sqrt(factorial(n));
        end

        function x = transform(~, x)
            % TRANSFORM Transform variable to polynomial's canonical domain.
            %
            %   x = transform(x) transforms @a x to the canonical domain
            %   for the specific polynomial family.
            
            % For Hermite polynomials, no transformation is needed 
            % as they are defined on the entire real line
        end
    end
end