classdef HermiteBasis < core.function.OrthogonalPolynomialBasis
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
    %   core.function.OrthogonalPolynomialBasis,
    %   core.function.LegendreBasis

    methods (Access = protected)
        function Y = first(~, X)
            % FIRST Compute the first term.

            Y = ones(1, size(X, 2));
        end

        function [a, b, c] = computeCoefficients(obj, n)
            % COMPUTECOEFFICIENTS Compute Hermite recurrence coefficients.

            if obj.IsNormalized
                a = 1 / sqrt(n+1);
                b = 0;
                c = -sqrt(n/(n + 1));
            else
                a = 1;
                b = 0;
                c = -n;
            end
        end
    end
end