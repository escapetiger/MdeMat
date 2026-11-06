classdef ChebyshevBasis < mdemat.common.spectral.OrthogonalPolynomialBasis
    % CHEBYSHEVBASIS Chebyshev orthogonal basis functions.
    %
    %   ChebyshevBasis generates Chebyshev polynomial basis functions of
    %   the first kind up to a specified maximum degree. The basis
    %   functions are orthogonal on the interval \f$[-1, 1]\f$ with weight
    %   function \f$ 1/\sqrt{1-x^2} \f$.
    %
    %   The Chebyshev polynomials of the first kind satisfy the recurrence:
    %
    %   \f[
    %       T_{n+1}(x) = 2x \cdot T_n(x) - T_{n-1}(x)
    %   \f]
    %
    % See also:
    %   mdemat.common.spectral.OrthogonalPolynomialBasis,
    %   mdemat.common.spectral.LegendreBasis

    methods (Access = protected)
        function Y = first(obj, X)
            % FIRST Compute the first term.
            
            if obj.IsNormalized
                Y = repmat(1/sqrt(pi), 1, size(X, 2));
            else
                Y = ones(1, size(X, 2));
            end
        end

        function [a, b, c] = computeCoefficients(obj, n)
            % COMPUTECOEFFICIENTS Compute Chebyshev recurrence
            % coefficients.

            if obj.IsNormalized
                if n == 0
                    a = sqrt(2);
                    b = 0;
                    c = 0;
                elseif n == 1
                    a = 2;
                    b = 0;
                    c = -sqrt(2);
                else
                    a = 2;
                    b = 0;
                    c = -1;
                end
            else
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
        end
    end
end
