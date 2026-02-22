classdef LegendreBasis < core.function.OrthogonalPolynomialBasis
    % LEGENDREBASIS Legendre orthogonal basis functions.
    %
    %   LegendreBasis generates Legendre polynomial basis functions up to
    %   a specified maximum degree. The basis functions are orthogonal on
    %   the interval \f$[-1, 1]\f$ with uniform weight function.
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
    %   core.function.OrthogonalPolynomialBasis,
    %   core.function.ChebyshevBasis
    
    methods (Access = protected)
        function Y = first(obj, X)
            % FIRST Compute the first term.
            
            if obj.IsNormalized
                Y = repmat(1/sqrt(2), 1, size(X, 2));
            else
                Y = ones(1, size(X, 2));
            end
        end
        
        function [a, b, c] = computeCoefficients(obj, n)
            % COMPUTECOEFFICIENTS Compute Legendre recurrence coefficients.
            
            if obj.IsNormalized
                a = sqrt((2 * n + 3) * (2 * n + 1)) / (n + 1);
                b = 0;
                if n >= 1 
                    c = -n / (n + 1) * sqrt((2 * n + 3) / (2 * n - 1));
                else
                    c = 0;
                end
            else
                a = (2 * n + 1) / (n + 1);
                b = 0;
                c = -n / (n + 1);
            end
        end
    end
end