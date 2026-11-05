classdef GaussLegendreRule < approx.integrate.UnivariateRule
    % GAUSSLEGENDRERULE Gauss-Legendre numerical integration rule.
    %
    %   GaussLegendreRule implements numerical integration using the
    %   Gauss-Legendre method, which provides optimal accuracy for
    %   polynomial integrands. For n nodes, this method exactly integrates
    %   polynomials of degree up to 2n-1.
    %
    % Examples:
    %   % Create a Gauss-Legendre rule
    %   rule = GaussLegendreRule();
    %   
    %   % Generate 5 nodes and weights for interval [0, 1]
    %   [x, w] = rule.generate(5);
    %   
    %   % Generate nodes and weights for custom interval
    %   [x, w] = rule.generate(3, -1, 1);
    %
    % Notes:
    %   The Gauss-Legendre quadrature formula achieves maximum polynomial
    %   exactness for a given number of nodes. The nodes are the roots of
    %   the Legendre polynomial of degree n.
    %
    % See also:
    %   approx.integrate.UnivariateRule, 
    %   approx.integrate.ClosedNewtonCotesRule,
    %   approx.integrate.GaussLobattoRule

    methods  
        function [x, w] = generate(~, n, a, b)
            % GENERATE Generate 1D Gauss-Legendre nodes and weights.
            %
            %   [x, w] = generate(obj, n) generates n Gauss-Legendre nodes
            %   and corresponding weights for the interval [0, 1].
            %
            %   [x, w] = generate(obj, n, a, b) generates nodes and
            %   weights for the interval [a, b].
            %
            % Inputs:
            %   n - Number of quadrature nodes (n >= 1)
            %   a - Lower integration bound (optional, default: 0)
            %   b - Upper integration bound (optional, default: 1)
            %
            % Outputs:
            %   x - Row vector of quadrature nodes within [a, b]
            %   w - Row vector of quadrature weights

            if nargin < 3, [a, b] = deal(0, 1); end

            core.except.assert(isscalar(a) && isscalar(b), ...
                'InvalidInput', 'a and b must be scalars.');
            core.except.assert(a <= b, 'InvalidInput', ...
                'a <= b is required.');

            if a == b
                [x, w] = deal(a, 1);
                return;
            end

            beta = 0.5 * (1 - (2 * (1:n - 1)).^(-2)).^(-0.5);
            T = diag(beta, 1) + diag(beta, -1);
            [V, D] = eig(T);
            x = (diag(D).' + 1) / 2;
            w = V(1, :).^2;
            x = a + (b - a) * x;
            w = w * (b - a);
        end
    end
end