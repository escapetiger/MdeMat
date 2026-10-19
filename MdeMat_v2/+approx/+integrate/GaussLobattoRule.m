classdef GaussLobattoRule < approx.integrate.UnivariateRule
    % GAUSSLOBATTORULE Gauss-Lobatto numerical integration rule.
    %
    %   GaussLobattoRule implements numerical integration using the
    %   Gauss-Lobatto method, which includes both interval endpoints.
    %   This method provides high accuracy while constraining the nodes
    %   to include the boundary points.
    %
    % Examples:
    %   % Create a Gauss-Lobatto rule
    %   rule = GaussLobattoRule();
    %   
    %   % Generate 5 nodes and weights for interval [0, 1]
    %   [x, w] = rule.generate(5);
    %   
    %   % Generate nodes and weights for custom interval
    %   [x, w] = rule.generate(4, false, -1, 1);
    %
    % Notes:
    %   The Gauss-Lobatto quadrature formula includes both endpoints of
    %   the integration interval. For n nodes, the method has polynomial
    %   exactness of degree 2n-3. Requires n >= 2 for meaningful results.
    %
    % See also:
    %   approx.integrate.UnivariateRule, 
    %   approx.integrate.GaussLegendreRule,
    %   approx.integrate.ClosedNewtonCotesRule

    methods        
        function [x, w] = generate(~, n, a, b)
            % GENERATE Generate 1D Gauss-Lobatto nodes and weights.
            %
            %   [x, w] = generate(obj, n) generates n Gauss-Lobatto nodes
            %   and corresponding weights for the interval [0, 1].
            %
            %   [x, w] = generate(obj, n, a, b) generates nodes and
            %   weights for the interval [a, b].
            %
            % Inputs:
            %   n - Number of quadrature nodes (n >= 2)
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

            if n == 1
                x = (a + b) / 2;
                w = 1;
                return;
            end

            if a == b
                [x, w] = deal(a, 1);
                return;
            end

            x = cos(pi * ((n - 1):-1:0) / (n - 1));
            P = zeros(n, n);
            x0 = 2 * ones(size(x));

            while max(abs(x-x0)) > eps(1)
                x0 = x;
                P(1, :) = 1;
                P(2, :) = x;
                for k = 3:n
                    P(k, :) = ((2 * k - 3) * x .* P(k-1, :) - (k - 2) * P(k-2, :)) / (k - 1);
                end
                x = x0 - (x .* P(n, :) - P(n-1, :)) ./ (n * P(n, :));
            end
            x = (x + 1) / 2;
            w = 1 ./ ((n - 1) * n * P(n, :).^2);
            x = a + (b - a) * x;
            w = w * (b - a);
        end
    end
end