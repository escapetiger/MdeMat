classdef ClosedNewtonCotesRule < approx.integrate.UnivariateRule
    % CLOSEDNEWTONCOTESRULE Closed Newton-Cotes numerical integration rule.
    %
    %   ClosedNewtonCotesRule implements numerical integration using the
    %   closed Newton-Cotes method, which includes the interval endpoints.
    %   This method creates equally spaced nodes and calculates weights
    %   using Lagrange interpolation polynomials.
    %
    % Examples:
    %   % Create a closed Newton-Cotes rule
    %   rule = ClosedNewtonCotesRule();
    %   
    %   % Generate 5 nodes and weights for interval [0, 1]
    %   [x, w] = rule.generate(5);
    %   
    %   % Generate nodes and weights for custom interval
    %   [x, w] = rule.generate(5, -1, 1);
    %
    % Notes:
    %   The closed Newton-Cotes formula includes both endpoints of the
    %   integration interval. For n nodes, the method has polynomial
    %   exactness of degree n-1 for odd n and n for even n.
    %
    % See also:
    %   approx.integrate.UnivariateRule, 
    %   approx.integrate.GaussLegendreRule,
    %   approx.integrate.GaussLobattoRule

    methods
        function [x, w] = generate(~, n, a, b)
            % GENERATE Generate 1D closed Newton-Cotes nodes and weights.
            %
            %   [x, w] = generate(obj, n) generates n equally spaced nodes
            %   and corresponding weights for the interval [0, 1].
            %
            %   [x, w] = generate(obj, n, a, b) generates nodes and weights
            %   for the interval [a, b].
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

            x = linspace(0, 1, n);
            w = zeros(1, n);
            for i = 1:n
                L = 1;
                for j = 1:n
                    if i ~= j
                        L = conv(L, [1, -x(j)]) / (x(i) - x(j));
                    end
                end
                I = polyint(L);
                w(i) = polyval(I, 1) - polyval(I, 0);
            end

            x = a + (b - a) * x;
            w = w * (b - a);
        end
    end
end