classdef PeriodicTrapezoidalRule < approx.integrate.UnivariateRule
    % PERIODICTRAPEZOIDAL Periodic trapezoidal numerical integration rule.
    %
    %   PeriodicTrapezoidalRule implements numerical integration using the
    %   trapezoidal method on periodic domains. This method is particularly
    %   efficient for integrating periodic functions, providing exponential
    %   convergence for smooth periodic integrands.
    %
    % Examples:
    %   % Create a periodic trapezoidal rule
    %   rule = PeriodicTrapezoidalRule();
    %   
    %   % Generate 8 nodes and weights for interval [0, 1]
    %   [x, w] = rule.generate(8);
    %   
    %   % Generate nodes for interval [0, 2π]
    %   [x, w] = rule.generate(16, false, 0, 2*pi);
    %
    % Notes:
    %   The periodic trapezoidal rule assumes the integrand is periodic
    %   over the integration interval. It excludes the right endpoint
    %   to avoid double-counting due to periodicity.
    %
    % See also:
    %   approx.integrate.UnivariateRule, 
    %   approx.integrate.ClosedNewtonCotesRule,
    %   approx.integrate.GaussLegendreRule

    methods
        function [x, w] = generate(~, n, a, b)
            % GENERATE Generate 1D periodic trapezoidal nodes and weights.
            %
            %   [x, w] = generate(obj, n) generates n equally spaced nodes
            %   and uniform weights for the interval [0, 1], excluding
            %   the right endpoint.
            %
            %   [x, w] = generate(obj, n, a, b) generates nodes and
            %   weights for the interval [a, b).
            %
            % Inputs:
            %   n - Number of quadrature nodes (n >= 1)
            %   a - Lower integration bound (optional, default: 0)
            %   b - Upper integration bound (optional, default: 1)
            %
            % Outputs:
            %   x - Row vector of quadrature nodes within [a, b)
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

            x = linspace(0, 1-1/n, n);
            w = ones(1, n) / n;
            x = a + (b - a) * x;
            w = w * (b - a);
        end
    end
end