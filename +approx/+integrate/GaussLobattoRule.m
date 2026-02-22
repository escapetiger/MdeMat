classdef GaussLobattoRule < approx.integrate.OrthotopeRule
    % GAUSSLOBATTORULE Gauss-Lobatto numerical integration rule.
    %
    %   GaussLobattoRule implements numerical integration using the
    %   Gauss-Lobatto method, which includes both interval endpoints.
    %   This method provides high accuracy while constraining the nodes
    %   to include the boundary points.
    %
    %   The Gauss-Lobatto quadrature formula includes both endpoints of the
    %   integration interval. For \f$n\f$ nodes, the method has polynomial
    %   exactness of degree \f$2n-3\f$.
    %
    % See also:
    %   approx.integrate.GaussLegendreRule,
    %   approx.integrate.ClosedNewtonCotesRule,
    %   approx.integrate.MonteCarloOrthotopeRule

    methods (Access = protected)
        function [x, w] = generate1d(~, np, options)
            % GENERATE1D Generate 1D Gauss-Lobatto nodes and weights.
            %
            %   [x, w] = generate1d(obj, np) generates @a np Gauss-Lobatto
            %   nodes @a x and corresponding weights @a w for the interval
            %   \f$[0, 1]\f$.
            %
            %   [x, w] = generate1d(obj, np, lower=a, upper=b) generates
            %   nodes and weights for the interval \f$[a, b]\f$.

            arguments
                ~
                np(1, 1) {mustBePositive, mustBeInteger}
                options.lower(1, 1) {mustBeReal} = 0
                options.upper(1, 1) {mustBeReal} = 1
            end

            a = options.lower;
            b = options.upper;

            core.except.assert(a <= b, 'InvalidInput', ...
                'a <= b is required.');

            if np == 1
                x = (a + b) / 2;
                w = 1;
                return;
            end

            if a == b
                [x, w] = deal(a, 1);
                return;
            end

            x = cos(pi * ((np - 1):-1:0) / (np - 1));
            P = zeros(np, np);
            x0 = 2 * ones(size(x));

            while max(abs(x-x0)) > eps(1)
                x0 = x;
                P(1, :) = 1;
                P(2, :) = x;
                for k = 3:np
                    P(k, :) = ((2 * k - 3) * x .* P(k-1, :) - (k - 2) * P(k-2, :)) / (k - 1);
                end
                x = x0 - (x .* P(np, :) - P(np-1, :)) ./ (np * P(np, :));
            end
            x = (x + 1) / 2;
            w = 1 ./ ((np - 1) * np * P(np, :).^2);
            x = a + (b - a) * x;
            w = w * (b - a);
        end
    end
end