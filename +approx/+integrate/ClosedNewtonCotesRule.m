classdef ClosedNewtonCotesRule < approx.integrate.OrthotopeRule
    % CLOSEDNEWTONCOTESRULE Closed Newton-Cotes numerical integration rule.
    %
    %   ClosedNewtonCotesRule implements numerical integration using the
    %   closed Newton-Cotes method, which includes the interval endpoints.
    %   This method creates equally spaced nodes and calculates weights
    %   using Lagrange interpolation polynomials.
    %
    %   The closed Newton-Cotes formula includes both endpoints of the
    %   integration interval. For \f$n\f$ nodes, the method has polynomial
    %   exactness of degree \f$n-1\f$ for odd \f$n\f$ and \f$n\f$ for even
    %   \f$n\f$. The nodes are equally distributed across the integration
    %   domain.
    %
    % See also:
    %   approx.integrate.GaussLegendreRule, 
    %   approx.integrate.GaussLobattoRule, 
    %   approx.integrate.MonteCarloOrthotopeRule

    methods (Access = protected)
        function [x, w] = generate1d(~, np, options)
            % GENERATE1D Generate 1D closed Newton-Cotes nodes and weights.
            %
            %   [x, w] = generate1d(obj, np) generates @a np equally spaced
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

            if a == b
                [x, w] = deal(a, 1);
                return;
            end

            x = linspace(0, 1, np);
            w = zeros(1, np);
            for i = 1:np
                L = 1;
                for j = 1:np
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