classdef GaussLegendreRule < approx.integrate.OrthotopeRule
    % GAUSSLEGENDRERULE Gauss-Legendre numerical integration rule.
    %
    %   GaussLegendreRule implements numerical integration using the
    %   Gauss-Legendre method, which provides optimal accuracy for
    %   polynomial integrands. For \f$n\f$ nodes, this method exactly
    %   integrates polynomials of degree up to \f$2n-1\f$.
    %
    %   The Gauss-Legendre quadrature formula achieves maximum polynomial
    %   exactness for a given number of nodes. The nodes are the roots of
    %   the Legendre polynomial of degree n, providing optimal quadrature
    %   points that do not include the interval endpoints.
    %
    % See also:
    %   approx.integrate.ClosedNewtonCotesRule,
    %   approx.integrate.GaussLobattoRule,
    %   approx.integrate.MonteCarloOrthotopeRule

    methods (Access = protected)
        function [x, w] = generate1d(~, np, options)
            % GENERATE1D Generate 1D Gauss-Legendre nodes and weights.
            %
            %   [x, w] = generate1d(obj, np) generates @a np Gauss-Legendre
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

            core.except.assert(a <= b, 'InvalidInput', 'a <= b is required.');

            if a == b
                [x, w] = deal(a, 1);
                return;
            end

            beta = 0.5 * (1 - (2 * (1:np - 1)).^(-2)).^(-0.5);
            T = diag(beta, 1) + diag(beta, -1);
            [V, D] = eig(T);
            x = (diag(D).' + 1) / 2;
            w = V(1, :).^2;
            x = a + (b - a) * x;
            w = w * (b - a);
        end
    end
end