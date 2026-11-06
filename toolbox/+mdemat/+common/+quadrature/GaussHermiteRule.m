classdef GaussHermiteRule < mdemat.common.quadrature.OrthotopeRule
    % GAUSSHERMITERULE Gauss-Hermite numerical integration rule.
    %
    %   GaussHermiteRule implements numerical integration using the
    %   Gauss-Hermite method for integrands weighted by the Gaussian
    %   function. This method is optimal for integrating functions over
    %   the entire real line with the weight function:
    %
    %   \f[
    %     w(x) = \frac{e^{-x^2/2}}{\sqrt{2\pi}}
    %   \f]
    %
    %   For \f$n\f$ nodes, this method exactly integrates polynomials
    %   weighted by the Gaussian function up to degree \f$2n-1\f$. The
    %   nodes are the roots of the physicist's Hermite polynomial
    %   \f$He_n(x)\f$ scaled for the weight \f$e^{-x^2/2}\f$.
    %
    %   This quadrature is particularly useful for probability integrals
    %   and Fourier transforms involving Gaussian functions.
    %
    % See also:
    %   mdemat.common.quadrature.GaussLegendreRule,
    %   mdemat.common.quadrature.GaussLobattoRule,
    %   mdemat.common.quadrature.IntegrationRule

    methods (Access = protected)
        function [x, w] = generate1d(~, np, options)
            % GENERATE1D Generate 1D Gauss-Hermite nodes and weights.
            %
            %   [x, w] = generate1d(obj, np) generates @a np Gauss-Hermite
            %   nodes @a x and corresponding weights @a w for the weight
            %   function e^{-x^2/2}/sqrt(2*pi) over the infinite domain.
            %
            %   [x, w] = generate1d(obj, np, lower=a, upper=b) generates
            %   nodes and weights, but ignores lower and upper bounds since
            %   Gauss-Hermite quadrature is defined over the infinite
            %   domain.

            arguments
                ~
                np(1, 1) {mustBePositive, mustBeInteger}
                options.lower(1, 1) {mustBeReal} = -inf
                options.upper(1, 1) {mustBeReal} = inf
            end

            if np == 1
                x = 0;
                w = 1;
                return;
            end

            i = 1:np - 1;
            beta = sqrt(i);
            J = diag(beta, 1) + diag(beta, -1);
            [V, D] = eig(J);
            x = diag(D);
            [~, idx] = sort(x);
            x = x(idx);
            V = V(:, idx);
            w = V(1, :).^2;
            % The Gaussian measure is even. Restore the exact +/- symmetry
            % that a floating-point eigensolve can perturb; odd moments
            % otherwise amplify unequal paired weights at high degree.
            x = (x-flipud(x))/2;
            w = (w+fliplr(w))/2;
            w = w/sum(w);
        end
    end
end
