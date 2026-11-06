classdef GaussTrapezoidalRule < mdemat.common.quadrature.SphereRule
    % GAUSSTRAPEZOIDALRULE Gauss-Trapezoidal integration over
    % sphere.
    %
    %   GaussTrapezoidalRule implements numerical integration over
    %   spheres using a product of Gauss-Lobatto rules for zenith
    %   angles and a periodic trapezoidal rule for the azimuthal angle.
    %   This method provides efficient integration for functions on
    %   spherical domains.
    %   Higher-dimensional zenith factors use the Jacobi weight
    %   (1-mu_i^2)^((d-i-2)/2). In 2D/3D the existing periodic/Lobatto
    %   rules are retained. cossph uses mu_i=-cos(theta_i).
    %
    %   Spherical coordinates map Cartesian coordinates
    %   \f$(y_1,...,y_d)\f$ as:
    %   \f[
    %     y_1 = r*\cos(x_1)
    %     y_i = r*\sin(x_1)*...*\sin(x_{i-1})*\cos(x_i), i = 2,...,d-1
    %     y_d = r*\sin(x_1)*...*\sin(x_{d-1})
    %   \f]
    %   where:
    %   \f[
    %     x_i \in [0, \pi] for i = 1,...,d-2 (zenith angles)
    %     x_{d-1} \in [0, 2\pi] (azimuthal angle)
    %   \f]
    %
    % See also:
    %   mdemat.common.quadrature.IntegrationRule,
    %   mdemat.common.quadrature.LebedevRule,
    %   mdemat.common.quadrature.MonteCarloSphereRule

    methods (Access = protected)
        function [X, w] = generateImpl(obj, geometry, np)
            % GENERATEIMPL Generate Gauss-Trapezoidal quadrature rule.
            %
            %   [X, w] = generateImpl(obj, geometry, np, options) generates
            %   integration points @a X and corresponding weights @a w for
            %   the sphere @a geometry using @a np points in each
            %   spherical coordinate direction.

            arguments
                obj mdemat.common.quadrature.GaussTrapezoidalRule
                geometry mdemat.data.geometry.Sphere
                np(1, :) {mustBePositive, mustBeInteger}
            end

            coord = obj.Coord;
            r = geometry.Radius;
            nd = geometry.NDims;

            if nd == 1
                assert(all(numel(np) == 1 && np == 2, 'all'), ...
                    'mdemat:common:quadrature:GaussTrapezoidalRule:InvalidInput', '%s', 'For 0-sphere, np must be 2.');
                X = geometry.Center(:) + [-r, r];
                w = [1, 1];
                return;
            end

            k = nd - 1;
            assert(all(numel(np) == k, 'all'), 'mdemat:common:quadrature:GaussTrapezoidalRule:InvalidInput', ...
                'np must have %d elements for (%d)-sphere.', k, nd-1);

            X = cell(1, k);
            w = cell(1, k);

            G1 = mdemat.data.geometry.Orthotope([-1, 1]);
            GL = mdemat.common.quadrature.GaussLobattoRule(1);
            for i = 1:k - 1
                alpha = (nd-i-2)/2;
                if alpha == 0
                    [X{i}, w{i}] = GL.generate(G1, np(i));
                else
                    [X{i}, w{i}] = symmetricJacobi(np(i), alpha);
                end
            end

            G2 = mdemat.data.geometry.Orthotope([0, 2*pi]);
            PT = mdemat.common.quadrature.PeriodicTrapezoidalRule(1);
            [X{k}, w{k}] = PT.generate(G2, np(k));
            X{k} = mod(X{k} + pi/6, 2*pi);

            [X{1:k}] = ndgrid(X{:});
            [w{1:k}] = ndgrid(w{:});
            X = reshape(cat(k+1, X{:}), [], k).';
            w = r^k * prod(reshape(cat(k+1, w{:}), [], k).', 1);
            w = w(:).';

            if ismember(coord, ["stdsph", "car"])
                X(1:nd-2, :) = acos(-X(1:nd-2, :));
            end

            if strcmpi(coord, "car")
                X = geometry.sphericalToCartesian(X);
            end
        end
    end
end

function [nodes, weights] = symmetricJacobi(nPoints, alpha)
% Gaussian rule for (1-mu^2)^alpha dmu on [-1,1], alpha > 0.
mass = sqrt(pi)*exp(gammaln(alpha+1)-gammaln(alpha+1.5));
if nPoints == 1
    nodes = 0;
    weights = mass;
    return;
end
k = 1:nPoints-1;
beta = sqrt(k.*(k+2*alpha)./((2*k+2*alpha).^2-1));
[vectors, values] = eig(diag(beta, 1)+diag(beta, -1));
[nodes, order] = sort(diag(values));
nodes = nodes.';
weights = mass*vectors(1, order).^2;
end
