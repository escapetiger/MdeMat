classdef GaussTrapezoidalRule < approx.integrate.SphereRule
    % GAUSSTRAPEZOIDALRULE Gauss-Trapezoidal integration over
    % sphere.
    %
    %   GaussTrapezoidalRule implements numerical integration over
    %   spheres using a product of Gauss-Lobatto rules for zenith
    %   angles and a periodic trapezoidal rule for the azimuthal angle.
    %   This method provides efficient integration for functions on
    %   spherical domains.
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
    %   approx.integrate.IntegrationRule,
    %   approx.integrate.LebedevRule,
    %   approx.integrate.MonteCarloSphereRule

    methods (Access = protected)
        function [X, w] = generateImpl(obj, geometry, np)
            % GENERATEIMPL Generate Gauss-Trapezoidal quadrature rule.
            %
            %   [X, w] = generateImpl(obj, geometry, np, options) generates
            %   integration points @a X and corresponding weights @a w for
            %   the sphere @a geometry using @a np points in each
            %   spherical coordinate direction.

            arguments
                obj approx.integrate.GaussTrapezoidalRule
                geometry core.geometry.Sphere
                np(1, :) {mustBePositive, mustBeInteger}
            end

            coord = obj.Coord;
            r = geometry.Radius;
            nd = geometry.NDims;

            if nd == 1
                core.except.assert(numel(np) == 1 && np == 2, ...
                    'InvalidInput', 'For 0-sphere, np must be 2.');
                X = [-r, r];
                w = [1, 1];
                return;
            end

            k = nd - 1;
            core.except.assert(numel(np) == k, 'InvalidInput', ...
                'np must have %d elements for (%d)-sphere.', k, nd-1);

            X = cell(1, k);
            w = cell(1, k);

            G1 = core.geometry.Orthotope([-1, 1]);
            GL = approx.integrate.GaussLobattoRule(1);
            for i = 1:k - 1
                [X{i}, w{i}] = GL.generate(G1, np(i));
            end

            G2 = core.geometry.Orthotope([0, 2*pi]);
            PT = approx.integrate.PeriodicTrapezoidalRule(1);
            [X{k}, w{k}] = PT.generate(G2, np(k));

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