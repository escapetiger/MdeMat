classdef LebedevRule < approx.integrate.IntegrationRule
    % LEBEDEVRULE Lebedev quadrature rule for spherical integration.
    %
    %   LebedevRule provides highly symmetric integration points on the
    %   surface of a sphere in 3D space. This method is particularly
    %   efficient for integrating functions with spherical symmetry and
    %   offers superior accuracy compared to standard product rules.
    %
    % Examples:
    %   % Create a Lebedev rule for 3D sphere
    %   rule = LebedevRule(3);
    %   
    %   % Generate 110 integration points on unit sphere
    %   [X, w] = rule.generate(110);
    %   
    %   % Custom sphere with center and radius
    %   [X, w] = rule.generate(194, false, [1, 0, 0], 2);
    %
    % Notes:
    %   Only specific numbers of integration points are available, as
    %   defined in the availableParams property. The method generates
    %   points with exact spherical symmetry.
    %
    % See also:
    %   approx.integrate.IntegrationRule, approx.integrate.SphericalMonteCarloRule

    properties (Constant)
        availableParams = [6, 14, 26, 38, 50, 86, 110, 146, 170, 194, ...
            302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, ...
            2354, 2702, 3074, 3470, 3890, 4334, 4802, 5810]; % Available numbers of Lebedev integration points
    end

    methods
        function [X, w] = generate(obj, n, c, r, s)
            % GENERATE Generate Lebedev quadrature nodes and weights.
            %
            %   [X, w] = generate(obj, n) generates n Lebedev integration
            %   points and weights on the unit sphere centered at origin.
            %
            %   [X, w] = generate(obj, n, c, r, s) generates points for
            %   a sphere with specified center, radius, and coordinate
            %   system.
            %
            % Inputs:
            %   obj - The LebedevRule object
            %   n - Number of integration points
            %   c - Centroid of the sphere (optional, default: [0, 0, 0])
            %   r - Radius of the sphere (optional, default: 1)
            %   s - Coordinate system (optional, default: 1)
            %
            % Notes:
            %   s = 1: Cartesian coordinates
            %   s = 2: Spherical coordinates
            %   s = 3: Modified spherical coordinates
            %
            % Outputs:
            %   X - 3×n matrix of integration nodes
            %   w - 1×n vector of integration weights

            if nargin < 4, c = [0, 0, 0]; end
            if nargin < 5, r = 1; end
            if nargin < 6, s = 1; end
            
            core.except.assert(ismember(n, obj.availableParams), ...
                'InvalidInput', ...
                ['Invalid number of points. Must be one of: ', ...
                sprintf('%d ', obj.availableParams)]);

            path = sprintf('%s/+msclab/+integrate/lebedev', pwd);
            addpath(path);
            [x, y, z, w] = feval(sprintf('ld%04d', n));
            rmpath(path);
            X = [x(:), y(:), z(:)].' * r + c(:);
            w = w(:).' * 4 * pi * r^2;

            if s >= 2
                G = msclab.geometry.Hypersphere(c, r);
                X = G.cartesianToSpherical(X);
            end

            if s >= 3, X(1, :) = -cos(X(1, :)); end
        end
    end
end