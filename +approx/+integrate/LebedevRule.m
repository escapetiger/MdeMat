classdef LebedevRule < approx.integrate.SphereRule
    % LEBEDEVRULE Lebedev quadrature rule for spherical integration.
    %
    %   LebedevRule provides highly symmetric integration points on the
    %   surface of a sphere in 3D space. This method is particularly
    %   efficient for integrating functions with spherical symmetry and
    %   offers superior accuracy compared to standard product rules.
    %
    %   The method generates points with exact spherical symmetry and
    %   supports only specific numbers of integration points that preserve
    %   the octahedral symmetry of the sphere.
    %
    % See also:
    %   approx.integrate.IntegrationRule, 
    %   approx.integrate.MonteCarloSphereRule

    properties (Constant)
        AvailableParams = [6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, ...
            230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, ...
            2354, 2702, 3074, 3470, 3890, 4334, 4802, 5810] % Available numbers of Lebedev integration points
    end

    methods
        function obj = LebedevRule(options)
            % LEBEDEVRULE Construct an instance of LebedevRule.
            %
            %   obj = LebedevRule() creates a Lebedev quadrature rule
            %   specifically for 3D spherical integration.

            arguments
                options.coord {mustBeMember(options.coord, ["car", "stdsph", "cossph"])} = "car"
            end

            obj@approx.integrate.SphereRule(3, coord = options.coord);
        end
    end

    methods (Access = protected)
        function [X, w] = generateImpl(obj, geometry, np)
            % GENERATEIMPL Implementation of Lebedev quadrature nodes and
            % weights generation.
            %
            %   [X, w] = generateImpl(obj, geometry, np) generates
            %   @a np Lebedev integration points @a X and corresponding
            %   weights @a w on the 3D sphere @a geometry.

            arguments
                obj approx.integrate.LebedevRule
                geometry core.geometry.Sphere
                np(1, 1) {mustBePositive, mustBeInteger}
            end
            
            core.except.assert(geometry.NDims == 3, 'InvalidInput', ...
                'Lebedev rule is only valid for 3D sphere.');

            core.except.assert(ismember(np, obj.AvailableParams), ...
                'InvalidInput', ...
                ['Invalid number of points. Must be one of: ', ...
                sprintf('%d ', obj.AvailableParams)]);

            lebedevPath = fullfile(fileparts(mfilename('fullpath')), 'lebedev');
            
            c = geometry.Center;
            r = geometry.Radius;

            addpath(lebedevPath);
            [x, y, z, w] = feval(sprintf('ld%04d', np));
            rmpath(lebedevPath);
            X = [x(:), y(:), z(:)].' * r + c(:);
            w = w(:).' * 4 * pi * r^2;

            if ismember(obj.Coord, ["stdsph", "cossph"])
                X = geometry.cartesianToSpherical(X);
            end

            if strcmpi(obj.Coord, "cossph")
                X(1, :) = -cos(X(1, :));
            end
        end
    end
end