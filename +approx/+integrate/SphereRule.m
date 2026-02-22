classdef SphereRule < approx.integrate.IntegrationRule
    % SPHERERULE Integration rule over the sphere.
    %
    %   SphereRule creates multidimensional integration rules on the
    %   sphere. This approach is efficient for integrating functions
    %   that can be separated into products of univariate functions or for
    %   rectangular integration domains.
    %
    % See also:
    %   approx.integrate.IntegrationRule

    properties
        Coord % Coordinate system for output points
    end

    methods
        function obj = SphereRule(nDims, options)
            % SPHERERULE Construct an instance of SphereRule.
            %
            %   obj = SphereRule(nDims) creates an integration rule
            %   instance that supports up to @a nDims dimensions.

            arguments
                nDims {mustBePositive, mustBeInteger}
                options.coord {mustBeMember(options.coord, ["car", "stdsph", "cossph"])} = "car"
            end
            obj@approx.integrate.IntegrationRule(nDims);
            obj.Coord = options.coord;
        end
    end

    methods (Abstract, Access = protected)
        % GENERATEIMPL Implementation of integration point generation.
        [X, w] = generateImpl(obj, geometry, np)
    end
end