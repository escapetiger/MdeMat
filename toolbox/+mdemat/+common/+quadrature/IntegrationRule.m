classdef IntegrationRule < handle
    % INTEGRATIONRULE Base class for all integration rules.
    %
    %   IntegrationRule provides a common interface for all numerical
    %   integration methods. This abstract class defines the basic
    %   structure and required methods that must be implemented by concrete
    %   integration rule subclasses.
    %
    %   The class defines two fundamental operations: interior domain
    %   integration and boundary integration. Concrete subclasses must
    %   implement the generate method for domain integration, while
    %   generateOnBoundary is optional depending on the specific rule.
    %
    % See also:
    %   mdemat.common.quadrature.OrthotopeRule, mdemat.common.quadrature.SphereRule

    properties
        NDims % Number of dimensions supported by this rule
    end

    methods
        function obj = IntegrationRule(nDims)
            % INTEGRATIONRULE Construct an instance of IntegrationRule.
            %
            %   obj = IntegrationRule(nDims) creates an integration rule
            %   instance that supports up to @a nDims dimensions.

            arguments
                nDims {mustBePositive, mustBeInteger}
            end

            obj.NDims = nDims;
        end

        function [X, w] = generate(obj, geometry, np)
            % GENERATE Generate integration points and weights.
            %
            %   [X, w] = generate(obj, geometry, np) generates integration
            %   points @a X and weights @a w within the geometric domain @a
            %   geometry using @a np points.

            arguments
                obj mdemat.common.quadrature.IntegrationRule
                geometry mdemat.data.geometry.Geometry
                np {mustBePositive, mustBeInteger}
            end

            assert(all(obj.NDims == geometry.NDims, 'all'), ...
                'mdemat:common:quadrature:IntegrationRule:InvalidInput', ['Integration rule supports %d', ...
                ' dimensions, but geometry has %d dimensions.'], ...
                obj.NDims, geometry.NDims);

            [X, w] = obj.generateImpl(geometry, np);
        end
    end

    methods (Abstract, Access = protected)
        % GENERATEIMPL Implementation of integration point generation.
        [X, w] = generateImpl(obj, geometry, np)
    end
end