classdef IntegrationRule < handle
    % INTEGRATIONRULE Base class for all integration rules.
    %
    %   IntegrationRule provides a common interface for all numerical
    %   integration methods. This abstract class defines the basic
    %   structure and required methods that must be implemented by concrete
    %   integration rule subclasses.
    %
    % See also:
    %   approx.integrate.UnivariateRule, approx.integrate.GaussLegendreRule,
    %   approx.integrate.ClosedNewtonCotesRule

    properties (Access = public)
        nDims % Number of dimensions for the integration domain
    end

    methods
        function obj = IntegrationRule(nDims)
            % INTEGRATIONRULE Constructor for IntegrationRule.
            %
            %   obj = IntegrationRule(nDims) creates a new IntegrationRule
            %   object with the specified number of dimensions.
            %
            % Inputs:
            %   nDims - Number of dimensions (positive integer)
            %
            % Outputs:
            %   obj - Constructed IntegrationRule object

            obj.nDims = nDims;
        end
    end

    methods (Abstract)
        % GENERATE Generate integration nodes and weights.
        [X, w] = generate(obj, varargin)
    end
end