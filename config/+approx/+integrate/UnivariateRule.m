classdef UnivariateRule < approx.integrate.IntegrationRule
    % UNIVARIATERULE Base class for all univariate integration rules.
    %
    %   UnivariateRule provides a specialized interface for one-dimensional
    %   numerical integration methods. It extends the general
    %   IntegrationRule class with specific validation and structure
    %   appropriate for univariate quadrature formulas.
    %
    % See also:
    %   approx.integrate.IntegrationRule, 
    %   approx.integrate.GaussLegendreRule,
    %   approx.integrate.ClosedNewtonCotesRule, 
    %   approx.integrate.GaussLobattoRule
    
    methods
        function obj = UnivariateRule()
            % UNIVARIATERULE Constructor for UnivariateRule.
            %
            %   obj = UnivariateRule() creates a univariate integration
            %   rule object with dimension set to 1.
            %
            % Outputs:
            %   obj - Constructed UnivariateRule object

            obj@approx.integrate.IntegrationRule(1);
        end

        [X, w] = generate(obj, n, varargin)
    end
end