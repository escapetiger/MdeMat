classdef SeparableRule < core.linalg.Separable & approx.integrate.IntegrationRule
    % SEPARABLERULE Separable multidimensional integration rule.
    %
    %   SeparableRule creates multidimensional integration rules by taking
    %   tensor products of one-dimensional rules. This approach is
    %   efficient for integrating functions that can be separated into
    %   products of univariate functions or for rectangular integration
    %   domains.
    %
    % Examples:
    %   % Create a 2D separable rule from two Gauss-Legendre rules
    %   rule1D = approx.integrate.GaussLegendreRule();
    %   rule2D = SeparableRule({rule1D, rule1D});
    %   
    %   % Generate integration points
    %   [X, w] = rule2D.generate({5}, {5});
    %   
    %   % Use static methods for common rules
    %   gaussRule3D = SeparableRule.gaussLegendre(3);
    %
    % Notes:
    %   The separable rule constructs multidimensional quadrature by
    %   taking tensor products of univariate rules. This is most efficient
    %   for functions that exhibit separability or low-rank structure.
    %
    % See also:
    %   core.linalg.SeparableObject, approx.integrate.IntegrationRule,
    %   approx.integrate.GaussLegendreRule

    methods
        function obj = SeparableRule(rules)
            % SEPARABLERULE Constructor for SeparableRule.
            %
            %   obj = SeparableRule() creates an empty separable rule.
            %
            %   obj = SeparableRule(rules) creates a separable rule from
            %   a cell array of univariate integration rules.
            %
            % Inputs:
            %   rules - Cell array of IntegrationRule objects (optional)
            %
            % Outputs:
            %   obj - Constructed SeparableRule object

            if nargin > 0
                cls = {'approx.integrate.IntegrationRule'};
                core.except.assert(core.validate.isAllClass(rules, cls), ...
                    'InvalidFactor', ...
                    'Factors must be one of {%s}.', strjoin(cls, ', '));
            end

            nDims = 0;
            if nargin > 0 && ~isempty(rules)
                for i = 1:length(rules)
                    nDims = nDims + rules{i}.nDims;
                end
            end

            obj@core.linalg.Separable(rules);
            obj@approx.integrate.IntegrationRule(nDims);
        end

        function [X, w] = generate(obj, varargin)
            % GENERATE Generate separable integration rule nodes and weights.
            %
            %   [X, w] = generate(obj, varargin) generates multidimensional
            %   integration nodes and weights by taking tensor products of
            %   the constituent univariate rules.
            %
            % Inputs:
            %   obj - The SeparableRule object
            %   varargin - Cell array of arguments for each factor rule
            %
            % Outputs:
            %   X - d×m matrix of integration nodes
            %   w - 1×m vector of integration weights

            d = obj.nDims;
            F = obj.factors;
            X = cell(1, d);
            w = cell(1, d);

            for i = 1:d
                if iscell(F)
                    Fi = F{i};
                else 
                    Fi = F(i);
                end
                [X{i}, w{i}] = Fi.generate(varargin{i}{:});
            end

            [X{1:d}] = ndgrid(X{:});
            X = reshape(cat(d+1, X{:}), [], d).';

            [w{1:d}] = ndgrid(w{:});
            w = reshape(prod(cat(d+1, w{:}), d+1), 1, []);
        end
    end

    methods (Static)
        function obj = closedNewtonCotes(nDims)
            % CLOSEDNEWTONCOTES Create multidimensional closed Newton-Cotes
            % rule.
            %
            %   obj = closedNewtonCotes(nDims) creates a separable rule
            %   using closed Newton-Cotes quadrature in each dimension.
            %
            % Inputs:
            %   nDims - Number of dimensions (positive integer)
            %
            % Outputs:
            %   obj - SeparableRule instance with closed Newton-Cotes factors

            rules = cell(1, nDims);
            for i = 1:nDims
                rules{i} = approx.integrate.ClosedNewtonCotesRule();
            end

            obj = approx.integrate.SeparableRule(rules);
        end

        function obj = gaussLegendre(nDims)
            % GAUSSLEGENDRE Create multidimensional Gauss-Legendre rule.
            %
            %   obj = gaussLegendre(nDims) creates a separable rule using
            %   Gauss-Legendre quadrature in each dimension.
            %
            % Inputs:
            %   nDims - Number of dimensions (positive integer)
            %
            % Outputs:
            %   obj - SeparableRule instance with Gauss-Legendre factors

            rules = cell(1, nDims);
            for i = 1:nDims
                rules{i} = approx.integrate.GaussLegendreRule();
            end

            obj = approx.integrate.SeparableRule(rules);
        end

        function obj = gaussLobatto(nDims)
            % GAUSSLOBATTO Create multidimensional Gauss-Lobatto rule.
            %
            %   obj = gaussLobatto(nDims) creates a separable rule using
            %   Gauss-Lobatto quadrature in each dimension.
            %
            % Inputs:
            %   nDims - Number of dimensions (positive integer)
            %
            % Outputs:
            %   obj - SeparableRule instance with Gauss-Lobatto factors

            rules = cell(1, nDims);
            for i = 1:nDims
                rules{i} = approx.integrate.GaussLobattoRule();
            end

            obj = approx.integrate.SeparableRule(rules);
        end
    end
end