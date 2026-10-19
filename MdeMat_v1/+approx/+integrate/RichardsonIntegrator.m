classdef RichardsonIntegrator < approx.integrate.Integrator
    % RICHARDSONINTEGRATOR Richardson extrapolation for numerical integration.
    %
    %   RichardsonIntegrator implements Richardson extrapolation to improve
    %   integration accuracy by combining quadrature rules at multiple
    %   resolution levels. This method leverages known convergence rates
    %   to achieve higher accuracy than the base integrator alone would
    %   provide at comparable computational cost.
    %
    % Examples:
    %   % Create Richardson integrator with 3 levels
    %   integrator = RichardsonIntegrator(3, 'gauss_legendre', 3);
    %   
    %   % Using existing rule object
    %   rule = approx.integrate.SeparableRule.gaussLobatto(2);
    %   integrator = RichardsonIntegrator(rule, 4);
    %
    % Notes:
    %   Richardson extrapolation assumes that the integration error can
    %   be expressed as a power series in the mesh size. The method
    %   recursively refines the domain and combines results to eliminate
    %   leading error terms.
    %
    % See also:
    %   approx.integrate.Integrator, approx.integrate.SeparableRule,
    %   approx.integrate.OrthotopeIntegrator

    properties (Access = public)
        nLevels % Number of refinement levels (positive integer)
    end

    methods
        function obj = RichardsonIntegrator(varargin)
            % RICHARDSONINTEGRATOR Constructor for RichardsonIntegrator.
            %
            %   obj = RichardsonIntegrator(rule, nLevels) creates a
            %   Richardson integrator using the specified rule and
            %   number of refinement levels.
            %
            %   obj = RichardsonIntegrator(nDims, ruleId, nLevels) creates
            %   a Richardson integrator with the specified parameters.
            %
            % Inputs:
            %   varargin - Input arguments
            %<   rule - SeparableRule object for the base integration
            %<   nLevels - Number of refinement levels (positive integer)
            %<   nDims - Number of dimensions (positive integer)
            %<   ruleId - Rule identifier string: {'gauss_legendre', 'gauss_lobatto'}
            %
            % Outputs:
            %   obj - Constructed RichardsonIntegrator object

            if nargin == 2
                [rule, nLevels] = varargin{:};
            elseif nargin == 3
                [nDims, ruleId, nLevels] = varargin{:};
                if strcmp(ruleId, 'gauss_legendre')
                    rule = approx.integrate.SeparableRule.gaussLegendre(nDims);
                elseif strcmp(ruleId, 'gauss_lobatto')
                    rule = approx.integrate.SeparableRule.gaussLobatto(nDims);
                else
                    error('Invalid rule identifier.');
                end
            else
                error('Invalid number of arguments.');
            end
            obj@approx.integrate.Integrator(rule);
            obj.nLevels = nLevels;
        end
        
        function varargout = setPoints(obj, n, f, a, b)
            % SETPOINTS Generate multi-level integration points.
            %
            %   setPoints(obj, n) generates integration points at multiple
            %   refinement levels with n points per dimension at the
            %   coarsest level for the unit orthotope.
            %
            %   setPoints(obj, n, f, a, b) generates points for a custom
            %   orthotope [a,b] with optional normalization.
            %
            % Inputs:
            %   obj - The RichardsonIntegrator object
            %   n - Number of points per dimension at coarsest level
            %   f - Normalization flag (optional, default: false)
            %   a - Lower bounds vector (optional, default: zeros)
            %   b - Upper bounds vector (optional, default: ones)
            %
            % Outputs:
            %   nodes - Cell array of integration nodes for each level
            %   weights - Cell array of integration weights for each level

            d = obj.nDims;
            if nargin < 3, f = false; end
            if nargin < 4, a = zeros(1, d); end
            if nargin < 5, b = ones(1, d); end

            if isscalar(n)
                n = repmat(n, 1, d);
            end

            m = obj.nLevels;
            X = cell(1, m);
            w = cell(1, m);
            for i = 1:m
                [X{i}, w{i}] = obj.refine(n, f, a, b, i);
            end
            if nargout == 0
                obj.nodes = X;
                obj.weights = w;
            else
                varargout{1} = X;
                varargout{2} = w;
            end
        end
    
        function newObj = copy(obj)
            % COPY Create a deep copy of the integrator.
            %
            %   newObj = copy(obj) creates a new RichardsonIntegrator
            %   object with the same rule, number of levels, nodes,
            %   and weights as the original object.
            %
            % Inputs:
            %   obj - The RichardsonIntegrator object to copy
            %
            % Outputs:
            %   newObj - Deep copy of the original integrator

            newObj = approx.integrate.RichardsonIntegrator(obj.rule, obj.nLevels);
            newObj.nodes = obj.nodes;
            newObj.weights = obj.weights;
        end
    end
    
    methods (Access = private)        
        function [X, w] = refine(obj, n, f, a, b, i)
            % REFINE Generate refined integration nodes and weights.

            d = obj.nDims;
            if i == 1
                args = arrayfun(@(k) {n(k), f, a(k), b(k)}, 1:d, 'Un', 0);
                [X, w] = obj.rule.generate(args{:});
            else
                m = 2^d;
                c = (a + b) / 2;
                Y = cell(1, m);
                v = cell(1, m);
                for j = 1:m
                    t = dec2bin(j - 1, d) - '0';
                    p = a .* (1 - t) + c .* t;
                    q = c .* (1 - t) + b .* t;
                    [Y{j}, v{j}] = obj.refine(n, f, p, q, i - 1);
                end
                X = [Y{:}];
                w = [v{:}];
            end
        end
    end
end