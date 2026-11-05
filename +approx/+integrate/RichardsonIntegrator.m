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
        
        function newObj = copy(obj)
            % COPY Create a deep copy of the integrator.
            %
            %   newObj = copy(obj) creates a new RichardsonIntegrator
            %   object with the same rule, number of levels, nodes, and
            %   weights as the original object.
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

        function obj = setPoints(obj, n, a, b)
            % SETPOINTS Generate multi-level integration points.
            %
            %   obj = setPoints(obj, n) generates integration points at
            %   multiple refinement levels with n points per dimension at
            %   the coarsest level for the unit orthotope.
            %
            %   obj = setPoints(obj, n, a, b) generates points for a custom
            %   orthotope [a,b] with optional normalization.
            %
            % Inputs:
            %   obj - The RichardsonIntegrator object
            %   n - Number of points per dimension at coarsest level
            %   a - Lower bounds vector (optional, default: zeros)
            %   b - Upper bounds vector (optional, default: ones)
            %
            % Outputs:
            %   obj - The RichardsonIntegrator object

            d = obj.nDims;
            if nargin < 3
                a = zeros(1, d); 
                b = ones(1, d);
            end

            if isscalar(n), n = repmat(n, 1, d); end

            m = obj.nLevels;
            obj.nodes = cell(1, m);
            obj.weights = cell(1, m);
            for i = 1:m
                [obj.nodes{i}, obj.weights{i}] = obj.refine(n, a, b, i);
            end
        end
    end
    
    methods (Access = private)
        function [X, w] = refine(obj, n, a, b, i)
            % REFINE Generate refined integration nodes and weights.

            d = obj.nDims;
            if i == 1
                args = arrayfun(@(k) {n(k), a(k), b(k)}, 1:d, 'Un', 0);
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
                    [Y{j}, v{j}] = obj.refine(n, p, q, i - 1);
                end
                X = [Y{:}];
                w = [v{:}];
            end
        end
    end
end