classdef OrthotopeIntegrator < approx.integrate.Integrator
    % ORTHOTOPEINTEGRATOR Numerical integration over orthotope domains.
    %
    %   OrthotopeIntegrator provides specialized integration capabilities
    %   for functions defined on multidimensional orthotopes
    %   (hyperrectangles). It uses separable integration rules that
    %   construct multidimensional quadrature by tensor products of
    %   univariate rules.
    %
    % Examples:
    %   % Create integrator with Gauss-Legendre rule for 3D orthotope
    %   integrator = OrthotopeIntegrator(3, 'gauss_legendre');
    %   
    %   % Create integrator with Gauss-Lobatto rule for 2D rectangle
    %   integrator = OrthotopeIntegrator(2, 'gauss_lobatto');
    %   
    %   % Using existing rule object
    %   rule = approx.integrate.SeparableRule.gaussLegendre(4);
    %   integrator = OrthotopeIntegrator(rule);
    %
    % Notes:
    %   Orthotopes are multidimensional generalizations of rectangles.
    %   The separable rules are most efficient for functions that can
    %   be expressed as products of univariate functions.
    %
    % See also:
    %   approx.integrate.Integrator, approx.integrate.SeparableRule,
    %   approx.integrate.OrthotopeFaceIntegrator

    methods
        function obj = OrthotopeIntegrator(varargin)
            % ORTHOTOPEINTEGRATOR Constructor for OrthotopeIntegrator.
            %
            %   obj = OrthotopeIntegrator(rule) creates an integrator
            %   using the specified separable integration rule.
            %
            %   obj = OrthotopeIntegrator(nDims, ruleId) creates an
            %   integrator for the specified dimension and rule type.
            %
            % Inputs:
            %   varargin - Input arguments
            %<   rule - SeparableRule object
            %<   nDims - Number of dimensions (positive integer)
            %<   ruleId - Rule identifier string: {'gauss_legendre', 'gauss_lobatto'}
            %
            % Outputs:
            %   obj - Constructed OrthotopeIntegrator object

            if nargin == 1
                rule = varargin{1};
            elseif nargin == 2
                [nDims, ruleId] = varargin{:};
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
        end

        function newObj = copy(obj)
            % COPY Create a deep copy of the orthotope integrator.
            %
            %   newObj = copy(obj) creates a new OrthotopeIntegrator
            %   object with the same rule, nodes, and weights as the
            %   original object.
            %
            % Inputs:
            %   obj - The OrthotopeIntegrator object to copy
            %
            % Outputs:
            %   newObj - Deep copy of the original integrator

            newObj = approx.integrate.OrthotopeIntegrator(obj.rule);
            newObj.nodes = obj.nodes;
            newObj.weights = obj.weights;
        end

        function obj = setPoints(obj, n, a, b)
            % SETPOINTS Generate integration points for the orthotope.
            %
            %   obj = setPoints(obj, n) generates integration points with
            %   @a n points per dimension for the unit orthotope [0,1]^d.
            %
            %   obj = setPoints(obj, n, a, b) generates points for a custom
            %   orthotope [@a a, @a b].
            %
            % Inputs:
            %   obj - The OrthotopeIntegrator object
            %   n - Number of points per dimension (vector or scalar)
            %   a - Lower bounds vector (optional, default: zeros)
            %   b - Upper bounds vector (optional, default: ones)
            %
            % Outputs:
            %   obj - The OrthotopeIntegrator object

            d = obj.nDims;
            if nargin < 3
                a = zeros(1, d);
                b = ones(1, d);
            end

            args = arrayfun(@(k) {n(k), a(k), b(k)}, 1:d, 'Un', 0);
            [obj.nodes, obj.weights] = obj.rule.generate(args{:});
        end
    end
end