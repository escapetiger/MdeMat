classdef HypersphereIntegrator < approx.integrate.Integrator
    % HYPERSPHEREINTEGRATOR Numerical integration over hypersphere domains.
    %
    %   HypersphereIntegrator provides specialized integration capabilities
    %   for functions defined on hyperspheres of arbitrary dimension. It
    %   supports multiple integration methods including Gauss-Trapezoidal,
    %   Lebedev (for 3D), and Monte Carlo approaches.
    %
    % Examples:
    %   % Create integrator with Gauss-Trapezoidal rule for 3D sphere
    %   integrator = HypersphereIntegrator(3, 'gauss_trapezoidal');
    %   
    %   % Create integrator with Lebedev rule for 3D sphere
    %   integrator = HypersphereIntegrator(3, 'lebedev');
    %   
    %   % Create integrator with Monte Carlo rule
    %   integrator = HypersphereIntegrator(4, 'monte_carlo', 12345);
    %
    % Notes:
    %   Lebedev rules are only available for 3D spheres and provide highly
    %   symmetric integration points. Monte Carlo methods are suitable for
    %   high-dimensional problems and discontinuous integrands.
    %
    % See also:
    %   approx.integrate.Integrator, 
    %   approx.integrate.SphericalGaussTrapezoidalRule,
    %   approx.integrate.LebedevRule, 
    %   approx.integrate.SphericalMonteCarloRule

    methods
        function obj = HypersphereIntegrator(varargin)
            % HYPERSPHEREINTEGRATOR Constructor for HypersphereIntegrator.
            %
            %   obj = HypersphereIntegrator(rule) creates an integrator
            %   using the specified integration rule object.
            %
            %   obj = HypersphereIntegrator(nDims, ruleId) creates an
            %   integrator for the specified dimension and rule type.
            %
            %   obj = HypersphereIntegrator(nDims, ruleId, seed) creates
            %   a Monte Carlo integrator with the specified random seed.
            %
            % Inputs:
            %   varargin - Input arguments
            %<   rule - IntegrationRule object
            %<   nDims - Number of dimensions (positive integer)
            %<   ruleId - Rule identifier string: {'gauss_trapezoidal', 'lebedev', 'monte_carlo'}
            %<   seed - Random seed for Monte Carlo (optional)
            %
            % Outputs:
            %   obj - Constructed HypersphereIntegrator object

            if nargin == 1
                rule = varargin{1};
            elseif nargin == 2
                [nDims, ruleId] = varargin{:};
                if strcmp(ruleId, 'gauss_trapezoidal')
                    rule = approx.integrate.SphericalGaussTrapezoidalRule(nDims);
                elseif nDims == 3 && strcmp(ruleId, 'lebedev')
                    rule = approx.integrate.LebedevRule(nDims);
                else
                    error('Invalid rule identifier or dimension for specified rule.');
                end
            elseif nargin == 3
                [nDims, ruleId, seed] = varargin{:};
                if strcmp(ruleId, 'monte_carlo')
                    rule = approx.integrate.SphericalMonteCarloRule(nDims, seed);
                else
                    error('Three-argument constructor only supports monte_carlo rule.');
                end
            else
                error('Invalid number of arguments.');
            end
            obj@approx.integrate.Integrator(rule);
        end

        function newObj = copy(obj)
            % COPY Create a deep copy of the integrator.
            %
            %   newObj = copy(obj) creates a new HypersphereIntegrator
            %   object with the same rule, nodes, and weights as the
            %   original object.
            %
            % Inputs:
            %   obj - The HypersphereIntegrator object to copy
            %
            % Outputs:
            %   newObj - Deep copy of the original integrator

            newObj = approx.integrate.HypersphereIntegrator(obj.rule);
            newObj.nodes = obj.nodes;
            newObj.weights = obj.weights;
        end
    end
end