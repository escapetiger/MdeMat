classdef Integrator < handle
    % INTEGRATOR Base class for all integrators.
    %
    %   Integrator provides a common interface for numerical integration
    %   across different domains and methods. It manages integration rules,
    %   nodes, weights, and provides the fundamental integration
    %   functionality that can be specialized by concrete subclasses.
    %
    % See also:
    %   approx.integrate.OrthotopeIntegrator, 
    %   approx.integrate.HypersphereIntegrator,
    %   approx.integrate.IntegrationRule
    
    properties (Access = public)
        rule    % Integration rule object
        nodes   % Integration nodes
        weights % Integration weights
    end

    properties (Dependent)
        nDims   % Number of dimensions
        nPoints % Number of integration points
    end

    methods
        function obj = Integrator(rule)
            % INTEGRATOR Constructor for Integrator.
            %
            %   obj = Integrator(rule) creates a numerical integrator
            %   using the specified integration rule.
            %
            % Inputs:
            %   rule - IntegrationRule object
            %
            % Outputs:
            %   obj - Constructed Integrator object

            obj.rule = rule;
            obj.nodes = [];
            obj.weights = [];
        end

        function n = get.nDims(obj)
            % GET.NDIMS Get the number of dimensions.
            %
            %   n = get.nDims(obj) returns the number of dimensions
            %   of the integration domain as defined by the rule.
            %
            % Inputs:
            %   obj - The Integrator object
            %
            % Outputs:
            %   n - Number of dimensions (positive integer)

            n = obj.rule.nDims;
        end

        function n = get.nPoints(obj)
            % GET.NPOINTS Get the number of integration points.
            %
            %   n = get.nPoints(obj) returns the total number of
            %   integration points currently stored in the integrator.
            %
            % Inputs:
            %   obj - The Integrator object
            %
            % Outputs:
            %   n - Number of integration points (non-negative integer)

            n = length(obj.weights);
        end

        function newObj = copy(obj)
            % COPY Create a deep copy of the integrator.
            %
            %   newObj = copy(obj) creates a new Integrator object
            %   with the same rule, nodes, and weights as the original.
            %
            % Inputs:
            %   obj - The Integrator object to copy
            %
            % Outputs:
            %   newObj - Deep copy of the original integrator

            newObj = approx.integrate.Integrator(obj.rule);
            newObj.nodes = obj.nodes;
            newObj.weights = obj.weights;
        end
    
        function obj = setPoints(obj, varargin)
            % SETPOINTS Set integration points using the underlying rule.
            %
            %   obj = setPoints(obj, varargin) generates integration nodes
            %   and weights using the underlying integration rule and
            %   stores them in the integrator object.
            %
            % Inputs:
            %   obj - The Integrator object
            %   varargin - Arguments passed to the integration rule
            %
            % Outputs:
            %   obj - The Integrator object

            [obj.nodes, obj.weights] = obj.rule.generate(varargin{:});
        end

        function obj = addWeightFunction(obj, f)
            % ADDWEIGHTS Add weight function for integral.
            %
            %   obj = addWeights(obj, f) scales the integration weights
            %   according to the specified weight function.
            %
            % Inputs:
            %   obj - The Integrator object
            %   f - Weight function handle
            %
            % Outputs:
            %   obj - The Integrator object

            omega = f(obj.nodes);
            obj.weights = obj.weights .* omega;
        end

        function I = integrate(obj, f)
            % INTEGRATE Evaluate the integral of function f.
            %
            %   I = integrate(obj, f) approximates the integral of
            %   function f over the integration domain using the
            %   precomputed nodes and weights.
            %
            % Inputs:
            %   obj - The Integrator object
            %   f - Function handle or cell array of functions to integrate
            %
            % Outputs:
            %   I - Approximate value of the integral
            
            X = obj.nodes;
            w = obj.weights;
            core.except.assert(~isempty(X) && ~isempty(w), ...
                'NotPrepared', 'Please set nodes and weights first.');

            if iscell(f)
                I = arrayfun(@(i) f{i}(X{i}) * w{i}(:), 1:numel(X), 'Un', 0);
            else
                I = f(X) * w(:);
            end
        end
    end
end