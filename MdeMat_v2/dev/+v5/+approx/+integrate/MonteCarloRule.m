classdef MonteCarloRule < approx.integrate.IntegrationRule
    % MONTECARLORULE Monte Carlo numerical integration rule.
    %
    %   MonteCarloRule provides a foundation for Monte Carlo integration
    %   methods. This class manages random number generation with
    %   reproducible seeding and defines the interface for Monte Carlo
    %   quadrature implementations.
    %
    % Examples:
    %   % Create a 2D Monte Carlo rule with specific seed
    %   rule = MonteCarloRule(2, 12345);
    %   
    %   % Generate 1000 random integration points
    %   [X, w] = rule.generate(1000);
    %
    % Notes:
    %   This is an abstract class that must be subclassed to implement
    %   specific Monte Carlo integration methods. The generateImpl method
    %   must be implemented by concrete subclasses.
    %
    % See also:
    %   approx.integrate.IntegrationRule, approx.integrate.SphericalMonteCarloRule

    properties (Access = public)
        seed % Random number generator seed for reproducibility
    end
    
    methods
        function obj = MonteCarloRule(nDims, seed)
            % MONTECARLORULE Constructor for MonteCarloRule.
            %
            %   obj = MonteCarloRule(nDims) creates a Monte Carlo rule
            %   for the specified number of dimensions with default seed.
            %
            %   obj = MonteCarloRule(nDims, seed) creates a Monte Carlo
            %   rule with the specified random number generator seed.
            %
            % Inputs:
            %   nDims - Number of dimensions (positive integer)
            %   seed - Random number generator seed (optional, default: current rng state)
            %
            % Outputs:
            %   obj - Constructed MonteCarloRule object

            if nargin < 2, seed = rng; end

            obj@approx.integrate.IntegrationRule(nDims);
            obj.seed = seed;
        end

        function [X, w] = generate(obj, n, varargin)
            % GENERATE Generate Monte Carlo integration nodes and weights.
            %
            %   [X, w] = generate(obj, n) generates n Monte Carlo
            %   integration points and corresponding weights.
            %
            %   [X, w] = generate(obj, n, varargin) passes additional
            %   arguments to the concrete implementation.
            %
            % Inputs:
            %   obj - The MonteCarloRule object
            %   n - Number of Monte Carlo points (positive integer)
            %   varargin - Additional arguments for specific implementations
            %
            % Outputs:
            %   X - d×n matrix of integration nodes
            %   w - 1×n vector of integration weights
            
            core.except.assert(msclab.utilities.isPositiveInteger(n), ...
                'InvalidInput', 'n must be a positive integer.');
            
            oldSeed = rng;
            rng(obj.seed);

            [X, w] = obj.generateImpl(n, varargin{:});

            obj.seed = rng;
            rng(oldSeed);
        end
    end

    methods (Abstract, Access = protected)
        % GENERATEIMPL Implementation of Monte Carlo point generation.
        [X, w] = generateImpl(obj, n, varargin)
    end
end