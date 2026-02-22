classdef MonteCarloOrthotopeRule < approx.integrate.OrthotopeRule
    % MONTECARLOORTHOTOPERULE Monte Carlo numerical integration rule.
    %
    %   MonteCarloOrthotopeRule implements Monte Carlo integration methods using
    %   random sampling for numerical quadrature. This class manages
    %   random number generation with reproducible seeding and provides
    %   uniform random sampling over orthotope domains.
    %
    %   Monte Carlo methods approximate integrals by random sampling and
    %   averaging function values. The convergence rate is typically
    %   O(1/sqrt(n)) regardless of dimension, making it effective for
    %   high-dimensional integration problems.
    %
    % See also:
    %   GaussLegendreRule, GaussLobattoRule, MonteCarloSphereRule

    properties
        Seed % Random number generator seed for reproducibility
    end

    methods
        function obj = MonteCarloOrthotopeRule(nDims, options)
            % MONTECARLOORTHOTOPERULE Constructor for MonteCarloOrthotopeRule.
            %
            %   obj = MonteCarloOrthotopeRule(nDims) creates a Monte Carlo rule that
            %   supports up to @a nDims dimensions with default random
            %   seed.
            %
            %   obj = MonteCarloOrthotopeRule(nDims, seed=seed) creates a Monte Carlo
            %   rule with the specified random number generator seed @a seed.

            arguments
                nDims {mustBePositive, mustBeInteger}
                options.seed = []
            end

            obj@approx.integrate.OrthotopeRule(nDims);
            obj.Seed = options.seed;
        end
    end

    methods (Access = protected)
        function [x, w] = generate1d(obj, np, options)
            % GENERATE1D Generate 1D Monte-Carlo nodes and weights.
            %
            %   [x, w] = generate1d(obj, np) generates @a np Monte-Carlo
            %   nodes @a x and corresponding weights @a w for the interval
            %   [0, 1].
            %
            %   [x, w] = generate1d(obj, np, lower=a, upper=b) generates
            %   Monte-Carlo nodes and weights for the interval [a, b].

            arguments
                obj approx.integrate.MonteCarloOrthotopeRule
                np(1, 1) {mustBePositive, mustBeInteger}
                options.lower(1, 1) {mustBeReal} = 0
                options.upper(1, 1) {mustBeReal} = 1
            end

            seed = obj.Seed;
            a = options.lower;
            b = options.upper;

            if ~isempty(seed)
                oldSeed = rng;
                rng(seed);
            end

            x = a + (b - a) * rand(1, np);
            w = ones(1, np) * (b - a) / np;

            if ~isempty(seed)
                rng(oldSeed);
            end
        end
    end
end