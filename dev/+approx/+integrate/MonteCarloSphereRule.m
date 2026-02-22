classdef MonteCarloSphereRule < approx.integrate.SphereRule
    % MONTECARLOSPHERERULE Monte Carlo integration rule on sphere.
    %
    %   MonteCarloSphereRule implements Monte Carlo integration over
    %   spheres of arbitrary dimension using the normalized Gaussian
    %   method for uniform point generation. This method is particularly
    %   useful for high-dimensional integration problems and integrands
    %   with discontinuities or singularities.
    %
    %   The implementation uses the normalized Gaussian method:
    %   multivariate normal random vectors are generated and normalized to
    %   unit length, ensuring uniform distribution on the sphere
    %   surface.
    %
    % See also:
    %   approx.integrate.MonteCarloOrthotopeRule,
    %   approx.integrate.LebedevRule,
    %   approx.integrate.GaussTrapezoidalRule

    properties
        Seed % Random number generator seed for reproducibility
    end

    methods
        function obj = MonteCarloSphereRule(nDims, options)
            % MONTECARLOSPHERERULE Constructor for
            % MonteCarloSphereRule.
            %
            %   obj = MonteCarloSphereRule(nDims) creates a Monte Carlo
            %   rule that supports up to @a nDims dimensions with default
            %   random seed.
            %
            %   obj = MonteCarloSphereRule(nDims, seed=seed) creates a
            %   Monte Carlo rule with the specified random number generator
            %   seed @a seed.

            arguments
                nDims {mustBePositive, mustBeInteger}
                options.seed = []
                options.coord {mustBeMember(options.coord, ["car", "stdsph", "cossph"])} = "car"
            end

            obj@approx.integrate.SphereRule(nDims, coord=options.coord);
            obj.Seed = options.seed;
        end
    end

    methods (Access = protected)
        function [X, w] = generateImpl(obj, geometry, np)
            % GENERATEIMPL Generate uniformly distributed nodes and
            % weights.
            %
            %   [X, w] = generateImpl(obj, geometry, np) generates
            %   @a np uniformly distributed points @a X and corresponding
            %   weights @a w on the sphere @a geometry.

            arguments
                obj approx.integrate.MonteCarloSphereRule
                geometry core.geometry.Sphere
                np(1, 1) {mustBePositive, mustBeInteger}
            end

            c = geometry.Center;
            r = geometry.Radius;
            nd = geometry.NDims;
            coord = obj.Coord;
            seed = obj.Seed;

            if ~isempty(seed)
                oldSeed = rng;
                rng(seed);
            end

            X = mvnrnd(zeros(1, nd), eye(nd, nd), np);
            X = (X ./ vecnorm(X, 2, 2)).';
            X = X * r + c(:);
            w = ones(1, np) / np;
            w = w * 2 * pi^(nd / 2) / gamma(nd/2) * r^(nd - 1);

            if ismember(coord, ["stdsph", "cossph"])
                X = geometry.cartesianToSpherical(X);
            end

            if strcmpi(coord, "cossph")
                X(1, :) = -cos(X(1, :));
            end

            if ~isempty(seed)
                rng(oldSeed);
            end
        end
    end
end