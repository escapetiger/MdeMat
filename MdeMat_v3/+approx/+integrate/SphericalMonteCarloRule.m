classdef SphericalMonteCarloRule < approx.integrate.MonteCarloRule
    % SPHERICALMONTECARLORULE Monte Carlo integration rule on hypersphere.
    %
    %   SphericalMonteCarloRule implements Monte Carlo integration over
    %   hyperspheres of arbitrary dimension using the normalized Gaussian
    %   method for uniform point generation. This method is particularly
    %   useful for high-dimensional integration problems and integrands
    %   with discontinuities or singularities.
    %
    % Examples:
    %   % Create rule for 3D sphere Monte Carlo integration
    %   rule = SphericalMonteCarloRule(3);
    %   
    %   % Generate 1000 random points on unit sphere
    %   [X, w] = rule.generate(1000);
    %   
    %   % Custom sphere with different parameters
    %   [X, w] = rule.generate(500, false, [1, 0, 0], 2, 2);
    %
    % Notes:
    %   Uses the normalized Gaussian method: generate multivariate normal
    %   random vectors and normalize to unit length. This ensures uniform
    %   distribution on the hypersphere surface.
    %
    % See also:
    %   approx.integrate.MonteCarloRule, approx.integrate.LebedevRule,
    %   approx.integrate.SphericalGaussTrapezoidalRule
    
    methods (Access = protected)
        function [X, w] = generateImpl(obj, n, c, r, s)
            % GENERATEIMPL Generate uniformly distributed nodes and weights.
            %
            %   [X, w] = generateImpl(obj, n) generates n uniformly
            %   distributed points on a unit hypersphere centered at origin.
            %
            %   [X, w] = generateImpl(obj, n, c, r, s) generates points for
            %   a hypersphere with specified parameters.
            %
            % Inputs:
            %   obj - The SphericalMonteCarloRule object
            %   n - Number of Monte Carlo points (positive integer)
            %   c - Centroid of the hypersphere (optional, default: zeros(1, d))
            %   r - Radius of the hypersphere (optional, default: 1)
            %   s - Coordinate system (optional, default: 1)
            %
            % Notes:
            %   s = 1: Cartesian coordinates
            %   s = 2: Spherical coordinates
            %   s = 3: Modified spherical coordinates
            %
            % Outputs:
            %   X - d×n matrix of integration nodes
            %   w - 1×n vector of integration weights

            
            d = obj.nDims;
            if nargin < 3, c = zeros(1, d); end
            if nargin < 4, r = 1; end
            if nargin < 5, s = 1; end

            X = mvnrnd(zeros(1, d), eye(d, d), n);
            X = (X ./ vecnorm(X, 2, 2)).';
            X = X * r + c(:);
            w = ones(1, n) / n;
            w = w * 2 * pi^(d/2) / gamma(d/2) * r^(d-1);

            if s >= 2
                G = msclab.geometry.Hypersphere(c, r);
                X = G.cartesianToSpherical(X);
            end

            if s >= 3, X(1:d-2, :) = -cos(X(1:d-2, :)); end
        end
    end
end