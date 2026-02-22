classdef SphericalGaussTrapezoidalRule < approx.integrate.IntegrationRule
    % SPHERICALGAUSSTRAPEZOIDALRULE Gauss-Trapezoidal integration over hypersphere.
    %
    %   SphericalGaussTrapezoidalRule implements numerical integration over
    %   hyperspheres using a product of Gauss-Lobatto rules for zenith angles
    %   and a periodic trapezoidal rule for the azimuthal angle. This method
    %   provides efficient integration for functions on spherical domains.
    %
    %   Hyperspherical coordinates map Cartesian coordinates 
    %   \f$(y_1,...,y_d)\f$ as:
    %   \f[  
    %     y_1 = r*\cos(x_1)
    %     y_i = r*\sin(x_1)*...*\sin(x_{i-1})*\cos(x_i), i = 2,...,d-1
    %     y_d = r*\sin(x_1)*...*\sin(x_{d-1})
    %   \f]
    %   where:
    %   \f[
    %     x_i \in [0, \pi] for i = 1,...,d-2 (zenith angles)
    %     x_{d-1} \in [0, 2\pi] (azimuthal angle)
    %   \f]
    %
    % Examples:
    %   % Create rule for 3D sphere integration
    %   rule = SphericalGaussTrapezoidalRule(3);
    %   
    %   % Generate integration points
    %   [X, w] = rule.generate([5, 8]);
    %   
    %   % Custom sphere with different center and radius
    %   [X, w] = rule.generate([3, 6], false, [1, 2, 0], 1.5);
    %
    % Notes:
    %   For d-dimensional hypersphere, requires d-1 parameters in n:
    %   first d-2 for zenith angles, last one for azimuthal angle.
    %
    % See also:
    %   approx.integrate.IntegrationRule, approx.integrate.LebedevRule,
    %   approx.integrate.SphericalMonteCarloRule

    methods        
        function [X, w] = generate(obj, n, c, r, s)
            % GENERATE Generate Gauss-Trapezoidal quadrature rule.
            %
            %   [X, w] = generate(obj, n) generates integration nodes and
            %   weights for a unit hypersphere centered at the origin.
            %
            %   [X, w] = generate(obj, n, c, r, s) generates nodes for a
            %   hypersphere with specified parameters and coordinate
            %   system.
            %
            % Inputs:
            %   obj - The SphericalGaussTrapezoidalRule object
            %   n - Vector of integration point numbers for each angle
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
            
            if d == 1
                core.except.assert(n == 2, 'InvalidInput', ...
                    'For 0-sphere, n must be 2.');
                X = [-r, r];
                w = [1, 1];
                return;
            end
            
            k = d - 1;
            X = cell(1, k);
            w = cell(1, k);
            GL = approx.integrate.GaussLobattoRule();
            for i = 1:k-1
                [X{i}, w{i}] = GL.generate(n(i), -1, 1);
            end
            
            PT = approx.integrate.PeriodicTrapezoidalRule();
            [X{k}, w{k}] = PT.generate(n(k), 0, 2*pi);

            [X{1:k}] = ndgrid(X{:});
            [w{1:k}] = ndgrid(w{:});
            X = reshape(cat(k+1, X{:}), [], k).';
            w = r^k * prod(reshape(cat(k+1, w{:}), [], k).', 1);
            w = w(:).';

            if s <= 2, X(1:d-2, :) = acos(-X(1:d-2, :)); end

            if s <= 1
                sphere = core.geometry.Hypersphere(c, r);
                X = sphere.sphericalToCartesian(X);
            end
        end
    end
end