classdef Hyperball < core.geometry.HypersphericalGeometry
    % HYPERBALL Geometry for n-ball domains.
    %
    %   Hyperball represents a generalized geometric domain that extends
    %   the concept of filled circles and solid spheres to arbitrary
    %   dimensions. It represents the region enclosed by an
    %   \f$(n-1)\f$-dimensional sphere in n-dimensional Euclidean space.
    %   Dimensional representations:
    %   - 1-ball: A line segment
    %   - 2-ball: A disk (filled circle)
    %   - 3-ball: A solid sphere
    %   - n-ball: A solid region in n-dimensional space
    %
    % Examples:
    %   % Create a 2-ball (disk) with radius 3 centered at [1,2]
    %   ball = core.geometry.Hyperball([1,2], 3);
    %   area = ball.magnitude();  % Returns pi*3^2 (disk area)
    %
    %   % Create a unit 3-ball centered at origin
    %   ball3D = core.geometry.Hyperball(3);
    %
    % See also:
    %   core.geometry.HypersphericalGeometry, core.geometry.Hypersphere

    properties (Dependent)
        boundary % Boundary of the hyperball
    end

    methods
        function obj = Hyperball(varargin)
            % HYPERBALL Constructor for Hyperball.
            %
            %   obj = Hyperball(nDims) creates a unit hyperball
            %   (n-dimensional solid ball) centered at the origin with
            %   radius 1. The dimension of the ball is specified by nDims.
            %
            %   obj = Hyperball(centroid, radius) creates a hyperball with
            %   the specified center point and radius. The dimension is
            %   determined by the length of the centroid vector.
            %
            % Inputs:
            %   varargin - Input arguments
            %<   nDims - Dimension of the n-ball
            %<   centroid - Centeroid of the n-ball (vector)
            %<   radius - Radius of the n-ball (positive scalar)
            %
            % Outputs:
            %   obj - Constructed Hyperball object

            switch nargin
                case 1
                    nDims = varargin{1};
                    core.except.assert(nDims >= 1, ...
                        'InvalidInput', ...
                        'Ball dimension must be a positive integer.');
                    nDims = varargin{1};
                    centroid = zeros(1, nDims);
                    radius = 1;
                case 2
                    [centroid, radius] = varargin{:};
            end

            obj@core.geometry.HypersphericalGeometry(centroid, radius);
        end

        function m = magnitude(obj)
            % MAGNITUDE Returns the volume of the hyperball.
            %
            %   m = magnitude(obj) computes the n-dimensional volume of the
            %   hyperball using the general formula for hyperballs. Special
            %   cases are optimized for common dimensions (1D, 2D, 3D) for
            %   computational efficiency and accuracy.
            %
            % Inputs:
            %   obj - The Hyperball object
            %
            % Outputs:
            %   m - Volume of the hyperball

            d = obj.nDims;
            r = obj.radius;
            switch d
                case 1
                    m = 2 * r;
                case 2
                    m = pi * r^2;
                case 3
                    m = 4 / 3 * pi * r^3;
                otherwise
                    m = pi^(d / 2) / gamma(d/2+1) * r^d;
            end
        end

        function TF = isInside(obj, X)
            % ISINSIDE Tests if points are inside the hyperball.
            %
            %   TF = isInside(obj, X) determines whether points lie
            %   strictly within the interior of the hyperball, excluding
            %   the boundary surface. Points are inside if their distance
            %   from the center is strictly less than the radius.
            %
            % Inputs:
            %   obj - The Hyperball object
            %   X - Coordinates (matrix)
            %
            % Outputs:
            %   TF - True if the point is inside the hyperball

            d = obj.nDims;
            n = size(X);
            core.except.assert(n(1) == d, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            r = obj.radius;
            c = obj.centroid(:);
            TF = sum((X - c).^2, 1) < r^2;
        end

        function TF = isOnBoundary(obj, X)
            % ISONBOUNDARY Tests if points are on the hyperball boundary.
            %
            %   TF = isOnBoundary(obj, X) determines whether points lie
            %   exactly on the surface of the hyperball (the corresponding
            %   hypersphere). Points are on the boundary if their distance
            %   from the center equals the radius within numerical
            %   tolerance.
            %
            % Inputs:
            %   obj - The Hyperball object
            %   X - Coordinates (matrix)
            %
            % Outputs:
            %   TF - True if the point is on the boundary

            d = obj.nDims;
            n = size(X);
            core.except.assert(n(1) == d, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            r = obj.radius;
            c = obj.centroid(:);
            tol = sqrt(eps);
            TF = abs(sum((X - c).^2, 1)-r^2) <= tol;
        end

        function sphere = get.boundary(obj)
            % GET.BOUNDARY Returns the surface of the hyperball.
            %
            %   sphere = get.boundary(obj) returns the corresponding
            %   Hypersphere object that represents the boundary surface of
            %   this hyperball. The boundary sphere has the same center and
            %   radius as the hyperball.
            %
            % Inputs:
            %   NULL
            %
            % Outputs:
            %   sphere - Hypersphere object representing the boundary

            sphere = core.geometry.Hypersphere(obj.centroid, obj.radius);
        end
    end
end