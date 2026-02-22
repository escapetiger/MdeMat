classdef Hypersphere < core.geometry.HypersphericalGeometry
    % HYPERSPHERE Geometry for (d-1)-sphere domains.
    %
    %   Hypersphere represents a generalized geometric surface that extends
    %   the concept of circles and spherical surfaces to arbitrary
    %   dimensions. It represents the set of points at a constant distance
    %   (radius) from a center point in d-dimensional Euclidean space.
    %   Dimensional representations:
    %   - 0-sphere: A pair of points on a line
    %   - 1-sphere: A circle in a plane
    %   - 2-sphere: A sphere surface in 3D space
    %   - n-sphere: A surface in d-dimensional space
    %
    % Examples:
    %   % Create a 1-sphere (circle) with radius 3 centered at [1,2]
    %   sphere = core.geometry.Hypersphere([1,2], 3);
    %   circumference = sphere.magnitude();  % Returns 2*pi*3
    %
    %   % Create a unit 2-sphere
    %   sphere2D = core.geometry.Hypersphere(2);
    %
    % See also: 
    %   core.geometry.HypersphericalGeometry, core.geometry.Hyperball

    properties (Dependent)
        nSphereDims % Number of spherical dimensions (nDims - 1)
    end
    
    methods
        function obj = Hypersphere(varargin)
            % HYPERSPHERE Constructor for Hypersphere.
            %
            %   obj = Hypersphere(nSphereDims) creates a unit hypersphere
            %   centered at the origin with unit radius. The dimension of
            %   the sphere surface is specified by nSphereDims, and the
            %   ambient space dimension is nSphereDims+1.
            %
            %   obj = Hypersphere(centroid, radius) creates a hypersphere
            %   with the specified center point and radius.
            %
            % Inputs:
            %   varargin - Input arguments
            %<   nSphereDims - Spherical dimension of the n-sphere
            %<   centroid - Centroid of the n-sphere (vector)
            %<   radius - Radius of the n-sphere (positive scalar)
            %
            % Outputs:
            %   obj - Constructed Hypersphere object

            switch nargin
                case 1
                    nSphereDims = varargin{1};
                    core.except.assert(nSphereDims >= 0, ...
                        'InvalidInput', ...
                        'Sphere dimension must be a nonnegative integer.');
                    nDims = nSphereDims + 1;
                    centroid = zeros(1, nDims);
                    radius = 1;
                case 2
                    [centroid, radius] = varargin{:};
            end
            
            obj@core.geometry.HypersphericalGeometry(centroid, radius);
        end
       
        function n = get.nSphereDims(obj)
            % GET.NSPHEREDIMS Returns the number of spherical dimensions.
            
            n = obj.nDims - 1;
        end

        function m = magnitude(obj)
            % MAGNITUDE Returns the surface area of the hypersphere.
            %
            %   m = magnitude(obj) computes the (n-1)-dimensional surface
            %   area of the hypersphere using the general formula for
            %   hyperspheres. Special cases are optimized for common
            %   dimensions (0-sphere, 1-sphere, 2-sphere) for computational
            %   efficiency and accuracy.
            %
            % Inputs:
            %   obj - The Hypersphere object
            %
            % Outputs:
            %   m - Surface area of the hypersphere

            n = obj.nSphereDims;
            
            r = obj.radius;
            switch n
                case 0
                    tol = sqrt(eps);
                    m = 2 * (abs(r) > tol);
                case 1
                    m = 2 * pi * r;
                case 2
                    m = 4 * pi * r^2;
                otherwise
                    m = 2 * pi^((n+1)/2) / gamma((n+1)/2) * r^n;
            end
        end
        
        function TF = isInside(obj, X)
            % ISINSIDE Tests if points are inside the hypersphere.
            %
            %   TF = isInside(obj, X) determines whether points lie inside
            %   the hypersphere. Since a hypersphere is a surface with no
            %   interior volume, this method always returns false for all
            %   points, regardless of their position.
            %
            % Inputs:
            %   obj - The Hypersphere object
            %   X - Coordinates (matrix)
            %
            % Outputs:
            %   TF - All false

            d = obj.nDims;
            n = size(X);
            core.except.assert(n(1) == d, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');
            
            TF = false(1, size(X, 2));
        end
        
        function TF = isOnBoundary(obj, X)
            % ISONBOUNDARY Tests if points are on the hypersphere surface.
            %
            %   TF = isOnBoundary(obj, X) determines whether points lie
            %   exactly on the hypersphere surface by checking if their
            %   distance from the center equals the radius within numerical
            %   tolerance. For a hypersphere, this is equivalent to testing
            %   if points lie on the sphere.
            %
            % Inputs:
            %   obj - The Hypersphere object
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
        
        function S = cartesianToSpherical(obj, C)
            % CARTESIANTOSPHERICAL Converts Cartesian to spherical
            % coordinates.
            %
            %   S = cartesianToSpherical(obj, C) transforms points from
            %   Cartesian to hyperspherical coordinates, returning only
            %   angular coordinates since all points on the sphere have the
            %   same radius. The input points must lie on the hypersphere
            %   surface.
            %
            % Inputs:
            %   obj - The Hypersphere object
            %   C - A matrix of size d×m of points on the sphere surface
            %
            % Outputs:
            %   S - A matrix of size (d-1)×m of angular coordinates

            d = obj.nDims;
            n = size(C);
            core.except.assert(n(1) == d, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            core.except.assert(obj.isOnBoundary(C), ...
                'InvalidInput', ...
                'Input points must be on the hypersphere surface.');

            S = cartesianToSpherical@core.geometry.HypersphericalGeometry(obj, C);
            S = reshape(S, d, []);
            S = S(2:d, :);
            S = reshape(S, [d-1, n(2:end)]);
        end
        
        function C = sphericalToCartesian(obj, S)
            % SPHERICALTOCARTESIAN Converts spherical to Cartesian
            % coordinates.
            %
            %   C = sphericalToCartesian(obj, S) transforms points from
            %   hyperspherical to Cartesian coordinates, receiving only
            %   angular coordinates and using the sphere's radius
            %   automatically. The resulting points will lie on the
            %   hypersphere surface.
            %
            % Inputs:
            %   obj - The Hypersphere object
            %   S - A matrix of size (d-1)×m of angular coordinates
            %
            % Outputs:
            %   C - A matrix of size d×m of points in Cartesian coordinates
            
            d = obj.nDims;
            n = size(S);
            core.except.assert(n(1) == d-1, 'DimensionMismatch', ...
                'Point dimension must match spherical dimension.');
            
            r = obj.radius;
            S = cat(1, r * ones([1, n(2:end)]), S);
            C = sphericalToCartesian@core.geometry.HypersphericalGeometry(obj, S);
        end
    end
end