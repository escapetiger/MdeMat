classdef HypersphericalGeometry < core.geometry.Geometry
    % HYPERSPHERICALGEOMETRY Base class for all hyperspherical geometries.
    %
    %   HypersphericalGeometry provides common functionality for
    %   hyperspherical geometries such as hyperball (solid sphere) and
    %   hypersphere (surface). It implements coordinate transformations
    %   between Cartesian and hyperspherical coordinates, and provides
    %   helper methods for computing distances between points in the
    %   domain.
    %
    % See also:
    %   core.geometry.Hyperball, core.geometry.Hypersphere
    
    properties
        centroid % Centeroid
        radius % Radius
    end
    
    methods
        function obj = HypersphericalGeometry(centroid, radius)
            % HYPERSPHERICALGEOMETRY Constructor for
            % HypersphericalGeometry.
            %
            %   obj = HypersphericalGeometry(centroid, radius) creates a
            %   hyperspherical geometry with the specified center and
            %   radius.
            %
            % Inputs:
            %   centroid - Centroid (vector)
            %   radius - Radius (positive scalar)
            %
            % Outputs:
            %   obj - Constructed HypersphericalGeometry object
            
            core.except.assert(nargin >= 2, ...
                'InvalidInput', 'Centroid and radius must be specified.');
            
            core.except.assert(isvector(centroid), ...
                'InvalidInput', 'Centroid must be a vector.');

            core.except.assert(isscalar(radius) && radius >= 0, ...
                'InvalidInput', 'Radius must be a nonnegative scalar.');
            
            obj@core.geometry.Geometry(length(centroid));
            obj.centroid = centroid(:).';
            obj.radius = radius;
        end
        
        function S = cartesianToSpherical(obj, C)
            % CARTESIANTOSPHERICAL Converts Cartesian to hyperspherical
            % coordinates.
            %
            %   S = cartesianToSpherical(obj, C) converts points from
            %   Cartesian coordinates to hyperspherical coordinates. The
            %   conversion uses the standard mathematical convention where
            %   the first coordinate is the radius, followed by (n-1)
            %   angles. This transformation is crucial for specialized
            %   algorithms that operate in spherical coordinates.
            %
            % Inputs:
            %   obj - The HypersphericalGeometry object
            %   C - Cartesian coordinates (dxm matrix) 
            %
            % Outputs:
            %   S - Hyperspherical coordinates (d×m matrix) 
            
            d = obj.nDims;
            if d == 1
                S = C;
                return;
            end

            s = size(C);
            core.except.assert(s(1) == d, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            C = reshape(C, s(1), []);
            
            core.except.assert( ...
                all(obj.isInside(C) | obj.isOnBoundary(C)), ...
                'InvalidInput', ...
                'Input points must be inside or on the boundary.');

            C = C - obj.centroid(:);
            S = zeros(size(C));
            r = sqrt(sum(C.^2, 1));
            S(1, :) = r;
            mr = (r ~= 0);
            if any(mr)
                for i = 2:(d - 1)
                    p = sqrt(sum(C((i-1):d, mr).^2, 1));
                    mp = (p ~= 0);
                    if any(mp)
                        ms = mr;
                        ms(mr) = mp;
                        S(i, ms) = acos(C(i-1, ms) ./ p(mp));
                    end
                end
                a = atan2(C(d, mr), C(d-1, mr));
                ma = (a < 0);
                a(ma) = a(ma) + 2 * pi;
                S(d, mr) = a;
            end
            S = reshape(S, s);
        end
        
        function C = sphericalToCartesian(obj, S)
            % SPHERICALTOCARTESIAN Converts hyperspherical to Cartesian
            % coordinates.
            %
            %   C = sphericalToCartesian(obj, S) converts points from
            %   hyperspherical coordinates back to Cartesian coordinates.
            %   This is the inverse of the cartesianToSpherical method. In
            %   hyperspherical coordinates, the first coordinate is the
            %   radius and the remaining coordinates are angles.
            %
            % Inputs:
            %   obj - The HypersphericalGeometry object
            %   S - Hyperspherical coordinates (matrix)
            %
            % Outputs:
            %   C - Cartesian coordinates (matrix)
            
            d = obj.nDims;
            if d == 1
                C = S;
                return;
            end

            s = size(S);
            
            core.except.assert(s(1) == d, ...
                'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');
            
            S = reshape(S, s(1), []);

            core.except.assert(all(S(1, :) >= 0), ...
                'InvalidInput', 'Radius must be non-negative.');

            core.except.assert( ...
                all(S(2:d-1, :) >= 0 & S(2:d-1, :) <= pi), ...
                'InvalidInput', ...
                'First d-2 angles must be in the range [0, pi].');
            
            core.except.assert( ...
                all(S(d, :) >= 0 & S(d, :) <= 2*pi), ...
                'InvalidInput', ...
                'Last angles must be in the range [0, 2*pi].');

            r = S(1, :);
            A = S(2:end, :);
            C = zeros(size(S));
            for i = 1:(d-1)
                u = ones(1, size(S, 2));
                for j = 1:(i-1)
                    u = u .* sin(A(j, :));
                end
                C(i, :) = r .* u .* cos(A(i, :));
            end
            u = ones(1, size(S, 2));
            for j = 1:(d-1)
                u = u .* sin(A(j, :));
            end
            C(d, :) = r .* u;
            C = C + obj.centroid(:);
            C = reshape(C, s);
        end
    end
end