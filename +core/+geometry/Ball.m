classdef Ball < core.geometry.Geometry
    % BALL Geometry for n-ball domains.
    %
    %   Ball represents a generalized geometric domain that extends
    %   the concept of filled circles and solid spheres to arbitrary
    %   dimensions. It represents the region enclosed by an
    %   \f$(n-1)\f$-dimensional sphere in n-dimensional Euclidean space.
    %   Dimensional representations:
    %   - 1-ball: A line segment
    %   - 2-ball: A disk (filled circle)
    %   - 3-ball: A solid sphere
    %   - n-ball: A solid region in n-dimensional space
    %
    % See also:
    %   core.geometry.Geometry, core.geometry.Sphere

    properties
        Center % Center
        Radius % Radius
    end

    properties (Dependent)
        BBox % Bounding box
        Boundary % Boundary of the ball
    end

    methods
        function obj = Ball(center, radius)
            % BALL Construct an instance of Ball.
            %
            %   obj = Ball(center, radius) creates a ball with
            %   specified @a center and @a radius.

            arguments
                center {mustBeVector, mustBeNumeric}
                radius (1,1) {mustBeNonnegative, mustBeNumeric}
            end

            obj@core.geometry.Geometry(length(center));
            obj.Center = center(:).';
            obj.Radius = radius;
        end

        function bbox = get.BBox(obj)
            % GET.BBOX Returns the bounding box.

            bbox = [obj.Center - obj.Radius, obj.Center + obj.Radius].';
        end

        function m = magnitude(obj)
            % MAGNITUDE Returns the volume of the ball.
            %
            %   m = magnitude(obj) computes the n-dimensional volume of the
            %   ball using the general formula for hyperballs.

            arguments
                obj core.geometry.Ball
            end

            nd = obj.NDims;
            r = obj.Radius;
            switch nd
                case 1
                    m = 2 * r;
                case 2
                    m = pi * r^2;
                case 3
                    m = 4 / 3 * pi * r^3;
                otherwise
                    m = pi^(nd / 2) / gamma(nd/2+1) * r^nd;
            end
        end

        function TF = isInside(obj, X)
            % ISINSIDE Tests if points are inside the ball.
            %
            %   TF = isInside(obj, X) determines whether points @a X lie
            %   strictly within the interior of the ball, excluding
            %   the boundary surface. @a X should be an nd-by-np matrix.

            arguments
                obj core.geometry.Ball
                X {mustBeNumeric}
            end

            nd = obj.NDims;
            nx = size(X);
            core.except.assert(nx(1) == nd, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            r = obj.Radius;
            c = obj.Center(:);
            TF = sum((X - c).^2, 1) < r^2;
        end

        function TF = isOnBoundary(obj, X)
            % ISONBOUNDARY Tests if points are on the ball boundary.
            %
            %   TF = isOnBoundary(obj, X) determines whether points @a X lie
            %   exactly on the surface of the ball. @a X should be an
            %   nd-by-np matrix.

            arguments
                obj core.geometry.Ball
                X {mustBeNumeric}
            end

            nd = obj.NDims;
            nx = size(X);
            core.except.assert(nx(1) == nd, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            r = obj.Radius;
            c = obj.Center(:);
            tol = sqrt(eps);
            TF = abs(sum((X - c).^2, 1)-r^2) <= tol;
        end

        function sphere = get.Boundary(obj)
            % GET.BOUNDARY Returns the surface of the ball.
            %
            %   sphere = get.Boundary(obj) returns the corresponding
            %   Sphere object that represents the boundary surface of
            %   this ball. The boundary sphere has the same center and
            %   radius as the ball.

            arguments
                obj core.geometry.Ball
            end

            sphere = core.geometry.Sphere(obj.Center, obj.Radius);
        end

        function S = cartesianToSpherical(obj, C)
            % CARTESIANTOSPHERICAL Converts Cartesian to spherical
            % coordinates.
            %
            %   S = cartesianToSpherical(obj, C) converts points from
            %   Cartesian coordinates to spherical coordinates using
            %   standard mathematical convention where the first coordinate
            %   is the radius, followed by (n-1) angles.

            arguments
                obj core.geometry.Ball
                C {mustBeNumeric}
            end

            nd = obj.NDims;
            if nd == 1
                S = C;
                return;
            end

            nx = size(C);
            core.except.assert(nx(1) == nd, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            C = reshape(C, nx(1), []);

            isInDomain = all(obj.isInside(C) | obj.isOnBoundary(C));
            core.except.assert(isInDomain, 'InvalidInput', ...
                'Input points must be inside or on the boundary.');

            C = C - obj.Center(:);
            S = zeros(size(C));
            r = sqrt(sum(C.^2, 1));
            S(1, :) = r;
            mr = (r ~= 0);
            if any(mr)
                for i = 2:(nd - 1)
                    p = sqrt(sum(C((i-1):nd, mr).^2, 1));
                    mp = (p ~= 0);
                    if any(mp)
                        ms = mr;
                        ms(mr) = mp;
                        S(i, ms) = acos(C(i-1, ms) ./ p(mp));
                    end
                end
                a = atan2(C(nd, mr), C(nd-1, mr));
                ma = (a < 0);
                a(ma) = a(ma) + 2 * pi;
                S(nd, mr) = a;
            end
            S = reshape(S, nx);
        end

        function C = sphericalToCartesian(obj, S)
            % SPHERICALTOCARTESIAN Converts spherical to Cartesian
            % coordinates.
            %
            %   C = sphericalToCartesian(obj, S) converts points from
            %   spherical coordinates back to Cartesian coordinates.
            %   This is the inverse of the cartesianToSpherical method.

            arguments
                obj core.geometry.Ball
                S {mustBeNumeric}
            end

            nd = obj.NDims;
            if nd == 1
                C = S;
                return;
            end

            nx = size(S);

            core.except.assert(nx(1) == nd, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            S = reshape(S, nx(1), []);

            core.except.assert(all(S(1, :) >= 0), 'InvalidInput', ...
                'Radius must be non-negative.');

            isAngleInRange = all(S(2:nd-1, :) >= 0 & S(2:nd-1, :) <= pi);
            core.except.assert(isAngleInRange, 'InvalidInput', ...
                'First nd-2 angles must be in the range [0, pi].');

            isAngleInRange = all(S(nd, :) >= 0 & S(nd, :) <= 2*pi);
            core.except.assert(isAngleInRange, 'InvalidInput', ...
                'Last angles must be in the range [0, 2*pi].');

            r = S(1, :);
            A = S(2:end, :);
            C = zeros(size(S));
            for i = 1:(nd-1)
                u = ones(1, size(S, 2));
                for j = 1:(i-1)
                    u = u .* sin(A(j, :));
                end
                C(i, :) = r .* u .* cos(A(i, :));
            end
            u = ones(1, size(S, 2));
            for j = 1:(nd-1)
                u = u .* sin(A(j, :));
            end
            C(nd, :) = r .* u;
            C = C + obj.Center(:);
            C = reshape(C, nx);
        end
    end

    methods (Static)
        function obj = unit(nDims)
            % UNIT Creates a unit ball centered at origin.
            %
            %   obj = unit(nDims) creates a unit ball centered at
            %   the origin with radius 1 and dimension @a nDims.

            core.except.assert(nDims >= 1, 'InvalidInput', ...
                'Dimension must be a positive integer.');
            center = zeros(1, nDims);
            radius = 1;
            obj = core.geometry.Ball(center, radius);
        end
    end
end