classdef Sphere < core.geometry.Geometry
    % SPHERE Geometry for (d-1)-sphere domains.
    %
    %   Sphere represents a generalized geometric surface that extends
    %   the concept of circles and spherical surfaces to arbitrary
    %   dimensions. It represents the set of points at a constant distance
    %   (radius) from a center point in d-dimensional Euclidean space.
    %   Dimensional representations:
    %   - 0-sphere: A pair of points on a line
    %   - 1-sphere: A circle in a plane
    %   - 2-sphere: A sphere surface in 3D space
    %   - n-sphere: A surface in d-dimensional space
    %
    % See also:
    %   core.geometry.Geometry, core.geometry.Ball

    properties
        Center % Center
        Radius % Radius
    end

    properties (Dependent)
        BBox % Bounding box
        NSphDims % Number of spherical dimensions (nd - 1)
    end

    methods
        function obj = Sphere(center, radius)
            % SPHERE Construct an instance of Sphere.
            %
            %   obj = Sphere(center, radius) creates a sphere
            %   with specified @a center and @a radius.

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

        function n = get.NSphDims(obj)
            % GET.NSPHDIMS Returns the number of spherical dimensions.

            n = obj.NDims - 1;
        end

        function m = magnitude(obj)
            % MAGNITUDE Returns the surface area of the sphere.
            %
            %   m = magnitude(obj) computes the (n-1)-dimensional surface
            %   area of the sphere using the general formula for
            %   hyperspheres.

            arguments
                obj core.geometry.Sphere
            end

            nsd = obj.NSphDims;

            r = obj.Radius;
            switch nsd
                case 0
                    tol = sqrt(eps);
                    m = 2 * (abs(r) > tol);
                case 1
                    m = 2 * pi * r;
                case 2
                    m = 4 * pi * r^2;
                otherwise
                    m = 2 * pi^((nsd + 1) / 2) / gamma((nsd + 1)/2) * r^nsd;
            end
        end

        function TF = isInside(obj, X)
            % ISINSIDE Tests if points are inside the sphere.
            %
            %   TF = isInside(obj, X) always returns false since a
            %   sphere is a surface with no interior volume.

            arguments
                obj core.geometry.Sphere
                X{mustBeNumeric}
            end

            nd = obj.NDims;
            nx = size(X);
            core.except.assert(nx(1) == nd, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            TF = false(1, nx(2));
        end

        function TF = isOnBoundary(obj, X)
            % ISONBOUNDARY Tests if points are on the sphere surface.
            %
            %   TF = isOnBoundary(obj, X) determines whether points @a X
            %   lie exactly on the sphere surface by checking if their
            %   distance from the center equals the radius.

            arguments
                obj core.geometry.Sphere
                X{mustBeNumeric}
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

        function S = cartesianToSpherical(obj, C)
            % CARTESIANTOSPHERICAL Converts Cartesian to spherical
            % coordinates.
            %
            %   S = cartesianToSpherical(obj, C) transforms points from
            %   Cartesian to spherical coordinates, returning only
            %   angular coordinates since all points have the same radius.

            arguments
                obj core.geometry.Sphere
                C{mustBeNumeric}
            end

            nd = obj.NDims;
            if nd == 1
                S = C;
                return;
            end

            nx = size(C);
            core.except.assert(nx(1) == nd, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            core.except.assert(obj.isOnBoundary(C), ...
                'InvalidInput', ...
                'Input points must be on the sphere surface.');

            C = reshape(C, nx(1), []);
            C = C - obj.Center(:);

            % Compute angular coordinates directly
            S = zeros(nd-1, size(C, 2));
            r = sqrt(sum(C.^2, 1));
            mr = (r ~= 0);

            if any(mr)
                for i = 1:(nd - 2)
                    p = sqrt(sum(C(i:nd, mr).^2, 1));
                    mp = (p ~= 0);
                    if any(mp)
                        ms = mr;
                        ms(mr) = mp;
                        S(i, ms) = acos(C(i, ms) ./ p(mp));
                    end
                end
                a = atan2(C(nd, mr), C(nd-1, mr));
                ma = (a < 0);
                a(ma) = a(ma) + 2 * pi;
                S(nd-1, mr) = a;
            end

            S = reshape(S, [nd - 1, nx(2:end)]);
        end

        function C = sphericalToCartesian(obj, S)
            % SPHERICALTOCARTESIAN Converts spherical to Cartesian
            % coordinates.
            %
            %   C = sphericalToCartesian(obj, S) transforms points from
            %   spherical to Cartesian coordinates, receiving only
            %   angular coordinates and using the sphere's radius
            %   automatically. The resulting points will lie on the
            %   sphere surface.

            arguments
                obj core.geometry.Sphere
                S{mustBeNumeric}
            end

            nd = obj.NDims;
            if nd == 1
                C = S;
                return;
            end

            nx = size(S);
            core.except.assert(nx(1) == nd-1, 'DimensionMismatch', ...
                'Point dimension must match spherical dimension.');

            S = reshape(S, nx(1), []);

            % Validate angle ranges
            isAngleInRange = all(S(1:nd-2, :) >= 0 & S(1:nd-2, :) <= pi);
            core.except.assert(isAngleInRange, 'InvalidInput', ...
                'First nd-2 angles must be in the range [0, pi].');

            isAngleInRange = all(S(nd-1, :) >= 0 & S(nd-1, :) <= 2*pi);
            core.except.assert(isAngleInRange, 'InvalidInput', ...
                'Last angles must be in the range [0, 2*pi].');

            r = obj.Radius;
            A = S;
            C = zeros(nd, size(S, 2));

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
            C = reshape(C, [nd, nx(2:end)]);
        end

    end

    methods (Static)
        function obj = unit(nDims)
            % UNIT Creates a unit sphere centered at origin.
            %
            %   obj = unit(nDims) creates a unit sphere centered at
            %   the origin with radius 1 and dimension @a nDims.

            core.except.assert(nDims >= 1, 'InvalidInput', ...
                'Dimension must be a positive integer.');
            center = zeros(1, nDims);
            radius = 1;
            obj = core.geometry.Sphere(center, radius);
        end
    end
end