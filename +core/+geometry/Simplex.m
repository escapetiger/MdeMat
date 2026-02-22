classdef Simplex < core.geometry.Polytope
    % SIMPLEX Geometry for n-dimensional simplex domains.
    %
    %   Simplex represents a generalized geometric domain that extends the
    %   concept of triangles and tetrahedra to arbitrary dimensions. It is
    %   the simplest possible convex polytope in any given dimensional
    %   space.
    %
    %   Dimensional Representations:
    %   - In 1D: A line segment (2 vertices)
    %   - In 2D: A triangle (3 vertices)
    %   - In 3D: A tetrahedron (4 vertices)
    %   - In nD: An n-simplex (n+1 vertices)
    %
    % See also:
    %   core.geometry.Polytope, core.geometry.Orthotope

    properties
        NFaces % Number of boundaries (nDims + 1)
        OuterNormals % Outer normals of boundaries
    end

    methods
        function obj = Simplex(V)
            % SIMPLEX Construct an instance of Simplex.
            %
            %   obj = Simplex(V) creates a simplex from the given vertices,
            %   with strict validation of vertex configuration. The
            %   vertices must be affinely independent and form a valid
            %   simplex. A nd-dimensional simplex requires exactly nd+1
            %   vertices.
            %
            %   The @a V parameter must be a matrix of size nd×(nd+1) where
            %   each column represents a vertex. Vertices must be affinely
            %   independent to ensure the simplex is non-degenerate.

            arguments
                V {mustBeNumeric}
            end

            core.except.assert(ismatrix(V), 'InvalidInput', ...
                'Vertices must be a matrix.');

            [nd, nv] = size(V);
            core.except.assert(nd+1 == nv, 'InvalidInput', ...
                'Vertices must be of shape (%d, %d).', nd, nd+1);

            E = V(:, 2:end) - V(:, 1);
            core.except.assert(rank(E) >= nd, 'InvalidInput', ...
                'Vertices must be affinely independent.');

            F = arrayfun(@(i) setdiff(1:nv, i), 1:nv, 'Un', 0);

            obj@core.geometry.Polytope(V, F);
        end

        function n = get.NFaces(obj)
            % GET.NFACES Returns the number of boundaries.

            n = obj.NDims + 1;
        end

        function N = get.OuterNormals(obj)
            % GET.OUTNORMALS Returns the outward normals of boundaries.

            V = obj.Vertices;
            nd = obj.NDims;
            normal = zeros(nd, nd+1);
            for i = 1:(nd + 1)
                idx = setdiff(1:(nd + 1), i);
                U = V(:, idx);
                U0 = U(:, 1);
                E = U(:, 2:end) - U0;
                ni = null(E.');
                if dot(U0-V(:, i), ni) < 0, ni = -ni; end
                normal(:, i) = ni / norm(ni);
            end
            N = normal;
        end

        function m = magnitude(obj)
            % MAGNITUDE Returns the volume of the simplex.
            %
            %   m = magnitude(obj) computes the n-dimensional volume of the
            %   simplex using the determinant formula. The volume is
            %   calculated as \f$|det(E)|/nd!\f$ where \f$E\f$ is the
            %   matrix of edge vectors from one vertex to all others.

            arguments
                obj core.geometry.Simplex
            end

            nd = obj.NDims;
            V = obj.Vertices;
            E = V(:, 2:end) - V(:, 1);
            m = abs(det(E)) / factorial(nd);
        end

        function TF = isInside(obj, X)
            % ISINSIDE Tests if points are inside the simplex.
            %
            %   TF = isInside(obj, X) uses barycentric coordinates to
            %   determine if points lie strictly within the interior of the
            %   simplex. Points are inside if all barycentric coordinates
            %   are positive (greater than tolerance).

            arguments
                obj core.geometry.Simplex
                X {mustBeNumeric}
            end

            nd = obj.NDims;
            nx = size(X);
            core.except.assert(nx(1) == nd, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            tol = sqrt(eps);
            B = obj.cartesianToBarycentric(X);
            TF = all(B > tol, 1);
        end

        function TF = isOnBoundary(obj, X)
            % ISONBOUNDARY Tests if points are on the simplex boundary.
            %
            %   TF = isOnBoundary(obj, X) uses barycentric coordinates to
            %   determine if points lie exactly on the boundary faces of
            %   the simplex. Points are on the boundary if at least one
            %   barycentric coordinate is zero (within tolerance).

            arguments
                obj core.geometry.Simplex
                X {mustBeNumeric}
            end

            nd = obj.NDims;
            nx = size(X);
            core.except.assert(nx(1) == nd, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            tol = sqrt(eps);
            B = obj.cartesianToBarycentric(X);
            TF = any(abs(B) <= tol, 1);
        end

        function B = cartesianToBarycentric(obj, C)
            % CARTESIANTOBARYCENTRIC Converts Cartesian to barycentric
            % coordinates.
            %
            %   B = cartesianToBarycentric(obj, C) transforms points from
            %   Cartesian coordinates to barycentric coordinates with
            %   respect to the simplex vertices. Barycentric coordinates
            %   sum to 1 and are non-negative for interior points. Each
            %   coordinate represents the weight of the corresponding
            %   vertex.

            arguments
                obj core.geometry.Simplex
                C {mustBeNumeric}
            end

            nd = obj.NDims;
            nx = size(C);
            core.except.assert(nx(1) == nd, 'DimensionMismatch', ...
                'Point dimension must match the simplex dimension.');

            C = reshape(C, nd, []);
            V = obj.Vertices;
            A = [V; ones(1, nd+1)];
            B = [C; ones(1, size(C, 2))];
            B = A \ B;
            nx(1) = nd + 1;
            B = reshape(B, nx);
        end

        function C = barycentricToCartesian(obj, B)
            % BARYCENTRICTOCARTESIAN Converts barycentric to Cartesian
            % coordinates.
            %
            %   C = barycentricToCartesian(obj, B) transforms points from
            %   barycentric coordinates back to Cartesian coordinates. The
            %   barycentric coordinates must sum to 1 and be non-negative
            %   for valid points within or on the simplex.

            arguments
                obj core.geometry.Simplex
                B {mustBeNumeric}
            end

            nd = obj.NDims;
            nx = size(B);
            core.except.assert(nx(1) == nd+1, 'DimensionMismatch', ...
                'Barycentric coordinates must be %d', nd+1);

            tol = sqrt(eps);
            B = reshape(B, nd+1, []);
            core.except.assert(all(abs(sum(B, 1)-1) <= tol), ...
                'InvalidInput', 'Barycentric coordinates must sum to 1.');
            core.except.assert(all(B >= -tol), 'InvalidInput', ...
                'Barycentric coordinates should be non-negative.');
            V = obj.Vertices;
            C = V * B;
            nx(1) = nd;
            C = reshape(C, nx);
        end

        function Y = transform(obj, X, V)
            % TRANSFORM Transform from reference coordinates to physical
            % coordinates.
            %
            %   Y = transform(obj, X, V) transforms points from reference
            %   Cartesian coordinates to physical Cartesian coordinates.

            arguments
                obj core.geometry.Simplex
                X {mustBeNumeric}
                V {mustBeNumeric}
            end

            nd = obj.NDims;
            X = reshape(X, nd, []);
            V0 = obj.Vertices;
            b0 = V0(:, 1);
            A0 = V0(:, 2:end) - b0;
            X = A0 \ (X - b0);
            V = reshape(V, nd, nd+1, []);
            b = V(:, 1, :);
            A = V(:, 2:end, :) - b;
            A = reshape(A, nd, nd, 1, []);
            X = reshape(X, 1, nd, [], 1);
            b = reshape(b, nd, 1, 1, []);
            Y = squeeze(sum(A .* X, 2) + b);
        end
    end
end