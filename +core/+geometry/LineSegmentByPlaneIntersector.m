classdef LineSegmentByPlaneIntersector < core.geometry.Intersector
    % LINESEGMENTBYPLANEINTERSECTOR Computes intersections between line segments and planes.
    %
    %   LineSegmentByPlaneIntersector computes intersection points between
    %   line segments and planes in n-dimensional space (\f$n \ge 2\f$).
    %   Line segments are defined by their endpoints, and planes are
    %   defined by normal vectors and offset constants.
    %
    %   The line segment is parameterized as:
    %   \f[
    %     L(t) = X1 + t*(X2 - X1),  t \in [0,1]
    %   \f]
    %   The plane is defined by the equation:
    %   \f[
    %     a^T * x + b = 0,
    %   \f] 
    %   where \f$a\f$ is the normal vector to the plane and \f$b\f$ is
    %   the offset constant.
    %
    % Examples:
    %   % Create a 3D line-plane intersector
    %   intersector = LineSegmentByPlaneIntersector(3);
    %
    %   % Define a line segment
    %   X1 = [0; 0; 0];   % First endpoint
    %   X2 = [1; 1; 1];   % Second endpoint
    %
    %   % Define a plane: x + y + z - 1 = 0
    %   a = [1; 1; 1];    % Normal vector
    %   b = -1;           % Offset constant
    %
    %   % Compute intersection
    %   [intersectionPoint, isValid] = intersector.intersect(X1, X2, a, b);
    %
    % Notes:
    %   The method handles parallel and coincident cases appropriately.
    %   If a line segment lies entirely on a plane, all points are
    %   considered intersection points.
    %
    % See also:
    %   core.geometry.Intersector,
    %   core.geometry.LineSegmentByLineSegmentIntersector

    methods
        function obj = LineSegmentByPlaneIntersector(nDims)
            % LINESEGMENTBYPLANEINTERSECTOR Constructor for intersector.
            %
            %   obj = LineSegmentByPlaneIntersector(nDims) creates an
            %   intersector for computing line segment-plane intersections
            %   in nDims-dimensional space.
            %
            % Inputs:
            %   nDims - Dimension of the space (must be >= 2)
            % 
            % Outputs:
            %   obj - The constructed LineSegmentByPlaneIntersector object
            %
            % Examples:
            %   % Create intersector for 2D space (line-line intersection)
            %   intersector = LineSegmentByPlaneIntersector(2);
            %
            %   % Create intersector for 3D space
            %   intersector = LineSegmentByPlaneIntersector(3);

            core.except.assert(nDims >= 2, 'InvalidInput', ...
                'Dimension must be no less than 2.');

            obj@core.geometry.Intersector(nDims);
        end

        function [X, TF] = intersect(obj, X1, X2, a, b, tol)
            % INTERSECT Computes intersection between line segments and planes.
            %
            %   [X, TF] = intersect(obj, X1, X2, a, b) computes the
            %   intersection between line segments defined by endpoints
            %   X1-X2 and planes defined by normal vectors a and offsets b.
            %
            %   [X, TF] = intersect(obj, X1, X2, a, b, tol) uses the
            %   specified numerical tolerance.
            %
            % Inputs:
            %   obj - The LineSegmentByPlaneIntersector object
            %   X1  - d×m matrix of first endpoints of line segments
            %   X2  - d×m matrix of second endpoints of line segments
            %   a   - d×k matrix of normal vectors of planes
            %   b   - 1×k vector of offset constants of planes
            %   tol - Numerical tolerance (default: sqrt(eps))
            %
            % Outputs:
            %   X  - d×p matrix of valid intersection points
            %   TF - 1×m logical vector indicating valid intersections
            %
            % Examples:
            %   % Basic intersection computation
            %   [point, hasIntersection] = intersector.intersect(X1, X2, a, b);
            %
            %   % With custom tolerance
            %   [point, hasIntersection] = intersector.intersect(X1, X2, a, b, 1e-10);
            %
            %   % Multiple line segments and planes
            %   X1 = [X1_seg1, X1_seg2, ...];
            %   X2 = [X2_seg1, X2_seg2, ...];
            %   a = [a_plane1, a_plane2, ...];
            %   b = [b_plane1, b_plane2, ...];
            %   [points, validIntersections] = intersector.intersect(X1, X2, a, b);

            if nargin < 6 || isempty(tol)
                tol = sqrt(eps);
            end

            core.except.assert(all(size(X1) == size(X2)), ...
                'InvalidInput', 'X1 and X2 must have the same size.');

            d = obj.nDims;

            n = size(X1);
            core.except.assert(n(1) == d, 'DimensionMismatch', ...
                'Segment dimension must match the intersector dimension.');
            X1 = reshape(X1, d, []);
            X2 = reshape(X2, d, []);

            m = size(a);
            core.except.assert(any(a ~= 0), ...
                'InvalidInput', 'Normal vector a cannot be zero.');
            core.except.assert(prod(m(2:end)) == numel(b), ...
                'InvalidInput', ...
                'Normal vectors a and offsets b must have compatible dimensions.');
            core.except.assert(m(1) == d, 'DimensionMismatch', ...
                'Plane dimension must match the intersector dimension.');
            a = reshape(a, d, []);
            b = reshape(b, 1, []);

            nL = size(X1, 2);
            nP = size(a, 2);
            isOneToOne = nL > 1 & nP > 1;
            D = X2 - X1;
            if isOneToOne
                U = sum(a.*X1+b, 1);
                V = sum(a.*D, 1);
            else
                D = repmat(D, 1, nP);
                X1 = repmat(X1, 1, nP);
                a = repelem(a, 1, nL);
                b = repelem(b, 1, nL);
                U = sum(a.*X1, 1) + b;
                V = sum(a.*D, 1);
            end
            P = -U ./ V;
            P(abs(U) <= tol) = 0;
            TF = (abs(V) > tol & P >= 0 & P <= 1) | (abs(U) <= tol & sum(D.^2, 1) < tol);
            if sum(TF) == 0
                X = [];
            else
                X = X1(:, TF) + D(:, TF) .* P(TF);
            end
        end
    end
end