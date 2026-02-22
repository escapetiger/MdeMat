classdef SegmentByPlaneIntersector < core.geometry.Geometry
    % SEGMENTBYPLANEINTERSECTOR Computes intersections between line
    % segments and planes.
    %
    %   SegmentByPlaneIntersector computes intersection points between
    %   line segments and planes in n-dimensional space (n >= 2).
    %   Line segments are defined by their endpoints, and planes are
    %   defined by normal vectors and offset constants.
    %
    %   The line segment is parameterized as:
    %   
    %   \f[
    %     L(t) = X1 + t*(X2 - X1),  t \in [0,1]
    %   \f]
    %   
    %   The plane is defined by the equation:
    %
    %   \f[
    %     a^T * x + b = 0,
    %   \f] 
    %   
    %   where \f$a\f$ is the normal vector to the plane and \f$b\f$ is
    %   the offset constant. The method handles parallel and coincident
    %   cases appropriately. If a line segment lies entirely on a plane,
    %   all points are considered intersection points.
    %
    % See also:
    %   core.geometry.Intersector,
    %   core.geometry.SegmentBySegmentIntersector

    methods
        function obj = SegmentByPlaneIntersector(nDims)
            % SEGMENTBYPLANEINTERSECTOR Construct an instance of
            % SegmentByPlaneIntersector.
            %
            %   obj = SegmentByPlaneIntersector(nDims) creates an
            %   intersector for computing line segment-plane intersections
            %   in @a nDims-dimensional space.

            arguments
                nDims {mustBeGreaterThan(nDims, 1)}
            end

            obj@core.geometry.Geometry(nDims);
        end

        function [X, TF] = intersect(obj, X1, X2, a, b, options)
            % INTERSECT Computes intersection between line segments and
            % planes.
            %
            %   [X, TF] = intersect(obj, X1, X2, a, b) computes the
            %   intersection between line segments defined by endpoints @a
            %   X1 - @a X2 and planes defined by normal vectors @a a and
            %   offsets @a b.
            %
            %   [X, TF] = intersect(obj, X1, X2, a, b, tol=tol) uses the
            %   specified numerical tolerance @a tol.

            arguments
                obj core.geometry.SegmentByPlaneIntersector
                X1 {mustBeNumeric}
                X2 {mustBeNumeric}
                a {mustBeNumeric}
                b {mustBeNumeric}
                options.tol {mustBeNumeric} = sqrt(eps)
            end

            core.except.assert(all(size(X1) == size(X2)), ...
                'InvalidInput', 'X1 and X2 must have the same size.');

            nd = obj.NDims;

            nx = size(X1);
            core.except.assert(nx(1) == nd, 'DimensionMismatch', ...
                'Segment dimension must match the intersector dimension.');
            X1 = reshape(X1, nd, []);
            X2 = reshape(X2, nd, []);

            m = size(a);
            core.except.assert(any(a ~= 0), ...
                'InvalidInput', 'Normal vector a cannot be zero.');
            core.except.assert(prod(m(2:end)) == numel(b), ...
                'InvalidInput', ...
                'Normal vectors a and offsets b must have compatible dimensions.');
            core.except.assert(m(1) == nd, 'DimensionMismatch', ...
                'Plane dimension must match the intersector dimension.');
            a = reshape(a, nd, []);
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
            P(abs(U) <= options.tol) = 0;
            TF = (abs(V) > options.tol & P >= 0 & P <= 1);
            TF = TF | (abs(U) <= options.tol & sum(D.^2, 1) < options.tol);
            if sum(TF) == 0
                X = [];
            else
                X = X1(:, TF) + D(:, TF) .* P(TF);
            end
        end
    end
end