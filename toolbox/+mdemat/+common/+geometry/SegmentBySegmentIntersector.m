classdef SegmentBySegmentIntersector < mdemat.data.geometry.Geometry
    % SEGMENTBYSEGMENTINTERSECTOR Computes intersections between
    % line segments.
    %
    %   SegmentBySegmentIntersector computes intersection points
    %   between pairs of line segments in n-dimensional space (n >= 2). The
    %   line segments are defined by their endpoints, and intersections are
    %   computed using a least-squares approach.
    %
    %   The line segments are parameterized as:
    % 
    %   \f[
    %     L1(s) = X1 + s*(X2 - X1)
    %     L2(t) = X3 + t*(X4 - X3)
    %   \f]
    %
    %   An intersection is considered valid only when:
    %     1. Two non-collinear line segments intersect at exactly one point
    %     2. Two collinear line segments intersect at exactly one point
    %
    %   If line segments would have more than one intersection point
    %   (e.g., they overlap), an error is thrown.
    %
    % See also:
    %   mdemat.common.geometry.SegmentByPlaneIntersector

    methods
        function obj = SegmentBySegmentIntersector(nDims)
            % SEGMENTBYSEGMENTINTERSECTOR Construct an instance of
            % SegmentBySegmentIntersector.
            %
            %   obj = SegmentBySegmentIntersector(nDims) creates
            %   an intersector for line segment intersections in @a nDims-
            %   dimensional space.

            arguments
                nDims {mustBeGreaterThan(nDims, 1)}
            end

            obj@mdemat.data.geometry.Geometry(nDims);
        end

        function [X, TF] = intersect(obj, X1, X2, X3, X4, options)
            % INTERSECT Computes intersection between line segments.
            %
            %   [X, TF] = intersect(obj, X1, X2, X3, X4) computes the
            %   intersection between line segments X1-X2 and X3-X4.
            %
            %   [X, TF] = intersect(obj, X1, X2, X3, X4, tol=tol) uses the
            %   specified numerical @a tol tolerance.

            arguments
                obj mdemat.common.geometry.SegmentBySegmentIntersector
                X1 {mustBeNumeric}
                X2 {mustBeNumeric}
                X3 {mustBeNumeric}
                X4 {mustBeNumeric}
                options.tol (1, 1) {mustBeReal, mustBeFinite, mustBePositive} = sqrt(eps)
            end

            assert(all(all(size(X1) == size(X2) ...
                & size(X2) == size(X3) & size(X3) == size(X4)), 'all'), ...
                'mdemat:common:geometry:SegmentBySegmentIntersector:InvalidInput', ...
                '%s', 'X1, X2, X3 and X4 must have the same size.');

            nd = obj.NDims;
            nx = size(X1);
            assert(all(nx(1) == nd, 'all'), 'mdemat:common:geometry:SegmentBySegmentIntersector:DimensionMismatch', ...
                '%s', 'Segment dimension must match the intersector dimension.');

            [X, TF] = obj.intersectND(X1, X2, X3, X4, options.tol);
        end
    end

    methods (Access = protected)
        function [X, TF] = intersectND(obj, X1, X2, X3, X4, tol)
            % INTERSECTND Internal method for n-dimensional intersection
            % computation.
            
            [d, n] = size(X1);
            X = nan(d, n);
            TF = false(1, n);
            for i = 1:n
                p1 = X1(:, i);
                p2 = X2(:, i);
                p3 = X3(:, i);
                p4 = X4(:, i);
                v1 = p2 - p1;
                v2 = p4 - p3;
                len1 = norm(v1);
                len2 = norm(v2);

                %< First, check if line segments coincide exactly (same endpoints)
                if (norm(p1-p3) < tol && norm(p2-p4) < tol) || (norm(p1-p4) < tol && norm(p2-p3) < tol)
                    %< Line segments have the same endpoints (in same or reversed order)
                    if len1 < tol && len2 < tol
                        %< Both are points - valid single intersection
                        TF(i) = true;
                        X(:, i) = p1;
                    else
                        %< At least one is not a point - infinite intersections
                        assert(all(0, 'all'), 'mdemat:common:geometry:SegmentBySegmentIntersector:CoincidentSegments', ...
                            '%s', ['Line segments coincide exactly, ' ...
                            'resulting in more than one intersection point.']);
                    end
                    continue;
                end

                %< Handle zero-length segments (points)
                if len1 < tol && len2 < tol
                    %< Both segments are points - check if they're the same point
                    if norm(p1-p3) < tol
                        TF(i) = true;
                        X(:, i) = p1;
                    end
                    continue;
                elseif len1 < tol
                    %< First segment is a point - check if it lies on second segment
                    onLine = obj.pointOnLineSegment(p1, p3, p4, tol);
                    if onLine
                        TF(i) = true;
                        X(:, i) = p1;
                    end
                    continue;
                elseif len2 < tol
                    %< Second segment is a point - check if it lies on first segment
                    onLine = obj.pointOnLineSegment(p3, p1, p2, tol);
                    if onLine
                        TF(i) = true;
                        X(:, i) = p3;
                    end
                    continue;
                end

                %< Handle regular line segments
                %< Check if lines are parallel or collinear

                %< Normalize direction vectors
                u1 = v1 / len1;
                u2 = v2 / len2;

                % The perpendicular endpoint displacement determines whether
                % the directions are indistinguishable at the distance tol.
                alignment = dot(u1, u2);
                parallelError = norm(u1 - alignment*u2);
                if parallelError <= tol/max(len1, len2)
                    separation = p3-p1;
                    separation = separation-dot(separation, u1)*u1;
                    if norm(separation) > tol
                        continue;
                    end

                    % Intersect the projected closed intervals before testing
                    % endpoints: a shared endpoint can also bound an overlap.
                    t3 = dot(p3-p1, u1)/len1;
                    t4 = dot(p4-p1, u1)/len1;
                    lower = max(0, min(t3, t4));
                    upper = min(1, max(t3, t4));
                    parameterTolerance = tol/len1;
                    if upper < lower-parameterTolerance
                        continue;
                    end
                    assert(all(upper-lower <= parameterTolerance, 'all'), ...
                        'mdemat:common:geometry:SegmentBySegmentIntersector:CollinearOverlap', ...
                        '%s', 'Collinear segments overlap in more than one point.');
                    parameter = min(1, max(0, (lower+upper)/2));
                    TF(i) = true;
                    X(:, i) = p1+parameter*v1;
                    continue;
                end

                %< Non-parallel lines - find closest points
                %< Set up the coefficient matrix for the least squares problem
                A = [v1, -v2];
                b = p3 - p1;

                %< Solve for parameters [s; t]
                % Full-column-rank least squares for nonparallel directions.
                theta = A \ b;
                s = theta(1);
                t = theta(2);

                %< Store parameters

                %< Compute residual (distance between closest points)
                Y1 = p1 + s * v1;
                Y2 = p3 + t * v2;

                %< Check if closest points are within tolerance and on segments
                if norm(Y1-Y2) <= tol && s >= 0 && s <= 1 && t >= 0 && t <= 1
                    TF(i) = true;

                    %< Use midpoint between closest points as intersection
                    X(:, i) = (Y1 + Y2) / 2;
                end
            end
            X = X(:, TF);
        end

        function [onLine, t] = pointOnLineSegment(~, point, lineStart, lineEnd, tol)
            % POINTONLINESEGMENT Determines if a point lies on a line
            % segment.

            %< If tolerance not provided, use default
            if nargin < 5 || isempty(tol)
                tol = sqrt(eps);
            end

            %< Get line direction
            lineDir = lineEnd - lineStart;
            lineLength = norm(lineDir);

            %< If line has zero length, check if point is at start
            if lineLength < tol
                onLine = norm(point-lineStart) < tol;
                t = 0;
                return;
            end

            %< Compute projection parameter
            t = dot(point-lineStart, lineDir) / (lineLength^2);

            %< Check if parameter is in [0,1]
            if t < 0 || t > 1
                onLine = false;
                return;
            end

            %< Compute closest point on line
            closestPoint = lineStart + t * lineDir;

            %< Check if point is close enough to line
            onLine = norm(point-closestPoint) < tol;
        end
    end
end
