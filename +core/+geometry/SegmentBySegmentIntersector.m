classdef SegmentBySegmentIntersector < core.geometry.Geometry
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
    %   core.geometry.Intersector,
    %   core.geometry.SegmentByPlaneIntersector

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

            obj@core.geometry.Geometry(nDims);
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
                obj core.geometry.SegmentBySegmentIntersector
                X1 {mustBeNumeric}
                X2 {mustBeNumeric}
                X3 {mustBeNumeric}
                X4 {mustBeNumeric}
                options.tol {mustBeNumeric} = sqrt(eps)
            end

            core.except.assert(all(size(X1) == size(X2) ...
                & size(X2) == size(X3) & size(X3) == size(X4)), ...
                'InvalidInput', ...
                'X1, X2, X3 and X4 must have the same size.');

            nd = obj.NDims;
            nx = size(X1);
            core.except.assert(nx(1) == nd, 'DimensionMismatch', ...
                'Segment dimension must match the intersector dimension.');

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
            P = nan(2, n);
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
                        P(:, i) = [0; 0];
                    else
                        %< At least one is not a point - infinite intersections
                        core.except.assert(0, 'CoincidentSegments', ...
                            ['Line segments coincide exactly, ' ...
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
                        P(:, i) = [0; 0];
                    end
                    continue;
                elseif len1 < tol
                    %< First segment is a point - check if it lies on second segment
                    [onLine, t] = obj.pointOnLineSegment(p1, p3, p4, tol);
                    if onLine
                        TF(i) = true;
                        X(:, i) = p1;
                        P(:, i) = [0; t];
                    end
                    continue;
                elseif len2 < tol
                    %< Second segment is a point - check if it lies on first segment
                    [onLine, t] = obj.pointOnLineSegment(p3, p1, p2, tol);
                    if onLine
                        TF(i) = true;
                        X(:, i) = p3;
                        P(:, i) = [t; 0];
                    end
                    continue;
                end

                %< Handle regular line segments
                %< Check if lines are parallel or collinear

                %< Normalize direction vectors
                u1 = v1 / len1;
                u2 = v2 / len2;

                %< Check if directions are parallel (dot product near ±1)
                a = abs(dot(u1, u2));

                if a > 1 - tol
                    %< Vectors are parallel - check if lines are collinear
                    %< Compute distance from p3 to the line through p1 and p2

                    %< For collinearity check, we need to see if p3 is on the line defined by p1 and p2
                    %< We do this by checking if the vector from p1 to p3 is parallel to v1

                    q = p3 - p1;

                    %< Project connect vector onto v1 and check orthogonal component
                    p = dot(q, u1) * u1;
                    c = q - p;

                    if norm(c) < tol
                        %< Lines are collinear - check for overlap

                        %< Project endpoints of second segment onto first segment
                        t3 = dot(p3-p1, u1) / len1;
                        t4 = dot(p4-p1, u1) / len1;

                        %< Check for single-point intersection cases

                        %< Case 1: p3 is at p1
                        if abs(t3) < tol
                            TF(i) = true;
                            X(:, i) = p1;
                            P(:, i) = [0; 0];
                            continue;
                        end

                        %< Case 2: p3 is at p2
                        if abs(t3-1) < tol
                            TF(i) = true;
                            X(:, i) = p2;
                            P(:, i) = [1; 0];
                            continue;
                        end

                        %< Case 3: p4 is at p1
                        if abs(t4) < tol
                            TF(i) = true;
                            X(:, i) = p1;
                            P(:, i) = [0; 1];
                            continue;
                        end

                        %< Case 4: p4 is at p2
                        if abs(t4-1) < tol
                            TF(i) = true;
                            X(:, i) = p2;
                            P(:, i) = [1; 1];
                            continue;
                        end

                        %< Check for overlap scenarios

                        %< Both endpoints of second segment inside first segment
                        bothInside = (t3 >= 0 && t3 <= 1 && t4 >= 0 && t4 <= 1);

                        %< Second segment contains first segment
                        contains = (t3 < 0 && t4 > 1) || (t3 > 1 && t4 < 0);

                        %< Segment spans across first segment
                        spans = (t3 < 0 && t4 > 0 && t4 <= 1) || (t3 >= 0 && t3 <= 1 && t4 > 1);
                        
                        %< Error case: Multiple intersection points
                        core.except.assert( ...
                            ~bothInside && ~contains && ~spans, ...
                            'CollinearOverlap', ...
                            ['Line segments are collinear and overlapping, ' ...
                            'resulting in more than one intersection point.'])
                    end

                    %< Parallel but not collinear, or collinear with no intersection
                    continue;
                end

                %< Non-parallel lines - find closest points
                %< Set up the coefficient matrix for the least squares problem
                A = [v1, -v2];
                b = p3 - p1;

                %< Solve for parameters [s; t]
                %< We use the Moore-Penrose pseudoinverse for the general case
                theta = pinv(A) * b;
                s = theta(1);
                t = theta(2);

                %< Store parameters
                P(:, i) = [s; t];

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