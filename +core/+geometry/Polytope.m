classdef Polytope < core.geometry.Geometry
    % POLYTOPE Base class for all polytopes.
    %
    %   Polytope defines a generalized n-dimensional polytope, representing
    %   geometric objects with flat sides that extend the concept of
    %   polygons (2D) and polyhedra (3D) to higher dimensional spaces. A
    %   polytope is fundamentally characterized by its vertices and
    %   optionally its faces.
    %
    % See also:
    %   core.geometry.Simplex, core.geometry.Orthotope

    properties
        Vertices % Polytope vertices (each column is a vertex)
        Faces % Cell array of faces (each cell contains vertex indices)
    end
    
    methods
        function obj = Polytope(V, F)
            % POLYTOPE Construct an instance of Polytope.
            %
            %   obj = Polytope(V, F) creates a polytope from given @a V
            %   vertices and @a F faces with validation of vertex
            %   configuration and face definitions.
            
            arguments
                V {mustBeNumeric}
                F {mustBeA(F, "cell")}
            end
            
            core.except.assert(ismatrix(V), 'InvalidInput', ...
                'Vertices must be a matrix.');
            
            nd = size(V, 1);
            nv = size(V, 2);
            core.except.assert(nv >= (nd + 1), 'InvalidInput', ...
                'A %d-dimensional polytope must have at least %d vertices.', nd, nd + 1)
            for i = 1:length(F)
                core.except.assert(all(ismember(F{i}, 1:nv)), ...
                    'InvalidInput', 'Face contains invalid vertex indices.');
                core.except.assert(length(F{i}) >= nd, 'InvalidInput', ...
                    'Face must have at least %d vertices.', nd);
            end

            obj@core.geometry.Geometry(nd);
            obj.Vertices = V;
            obj.Faces = F;
        end
        
        function TF = isInside(obj, X)
            % ISINSIDE Tests if points are inside the polytope.
            %
            %   TF = isInside(obj, X) uses barycentric coordinate
            %   computation to determine if points @a X lie within the convex
            %   hull of the polytope vertices. @a X should be an nd-by-np matrix.
            
            arguments
                obj core.geometry.Polytope
                X {mustBeNumeric}
            end
            
            nd = obj.NDims;
            nx = size(X);
            core.except.assert(nx(1) == nd, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            X = reshape(X, nx(1), []);
            tol = sqrt(eps);
            V = obj.Vertices;
            A = [V; ones(1, size(V, 2))];
            B = [X; ones(1, size(X, 2))];
            W = pinv(A) * B;
            TF = all(W >= -tol, 1) & (abs(sum(W, 1) - 1) < tol) & (sum((A*W-B).^2, 1) < tol);
            if length(nx) > 2
                TF = reshape(TF, nx(2:end));
            else
                TF = reshape(TF, 1, nx(2));
            end
        end
        
        function TF = isOnBoundary(obj, X)
            % ISONBOUNDARY Tests if points are on the polytope boundary.
            %
            %   TF = isOnBoundary(obj, X) checks if points @a X lie on any of
            %   the defined faces of the polytope. @a X should be an nd-by-np matrix.
            
            arguments
                obj core.geometry.Polytope
                X {mustBeNumeric}
            end
            
            nd = obj.NDims;
            nx = size(X);
            core.except.assert(nx(1) == nd, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');
            core.except.assert(~isempty(obj.Faces), 'MissingFaces', ...
                'Faces must be defined to determine boundary points.')

            V = obj.Vertices;
            X = reshape(X, nx(1), []);
            TF = false(1, size(X, 2));
            for j = 1:length(obj.Faces)
                TF = TF | obj.isOnFace(X, V(:, obj.Faces{j}));    
            end
            if length(nx) > 2
                TF = reshape(TF, nx(2:end));
            else
                TF = reshape(TF, 1, nx(2));
            end
        end
    end
    
    methods (Access = protected)                
        function TF = isOnFace(obj, X, V)
            % ISONFACE Tests if points lie on a specific face.
            
            nd = obj.NDims;
            tol = sqrt(eps);
            switch nd
                case 2
                    TF = obj.isOnFace2D(X, V, tol);
                case 3
                    TF = obj.isOnFace3D(X, V, tol);
                otherwise
                    TF = obj.isOnFaceND(X, V, tol);
            end
        end
        
        function TF = isOnFace2D(~, X, V, tol)
            % ISONFACE2D Tests if points lie on a 2D face (line segment).
            
            np = size(X, 2);
            TF = false(1, np);
            
            v1 = V(:, 1);
            v2 = V(:, 2);
            
            v = v2 - v1;
            
            vLen = norm(v);
            if vLen < eps
                return;
            end
            
            vUnit = v / vLen;
            
            vps = X - v1;
            
            projs = vUnit' * vps;
            
            withinBounds = (projs >= 0) & (projs <= vLen);
            
            projPoints = v1 + vUnit * projs;
            distToLine = sqrt(sum((X - projPoints).^2, 1));
            
            onLine = distToLine < tol;
            
            TF = withinBounds & onLine;
        end
        
        function TF = isOnFace3D(~, X, V, tol)
            % ISONFACE3D Tests if points lie on a 3D face (polygon).
            
            np = size(X, 2);
            TF = false(1, np);
            
            if size(V, 2) < 3
                return;
            end
            
            v1 = V(:, 1);
            v2 = V(:, 2);
            v3 = V(:, 3);
            
            normal = cross(v2 - v1, v3 - v1);
            normalLen = norm(normal);
            
            if normalLen < eps
                return;
            end
            
            normal = normal / normalLen;
            
            planeConst = -dot(normal, v1);
            
            distToPlane = abs(normal' * X + planeConst);
            onPlane = distToPlane < tol;
            
            TF = onPlane;
        end
        
        function TF = isOnFaceND(obj, X, V, tol)
            % ISONFACEND Tests if points lie on an n-dimensional face.
            
            nd = obj.NDims;
            np = size(X, 2);
            TF = false(1, np);
            
            if size(V, 2) < nd
                return;
            end
            
            v1 = V(:, 1);
            
            edges = V(:, 2:nd) - v1;
            
            normal = null(edges');
            
            if size(normal, 2) > 1
                normal = normal(:, 1);
            end
            
            normal = normal / norm(normal);
            
            planeConst = -dot(normal, v1);
            
            distToPlane = abs(normal' * X + planeConst);
            onPlane = distToPlane < tol;
            
            TF = onPlane;
        end
    end
end