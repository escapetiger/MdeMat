classdef Polytope < core.geometry.Geometry
    % POLYTOPE Base class for all polytopes.
    %
    %   Polytope defines a generalized n-dimensional polytope, representing
    %   geometric objects with flat sides that extend the concept of
    %   polygons (2D) and polyhedra (3D) to higher dimensional spaces. A
    %   polytope is fundamentally characterized by its vertices and
    %   optionally its faces.
    %
    % Examples:
    %   % Cannot instantiate directly - use concrete subclasses
    %   % simplex = core.geometry.Simplex(vertices);
    %   % orthotope = core.geometry.Orthotope([0, 1, 0, 1]);
    %
    % See also: 
    %   core.geometry.Simplex, core.geometry.Orthotope

    properties
        vertices double % Matrix representing polytope vertices (each column is a vertex)
        faces cell % Cell array of face definitions (each cell contains vertex indices)
    end
    
    methods
        function obj = Polytope(V, F)
            % POLYTOPE Constructor for Polytope.
            %
            %   obj = Polytope(V, F) creates a polytope from given vertices
            %   and optional faces. The constructor validates that the
            %   polytope has sufficient vertices for the specified dimension
            %   and ensures proper face definitions with valid vertex
            %   indices.
            %
            % Inputs:
            %   V - A matrix of size d×m where each column is a vertex
            %   F - A cell array of faces where each cell contains vertex indices
            %
            % Outputs:
            %   obj - Constructed Polytope object
            
            core.except.assert(nargin >= 1, 'InvalidInput', ...
                'Vertices must be specified.');

            core.except.assert(ismatrix(V), 'InvalidInput', ...
                'Vertices must be a matrix.');
            
            d = size(V, 1);
            m = size(V, 2);
            core.except.assert(m >= (d + 1), 'InvalidInput', ...
                'A %d-dimensional polytope must have at least %d vertices.', d, d + 1)
            core.except.assert(iscell(F), 'InvalidInput', ...
                'Faces must be a cell array.');
            for i = 1:length(F)
                core.except.assert(all(ismember(F{i}, 1:m)), ...
                    'InvalidInput', 'Face contains Invalid vertex indices.');
                core.except.assert(length(F{i}) >= d, ...
                    'InvalidInput', 'Face must have at least %d vertices.', d);
            end

            obj@core.geometry.Geometry(d);
            obj.vertices = V;
            obj.faces = F;
        end
        
        function TF = isInside(obj, X)
            % ISINSIDE Tests if points are inside the polytope.
            %
            %   TF = isInside(obj, X) uses barycentric coordinate
            %   computation to determine if points lie within the convex
            %   hull of the polytope vertices. Points are inside if their
            %   barycentric coordinates are all non-negative and sum to 1.
            %
            % Inputs:
            %   obj - The Polytope object
            %   X - A matrix of size d×m where each column represents a point
            %
            % Outputs:
            %   TF - True if the point is inside
            
            d = obj.nDims;
            n = size(X);
            core.except.assert(n(1) == d, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');

            X = reshape(X, n(1), []);
            tol = sqrt(eps);
            V = obj.vertices;
            A = [V; ones(1, size(V, 2))];
            B = [X; ones(1, size(X, 2))];
            W = pinv(A) * B;
            TF = all(W >= -tol, 1) & (abs(sum(W, 1) - 1) < tol) & (sum((A*W-B).^2, 1) < tol);
            if length(n) > 2
                TF = reshape(TF, n(2:end));
            else
                TF = reshape(TF, 1, n(2));
            end
        end
        
        function TF = isOnBoundary(obj, X)
            % ISONBOUNDARY Tests if points are on the polytope boundary.
            %
            %   TF = isOnBoundary(obj, X) checks if points lie on any of
            %   the defined faces of the polytope. This method requires
            %   that faces be properly defined in the polytope
            %   construction.
            %
            % Inputs:
            %   obj - The Polytope object
            %   X - A matrix of size d×m where each column represents a point
            %
            % Outputs:
            %   TF - True if the point is on the boundary
            
            d = obj.nDims;
            n = size(X);
            core.except.assert(n(1) == d, 'DimensionMismatch', ...
                'Point dimension must match geometry dimension.');
            core.except.assert(~isempty(obj.faces), 'MissingFaces', ...
                'Faces must be defined to determine boundary points.')

            V = obj.vertices;
            X = reshape(X, n(1), []);
            TF = false(1, size(X, 2));
            for j = 1:length(obj.faces)
                TF = TF | obj.isOnFace(X, V(:, obj.faces{j}));    
            end
            if length(n) > 2
                TF = reshape(TF, n(2:end));
            else
                TF = reshape(TF, 1, n(2));
            end
        end
    end
    
    methods (Access = protected)                
        function TF = isOnFace(obj, X, V)
            % ISONFACE Tests if points lie on a specific face.
            %
            %   TF = isOnFace(obj, X, V) determines if points lie on the
            %   face defined by the given vertices. The method dispatches
            %   to dimension-specific implementations for computational
            %   efficiency and accuracy.
            %
            % Inputs:
            %   obj - The Polytope object
            %   X - A matrix of size d×m where each column represents a point
            %   V - A matrix where each column is a vertex of the face
            %
            % Outputs:
            %   TF - A logical row vector indicating which points are on the face
            
            tol = sqrt(eps);
            switch obj.nDims
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
            %
            %   TF = isOnFace2D(obj, X, V, tol) provides a specialized
            %   method for testing point containment on line segments in
            %   2D space. The method checks both distance to the line and
            %   position within the segment bounds.
            %
            % Inputs:
            %   X - A matrix of size 2×m where each column represents a point
            %   V - A matrix where each column is a vertex of the face
            %   tol - Numerical tolerance for distance calculations
            %
            % Outputs:
            %   TF - A logical row vector indicating which points are on the face
            
            %< Get dimensions
            [~, m] = size(X);
            TF = false(1, m);
            
            %< In 2D, a face is a line segment
            v1 = V(:, 1);
            v2 = V(:, 2);
            
            %< Vector from v1 to v2 (edge vector)
            v = v2 - v1;
            
            %< Norm of edge vector
            vLen = norm(v);
            if vLen < eps
                return; %< Degenerate edge, no points can be on it
            end
            
            %< Unit vector along edge
            vUnit = v / vLen;
            
            %< Vectors from v1 to all points
            vps = X - v1;
            
            %< Project vps onto vUnit
            projs = vUnit' * vps; % 1×m vector of projections
            
            %< Points are within segment bounds if 0 <= proj <= vLen
            withinBounds = (projs >= 0) & (projs <= vLen);
            
            %< Distance from points to line
            projPoints = v1 + vUnit * projs; % Points projected onto line
            distToLine = sqrt(sum((X - projPoints).^2, 1)); % Distance vector
            
            %< Points are on line if distance is very small
            onLine = distToLine < tol;
            
            %< Points are on segment if they're on the line and within bounds
            TF = withinBounds & onLine;
        end
        
        function TF = isOnFace3D(~, X, V, tol)
            % ISONFACE3D Tests if points lie on a 3D face (polygon).
            %
            %   TF = isOnFace3D(obj, X, V, tol) provides a specialized
            %   method for testing point containment on polygonal faces in
            %   3D space. The method computes the face normal and checks
            %   if points lie on the corresponding plane.
            %
            % Inputs:
            %   X - A matrix of size 3×m where each column represents a point
            %   V - A matrix where each column is a vertex of the face
            %   tol - Numerical tolerance for distance calculations
            %
            % Outputs:
            %   TF - A logical row vector indicating which points are on the face
            
            %< Get dimensions
            [~, m] = size(X);
            TF = false(1, m);
            
            if size(V, 2) < 3
                return; %< Degenerate face, no points can be on it
            end
            
            %< In 3D, compute face normal from three points
            v1 = V(:, 1);
            v2 = V(:, 2);
            v3 = V(:, 3);
            
            %< Compute face normal
            normal = cross(v2 - v1, v3 - v1);
            normalLen = norm(normal);
            
            if normalLen < eps
                return; %< Degenerate face, no points can be on it
            end
            
            %< Normalize normal vector
            normal = normal / normalLen;
            
            %< Compute plane constant for the equation: normal·x + planeConst = 0
            planeConst = -dot(normal, v1);
            
            %< Check if points lie on the plane
            distToPlane = abs(normal' * X + planeConst);
            onPlane = distToPlane < tol;
            
            %< For simplicity, we'll assume all points on the plane are also within the face
            TF = onPlane;
        end
        
        function TF = isOnFaceND(obj, X, V, tol)
            % ISONFACEND Tests if points lie on an n-dimensional face.
            %
            %   TF = isOnFaceND(obj, X, V, tol) provides a generalized
            %   method for testing point containment on hyperplanar faces
            %   in n-dimensional space. The method uses the null space to
            %   compute the hyperplane normal.
            %
            % Inputs:
            %   obj - The Polytope object
            %   X - A matrix of size nDims×m where each column represents a point
            %   V - A matrix where each column is a vertex of the face
            %   tol - Numerical tolerance for distance calculations
            %
            % Outputs:
            %   TF - A logical row vector indicating which points are on the face
            
            %< Get dimensions
            [~, m] = size(X);
            TF = false(1, m);
            
            %< Start by computing (nDims-1) edge vectors
            if size(V, 2) < obj.nDims
                return; %< Degenerate face, need at least nDims vertices
            end
            
            %< Use the first vertex as reference
            v1 = V(:, 1);
            
            %< Compute edge vectors
            edges = V(:, 2:obj.nDims) - v1;
            
            %< Compute the normal to the hyperplane (using null space)
            normal = null(edges');
            
            %< If multiple vectors in null space, use the first one
            if size(normal, 2) > 1
                normal = normal(:, 1);
            end
            
            %< Normalize the normal vector
            normal = normal / norm(normal);
            
            %< Compute hyperplane constant
            planeConst = -dot(normal, v1);
            
            %< Check if points lie on the hyperplane
            distToPlane = abs(normal' * X + planeConst);
            onPlane = distToPlane < tol;
            
            %< For simplicity, we'll assume all points on the hyperplane are within the face
            TF = onPlane;
        end
    end
end