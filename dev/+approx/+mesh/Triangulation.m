classdef Triangulation < approx.mesh.PolytopalMesh
    % TRIANGULATION Computational mesh using simplicial elements.
    %
    %   Triangulation provides a concrete implementation of PolytopalMesh
    %   for simplicial elements (line segments in 1D, triangles in 2D,
    %   tetrahedra in 3D).
    %
    %   Element types by dimension:
    %   - 1D: Line segments (2 vertices, 2 faces as endpoints)
    %   - 2D: Triangles (3 vertices, 3 edges as faces)
    %   - 3D: Tetrahedra (4 vertices, 4 triangular faces)
    %   - nD: n-simplices (n+1 vertices, n+1 (n-1)-simplex faces)
    %
    %   The class uses affine coordinate transformations between reference
    %   and physical coordinates, where the reference simplex has vertices
    %   at the origin and unit coordinate vectors: {0, e_1, e_2, ..., e_n}.
    %
    %   Reference coordinates use the standard simplex with vertices at {0,
    %   e_1, e_2, ..., e_n} where e_i are unit coordinate vectors. This
    %   differs from some references that use centered or symmetric
    %   simplices.
    %
    %   Mesh refinement follows standard subdivision rules:
    %   - 1D: Each edge splits into 2 segments
    %   - 2D: Each triangle splits into 4 triangles
    %   - 3D: Each tetrahedron splits into 8 tetrahedra
    %   - General: Each k-level refinement produces 2^(k×d) sub-elements
    %
    % See also:
    %   approx.mesh.PolytopalMesh, approx.mesh.PolytopalMeshPlotter

    methods
        function graphObj = graphify(obj)
            % GRAPHIFY Generate graph representation from simplex
            % connectivity.
            %
            %   graphObj = graphify(obj) creates a graph representation by
            %   extracting all edges (vertex pairs) from the simplex
            %   elements and eliminating duplicates.

            arguments
                obj approx.mesh.Triangulation
            end

            E = obj.Elements.';
            T = nchoosek(1:obj.NVerticesPerElement, 2).';
            edges = reshape(E(T(:), :), 2, []);
            edges = unique(sort(edges, 1).', 'rows', 'stable');
            graphObj = approx.mesh.Graph(obj.Vertices, edges);
        end

        function Y = collocate(obj, X, EI)
            % COLLOCATE Map reference coordinates to physical coordinates.
            %
            %   Y = collocate(obj, X, EI) transforms reference Cartesian
            %   coordinates to physical Cartesian coordinates using affine
            %   transformations defined by simplex element vertices.

            arguments
                obj approx.mesh.Triangulation
                X {mustBeNumeric}
                EI {mustBeInteger} = []
            end

            if isempty(EI)
                EI = 1:obj.NElements;
            end
            nd = obj.NDims;
            V = obj.findElementVertices(EI);
            X = reshape(X, nd, []);
            V = reshape(V, nd, nd+1, []);
            b = V(:, 1, :);
            A = V(:, 2:end, :) - b;
            A = reshape(A, nd, nd, 1, []);
            X = reshape(X, 1, nd, [], 1);
            b = reshape(b, nd, 1, 1, []);
            Y = squeeze(sum(A .* X, 2) + b);
        end

        function newObj = refine(obj, nLevels)
            % REFINE Perform uniform subdivision of simplicial mesh.
            %
            %   newObj = refine(obj, nLevels) creates a refined mesh by
            %   uniformly subdividing each simplex element. The number of
            %   refinement levels determines the subdivision density.

            arguments
                obj approx.mesh.Triangulation
                nLevels {mustBeInteger, mustBeNonnegative} = 0
            end

            if nLevels == 0
                newObj = obj;
                return;
            end

            V = obj.Vertices;
            E = obj.Elements;
            for i = 1:nLevels
                switch obj.NDims
                    case 1
                        [V, E] = refineLines(V, E);
                    case 2
                        [V, E] = refineTriangles(V, E);
                    case 3
                        [V, E] = refineTetrahedra(V, E);
                end
            end
            newObj = approx.mesh.Triangulation(V, E);
        end

        function h = computeMeasure(obj)
            % COMPUTEMEASURE Compute characteristic mesh size.
            %
            %   h = computeMeasure(obj) returns the minimum edge length
            %   among all edges in the simplicial mesh, which serves as
            %   a characteristic measure of mesh resolution.

            arguments
                obj approx.mesh.Triangulation
            end

            nd = obj.NDims;
            T = nchoosek(1:(nd + 1), 2).';
            E = obj.Elements(:, T(:)).';
            V = reshape(obj.Vertices(E(:), :), 2, [], nd);
            edgeLengths = sqrt(sum(squeeze(diff(V, 1, 1)).^2, 2));
            h = min(edgeLengths, [], 'all');
        end

        function EJac = computeAllElementJacobians(obj)
            % COMPUTEALLELEMENTJACOBIANS Compute Jacobian matrices for all
            % elements.
            %
            %   EJac = computeAllElementJacobians(obj) computes the Jacobian
            %   matrix for coordinate transformation from reference to
            %   physical coordinates for every element in the mesh.

            arguments
                obj approx.mesh.Triangulation
            end

            V = obj.findElementVertices();
            EJac = permute(V(:, 2:end, :)-V(:, 1, :), [2, 1, 3]);
        end

        function detEJac = computeAllElementJacobianDeterminants(obj)
            % COMPUTEALLELEMENTJACOBIANDETERMINANTS Compute Jacobian
            % determinants for all elements.
            %
            %   detEJac = computeAllElementJacobianDeterminants(obj) computes
            %   the determinant of the Jacobian matrix for each simplex
            %   element, used in integration and volume calculations.

            arguments
                obj approx.mesh.Triangulation
            end

            EJac = obj.computeElementJacobians();
            switch obj.NDims
                case 1
                    detEJac = EJac;
                case 2
                    EJac11 = squeeze(EJac(1, 1, :));
                    EJac12 = squeeze(EJac(1, 2, :));
                    EJac21 = squeeze(EJac(2, 1, :));
                    EJac22 = squeeze(EJac(2, 2, :));
                    detEJac = EJac11 .* EJac22 - EJac12 .* EJac21;
                case 3
                    EJac11 = squeeze(EJac(1, 1, :));
                    EJac12 = squeeze(EJac(1, 2, :));
                    EJac13 = squeeze(EJac(1, 3, :));
                    EJac21 = squeeze(EJac(2, 1, :));
                    EJac22 = squeeze(EJac(2, 2, :));
                    EJac23 = squeeze(EJac(2, 3, :));
                    EJac31 = squeeze(EJac(3, 1, :));
                    EJac32 = squeeze(EJac(3, 2, :));
                    EJac33 = squeeze(EJac(3, 3, :));
                    detEJac = EJac11 .* (EJac22 .* EJac33 - EJac23 .* EJac32) ...
                        -EJac12 .* (EJac21 .* EJac33 - EJac23 .* EJac31) ...
                        +EJac13 .* (EJac21 .* EJac32 - EJac22 .* EJac31);

                otherwise
                    detEJac = arrayfun(@(i) det(EJac(:, :, i)), 1:size(EJac, 3)).';
            end
        end

        function invEJac = computeAllElementInverseJacobians(obj)
            % COMPUTEALLELEMENTINVERSEJACOBIANS Compute inverse Jacobian
            % matrices for all elements.
            %
            %   invEJac = computeAllElementInverseJacobians(obj) computes the
            %   inverse of the Jacobian matrix for each element, used in
            %   gradient computations and coordinate transformations.

            arguments
                obj approx.mesh.Triangulation
            end

            switch obj.NDims
                case 1
                    detEJac = obj.computeElementJacobianDeterminants();
                    invEJac = 1 ./ detEJac;

                case 2
                    EJac = obj.computeElementJacobians();
                    detEJac = obj.computeElementJacobianDeterminants();

                    EJac11 = squeeze(EJac(1, 1, :));
                    EJac12 = squeeze(EJac(1, 2, :));
                    EJac21 = squeeze(EJac(2, 1, :));
                    EJac22 = squeeze(EJac(2, 2, :));

                    invEJac = zeros(size(EJac));
                    invEJac(1, 1, :) = EJac22 ./ detEJac;
                    invEJac(1, 2, :) = -EJac12 ./ detEJac;
                    invEJac(2, 1, :) = -EJac21 ./ detEJac;
                    invEJac(2, 2, :) = EJac11 ./ detEJac;

                case 3
                    EJac = obj.computeElementJacobians();
                    detEJac = obj.computeElementJacobianDeterminants();

                    EJac11 = squeeze(EJac(1, 1, :));
                    EJac12 = squeeze(EJac(1, 2, :));
                    EJac13 = squeeze(EJac(1, 3, :));
                    EJac21 = squeeze(EJac(2, 1, :));
                    EJac22 = squeeze(EJac(2, 2, :));
                    EJac23 = squeeze(EJac(2, 3, :));
                    EJac31 = squeeze(EJac(3, 1, :));
                    EJac32 = squeeze(EJac(3, 2, :));
                    EJac33 = squeeze(EJac(3, 3, :));

                    invEJac = zeros(size(EJac));
                    invEJac(1, 1, :) = (EJac22 .* EJac33 - EJac23 .* EJac32) ./ detEJac;
                    invEJac(1, 2, :) = (EJac13 .* EJac32 - EJac12 .* EJac33) ./ detEJac;
                    invEJac(1, 3, :) = (EJac12 .* EJac23 - EJac13 .* EJac22) ./ detEJac;
                    invEJac(2, 1, :) = (EJac23 .* EJac31 - EJac21 .* EJac33) ./ detEJac;
                    invEJac(2, 2, :) = (EJac11 .* EJac33 - EJac13 .* EJac31) ./ detEJac;
                    invEJac(2, 3, :) = (EJac13 .* EJac21 - EJac11 .* EJac23) ./ detEJac;
                    invEJac(3, 1, :) = (EJac21 .* EJac32 - EJac22 .* EJac31) ./ detEJac;
                    invEJac(3, 2, :) = (EJac12 .* EJac31 - EJac11 .* EJac32) ./ detEJac;
                    invEJac(3, 3, :) = (EJac11 .* EJac22 - EJac12 .* EJac21) ./ detEJac;

                otherwise
                    EJac = obj.computeElementJacobians();
                    invEJac = zeros(size(EJac));
                    for i = 1:size(EJac, 3)
                        invEJac(:, :, i) = inv(EJac(:, :, i));
                    end
            end
        end

        function detFJac = computeAllFaceJacobianDeterminants(obj)
            % COMPUTEALLFACEJACOBIANDETERMINANTS Compute face Jacobian
            % determinants for all faces.
            %
            %   detFJac = computeAllFaceJacobianDeterminants(obj) computes
            %   Jacobian determinants for coordinate transformation on all
            %   element faces, used in boundary and interface integration.

            arguments
                obj approx.mesh.Triangulation
            end

            nd = obj.NDims;
            detFJac = cell(1, nd+1);

            for i = 1:nd + 1
                V = obj.findElementVertices();
                U = V(:, setdiff(1:nd+1, i), :);
                G = U(:, 2:end, :) - U(:, 1, :);
                G = pagemtimes(permute(G, [2, 1, 3]), G);
                gramDets = arrayfun(@(k) det(G(:, :, k)), 1:obj.NElements);
                detFJac{i} = reshape(sqrt(gramDets), [], 1);
            end
        end

        function normal = computeAllOutwardNormals(obj)
            % COMPUTEALLOUTWARDNORMALS Compute outward normal vectors for
            % all faces.
            %
            %   normal = computeAllOutwardNormals(obj) computes unit
            %   outward normal vectors for all element faces, essential for
            %   boundary conditions and flux computations.

            arguments
                obj approx.mesh.Triangulation
            end

            nd = obj.NDims;
            normal = cell(1, nd+1);

            for i = 1:nd + 1
                V = obj.findElementVertices();
                U = V(:, setdiff(1:nd+1, i), :);
                E = U(:, 2:end, :) - U(:, 1, :);

                normal{i} = zeros(nd, obj.NElements);
                for j = 1:obj.NElements
                    nullVec = null(squeeze(E(:, :, j)).');
                    if size(nullVec, 2) > 1
                        nullVec = nullVec(:, 1);
                    end
                    toVertex = V(:, i, j) - U(:, 1, j);
                    if dot(toVertex, nullVec) > 0
                        nullVec = -nullVec;
                    end
                    normal{i}(:, j) = nullVec / norm(nullVec);
                end
            end
        end
    end

    methods (Access = protected)
        function obj = setFaces(obj)
            % SETFACES Generate face connectivity from simplex elements.
            %
            %   obj = setFaces(obj) extracts face definitions from element
            %   connectivity and builds the face-to-element mapping table
            %   for simplicial meshes.

            nd = obj.NDims;
            E = obj.Elements.';
            T = nchoosek(1:nd+1, nd).';
            F = reshape(E(T(:), :), nd, []);
            F = sort(F, 1).';
            [obj.Faces, ~, rows] = unique(F, 'rows', 'stable');
            cols = kron((1:obj.NElements).', ones(nd+1, 1));
            vals = mod(0:(length(rows) - 1), nd+1) + 1;
            obj.FaceToElementTable = sparse(rows, cols, vals, obj.NFaces, obj.NElements);
            obj.Boundary = find(accumarray(rows, 1) == 1);
        end

        function obj = forcePositiveOrientation(obj)
            % FORCEPOSITIVEORIENTATION Ensure positive element
            % orientations.
            %
            %   obj = forcePositiveOrientation(obj) reorders element
            %   vertices to ensure positive Jacobian determinants, which is
            %   essential for proper integration and geometric
            %   computations.

            detEJac = obj.computeElementJacobianDeterminants();
            EI = find(detEJac < 0);

            if ~isempty(EI)
                switch obj.NDims
                    case 1
                        obj.Elements(EI, :) = obj.Elements(EI, [2, 1]);
                    case 2
                        obj.Elements(EI, :) = obj.Elements(EI, [1, 3, 2]);

                    case 3
                        obj.Elements(EI, :) = obj.Elements(EI, [1, 3, 2, 4]);

                    otherwise
                        nVerts = obj.NVerticesPerElement;
                        obj.Elements(EI, :) = obj.Elements(EI, [1:nVerts - 2, nVerts, nVerts - 1]);
                end
                obj.resetTransform();
            end
        end
    end
end

% ========================================================================
% REFINEMENT HELPER FUNCTIONS
% ========================================================================

function [V, E] = refineLines(V0, E0)
% REFINELINES Refine 1D line segments by midpoint subdivision.
%
%   [V, E] = refineLines(V0, E0) subdivides each line segment into
%   two segments by adding midpoint vertices.

%< Start with original vertices
V = V0;

%< Get vertex coordinates for line endpoints
V1 = V(E0(:, 1), :); % Start points [nElements × nDims]
V2 = V(E0(:, 2), :); % End points [nElements × nDims]

%< Compute midpoints for each line segment
midpoints = (V1 + V2) / 2; % [nElements × nDims]

%< Append midpoints to vertex list
V = [V; midpoints];

%< Indices of new midpoint vertices
midpointIndices = size(V0, 1) + (1:size(E0, 1))';

%< Create refined elements: each original line becomes two lines
E = zeros(2*size(E0, 1), 2);
E(1:2:end, :) = [E0(:, 1), midpointIndices]; % First half: start to midpoint
E(2:2:end, :) = [midpointIndices, E0(:, 2)]; % Second half: midpoint to end
end

function [V, E] = refineTriangles(V0, E0)
% REFINETRIANGLES Refine triangular elements by edge subdivision.
%
%   [V, E] = refineTriangles(V0, E0) subdivides each triangle into
%   four triangles by adding midpoint vertices on all edges.
%
% Notes:
%   Each triangle is subdivided into 4 sub-triangles:
%   - 3 corner triangles (one at each original vertex)
%   - 1 central triangle (connecting the 3 edge midpoints)

nV = size(V0, 1);
nE = size(E0, 1);

%< Extract all triangle edges (each triangle contributes 3 edges)
allEdges = [E0(:, [1, 2]); E0(:, [2, 3]); E0(:, [3, 1])]; % [3×nE × 2]

%< Sort edge vertices for canonical representation
sortedEdges = sort(allEdges, 2);

%< Find unique edges and their midpoints
[uniqueEdges, ~, edgeID] = unique(sortedEdges, 'rows');
nUniqueEdges = size(uniqueEdges, 1);

%< Compute edge midpoints
midpoints = (V0(uniqueEdges(:, 1), :) + V0(uniqueEdges(:, 2), :)) / 2;

%< Append midpoints to vertex list
V = [V0; midpoints];

%< Map edge IDs to midpoint vertex indices
midpointIndices = nV + (1:nUniqueEdges).';
edgeMidpoints = midpointIndices(edgeID);

%< Extract midpoint indices for each triangle's three edges
I12 = edgeMidpoints(1:nE); % Midpoint of edge v1-v2
I23 = edgeMidpoints(nE+1:2*nE); % Midpoint of edge v2-v3
I31 = edgeMidpoints(2*nE+1:end); % Midpoint of edge v3-v1

%< Original triangle vertex indices
v1 = E0(:, 1);
v2 = E0(:, 2);
v3 = E0(:, 3);

%< Create 4 sub-triangles per original triangle
E = zeros(4*nE, 3);
E(4*(0:nE - 1)+1, :) = [v1, I12, I31]; % Corner triangle at v1
E(4*(0:nE - 1)+2, :) = [v2, I23, I12]; % Corner triangle at v2
E(4*(0:nE - 1)+3, :) = [v3, I31, I23]; % Corner triangle at v3
E(4*(0:nE - 1)+4, :) = [I12, I23, I31]; % Central triangle
end

function [V, E] = refineTetrahedra(V0, E0)
% REFINETETRAHEDRA Refine tetrahedral elements by edge subdivision.
%
%   [V, E] = refineTetrahedra(V0, E0) subdivides each tetrahedron into
%   eight tetrahedra by adding midpoint vertices on all edges.
%
% Notes:
%   Each tetrahedron has 6 edges, creating 6 new midpoint vertices.
%   The 8 sub-tetrahedra consist of:
%   - 4 corner tetrahedra (one at each original vertex)
%   - 4 interior tetrahedra (connecting midpoints)

nV = size(V0, 1);
nE = size(E0, 1);

%< Define the 6 edges of a tetrahedron (vertex index pairs)
edgePairs = [1, 2; 1, 3; 1, 4; 2, 3; 2, 4; 3, 4];

%< Extract all edges from all tetrahedra
allEdges = reshape(E0(:, edgePairs), [], 2); % [6×nE × 2]

%< Sort edge vertices for canonical representation
sortedEdges = sort(allEdges, 2);

%< Find unique edges and compute their midpoints
[uniqueEdges, ~, edgeID] = unique(sortedEdges, 'rows');
midpoints = (V0(uniqueEdges(:, 1), :) + V0(uniqueEdges(:, 2), :)) / 2;

%< Append midpoints to vertex list
V = [V0; midpoints];

%< Map edges to midpoint indices for each tetrahedron
midpointIndices = nV + (1:size(uniqueEdges, 1));
edgeMidpoints = reshape(midpointIndices(edgeID), nE, 6); % [nE × 6]

%< Extract midpoint indices for each tetrahedron's 6 edges
m12 = edgeMidpoints(:, 1); % Midpoint of edge v1-v2
m13 = edgeMidpoints(:, 2); % Midpoint of edge v1-v3
m14 = edgeMidpoints(:, 3); % Midpoint of edge v1-v4
m23 = edgeMidpoints(:, 4); % Midpoint of edge v2-v3
m24 = edgeMidpoints(:, 5); % Midpoint of edge v2-v4
m34 = edgeMidpoints(:, 6); % Midpoint of edge v3-v4

%< Original tetrahedron vertex indices
v1 = E0(:, 1);
v2 = E0(:, 2);
v3 = E0(:, 3);
v4 = E0(:, 4);

%< Create 8 sub-tetrahedra per original tetrahedron
E = zeros(8*nE, 4);

%< Corner tetrahedra (one at each original vertex)
E(8*(0:nE - 1)+1, :) = [v1, m12, m13, m14]; % Corner at v1
E(8*(0:nE - 1)+2, :) = [v2, m12, m23, m24]; % Corner at v2
E(8*(0:nE - 1)+3, :) = [v3, m13, m23, m34]; % Corner at v3
E(8*(0:nE - 1)+4, :) = [v4, m14, m24, m34]; % Corner at v4

%< Interior tetrahedra (connecting midpoints)
E(8*(0:nE - 1)+5, :) = [m12, m13, m23, m14]; % Interior tetrahedron 1
E(8*(0:nE - 1)+6, :) = [m12, m23, m24, m14]; % Interior tetrahedron 2
E(8*(0:nE - 1)+7, :) = [m13, m23, m34, m14]; % Interior tetrahedron 3
E(8*(0:nE - 1)+8, :) = [m23, m24, m34, m14]; % Interior tetrahedron 4
end