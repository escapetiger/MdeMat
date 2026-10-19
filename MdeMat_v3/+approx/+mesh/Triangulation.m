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
    % Examples:
    %   % Create 2D triangular mesh
    %   vertices2D = [0, 0; 1, 0; 0.5, 1; 1, 1];
    %   elements2D = [1, 2, 3; 2, 4, 3];
    %   triMesh = Triangulation(vertices2D, elements2D);
    %
    %   % Create 3D tetrahedral mesh
    %   vertices3D = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];
    %   elements3D = [1, 2, 3, 4];
    %   tetMesh = Triangulation(vertices3D, elements3D);
    %
    %   % Map reference coordinates to physical coordinates
    %   xRef = [0.2; 0.3];  % Reference coordinates in triangle
    %   xPhys = triMesh.collocate(xRef, 1);  % Map to element 1
    %
    %   % Perform uniform mesh refinement
    %   refinedMesh = triMesh.refine(2);  % 2 levels: 4^2 = 16× more triangles
    %
    %   % Compute geometric quantities
    %   h = triMesh.computeMeasure();           % Minimum edge length
    %   detJ = triMesh.computeElementJacobianDeterminants([1, 2]);
    %   graph = triMesh.graphify();             % Extract connectivity graph
    %
    %   % Access mesh topology
    %   nVerts = triMesh.nVertices;       % Number of mesh vertices
    %   nElems = triMesh.nElements;       % Number of elements
    %   V = triMesh.findElementVertices([1, 3]);  % Vertices for elements 1,3
    %
    % Notes:
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
            %
            % Inputs:
            %   obj - The Triangulation object
            %
            % Outputs:
            %   graphObj - Graph object with mesh vertices and edges

            E = obj.elements.';
            T = nchoosek(1:obj.nVerticesPerElement, 2).';
            edges = reshape(E(T(:), :), 2, []);
            edges = unique(sort(edges, 1).', 'rows', 'stable');
            graphObj = approx.mesh.Graph(obj.vertices, edges);
        end

        function Y = collocate(obj, X, IE)
            % COLLOCATE Map reference coordinates to physical coordinates.
            %
            %   Y = collocate(obj, X, IE) transforms reference Cartesian
            %   coordinates to physical Cartesian coordinates using affine
            %   transformations defined by simplex element vertices.
            %
            % Inputs:
            %   obj - The Triangulation object
            %   X - Reference coordinates (nDims × nPoints)
            %   IE - Element linear indices (nElements × 1)
            %
            % Outputs:
            %   Y - Physical coordinates (nDims × nPoints × nElements)

            if nargin < 3 || isempty(IE)
                IE = 1:obj.nElements;
            end
            d = obj.nDims;
            V = obj.findElementVertices(IE);
            X = reshape(X, d, []);
            V = reshape(V, d, d+1, []);
            b = V(:, 1, :);
            A = V(:, 2:end, :) - b;
            A = reshape(A, d, d, 1, []);
            X = reshape(X, 1, d, [], 1);
            b = reshape(b, d, 1, 1, []);
            Y = squeeze(sum(A .* X, 2) + b);
        end

        function newObj = refine(obj, nLevels)
            % REFINE Perform uniform subdivision of simplicial mesh.
            %
            %   newObj = refine(obj, nLevels) creates a refined mesh by
            %   uniformly subdividing each simplex element. The number of
            %   refinement levels determines the subdivision density.
            %
            % Inputs:
            %   obj - The Triangulation object
            %   nLevels - Number of refinement levels (non-negative integer)
            %
            % Outputs:
            %   newObj - Refined Triangulation object

            if nLevels == 0
                newObj = obj;
                return;
            end

            V = obj.vertices;
            E = obj.elements;
            for i = 1:nLevels
                switch obj.nDims
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
            %
            % Inputs:
            %   obj - The Triangulation object
            %
            % Outputs:
            %   h - Minimum edge length (positive scalar)

            d = obj.nDims;
            T = nchoosek(1:(d + 1), 2).';
            E = obj.elements(:, T(:)).';
            V = reshape(obj.vertices(E(:), :), 2, [], d);
            edgeLengths = sqrt(sum(squeeze(diff(V, 1, 1)).^2, 2));
            h = min(edgeLengths, [], 'all');
        end

        function EJ = computeAllElementJacobians(obj)
            % COMPUTEALLELEMENTJACOBIANS Compute Jacobian matrices for all
            % elements.
            %
            %   EJ = computeAllElementJacobians(obj) computes the Jacobian
            %   matrix for coordinate transformation from reference to
            %   physical coordinates for every element in the mesh.
            %
            % Inputs:
            %   obj - The Triangulation object
            %
            % Outputs:
            %   EJ - Jacobian matrices (nDims × nDims × nElements)

            V = obj.findElementVertices();
            EJ = permute(V(:, 2:end, :)-V(:, 1, :), [2, 1, 3]);
        end

        function detEJ = computeAllElementJacobianDeterminants(obj)
            % COMPUTEALLELEMENTJACOBIANDETERMINANTS Compute Jacobian
            % determinants for all elements.
            %
            %   detEJ = computeAllElementJacobianDeterminants(obj) computes
            %   the determinant of the Jacobian matrix for each simplex
            %   element, used in integration and volume calculations.
            %
            % Inputs:
            %   obj - The Triangulation object
            %
            % Outputs:
            %   detEJ - Jacobian determinants (nElements × 1)

            EJ = obj.computeElementJacobians();
            switch obj.nDims
                case 1
                    detEJ = EJ;
                case 2
                    EJ11 = squeeze(EJ(1, 1, :));
                    EJ12 = squeeze(EJ(1, 2, :));
                    EJ21 = squeeze(EJ(2, 1, :));
                    EJ22 = squeeze(EJ(2, 2, :));
                    detEJ = EJ11 .* EJ22 - EJ12 .* EJ21;
                case 3
                    EJ11 = squeeze(EJ(1, 1, :));
                    EJ12 = squeeze(EJ(1, 2, :));
                    EJ13 = squeeze(EJ(1, 3, :));
                    EJ21 = squeeze(EJ(2, 1, :));
                    EJ22 = squeeze(EJ(2, 2, :));
                    EJ23 = squeeze(EJ(2, 3, :));
                    EJ31 = squeeze(EJ(3, 1, :));
                    EJ32 = squeeze(EJ(3, 2, :));
                    EJ33 = squeeze(EJ(3, 3, :));
                    detEJ = EJ11 .* (EJ22 .* EJ33 - EJ23 .* EJ32) ...
                        -EJ12 .* (EJ21 .* EJ33 - EJ23 .* EJ31) ...
                        +EJ13 .* (EJ21 .* EJ32 - EJ22 .* EJ31);

                otherwise
                    detEJ = arrayfun(@(i) det(EJ(:, :, i)), 1:size(EJ, 3));
            end
        end

        function invEJ = computeAllElementInverseJacobians(obj)
            % COMPUTEALLELEMENTINVERSEJACOBIANS Compute inverse Jacobian
            % matrices for all elements.
            %
            %   invEJ = computeAllElementInverseJacobians(obj) computes the
            %   inverse of the Jacobian matrix for each element, used in
            %   gradient computations and coordinate transformations.
            %
            % Inputs:
            %   obj - The Triangulation object
            %
            % Outputs:
            %   invEJ - Inverse Jacobian matrices (nDims × nDims × nElements)

            switch obj.nDims
                case 1
                    detEJ = obj.computeElementJacobianDeterminants();
                    invEJ = 1 ./ detEJ;

                case 2
                    EJ = obj.computeElementJacobians();
                    detEJ = obj.computeElementJacobianDeterminants();

                    EJ11 = squeeze(EJ(1, 1, :));
                    EJ12 = squeeze(EJ(1, 2, :));
                    EJ21 = squeeze(EJ(2, 1, :));
                    EJ22 = squeeze(EJ(2, 2, :));

                    invEJ = zeros(size(EJ));
                    invEJ(1, 1, :) = EJ22 ./ detEJ;
                    invEJ(1, 2, :) = -EJ12 ./ detEJ;
                    invEJ(2, 1, :) = -EJ21 ./ detEJ;
                    invEJ(2, 2, :) = EJ11 ./ detEJ;

                case 3
                    EJ = obj.computeElementJacobians();
                    detEJ = obj.computeElementJacobianDeterminants();

                    EJ11 = squeeze(EJ(1, 1, :));
                    EJ12 = squeeze(EJ(1, 2, :));
                    EJ13 = squeeze(EJ(1, 3, :));
                    EJ21 = squeeze(EJ(2, 1, :));
                    EJ22 = squeeze(EJ(2, 2, :));
                    EJ23 = squeeze(EJ(2, 3, :));
                    EJ31 = squeeze(EJ(3, 1, :));
                    EJ32 = squeeze(EJ(3, 2, :));
                    EJ33 = squeeze(EJ(3, 3, :));

                    invEJ = zeros(size(EJ));
                    invEJ(1, 1, :) = (EJ22 .* EJ33 - EJ23 .* EJ32) ./ detEJ;
                    invEJ(1, 2, :) = (EJ13 .* EJ32 - EJ12 .* EJ33) ./ detEJ;
                    invEJ(1, 3, :) = (EJ12 .* EJ23 - EJ13 .* EJ22) ./ detEJ;
                    invEJ(2, 1, :) = (EJ23 .* EJ31 - EJ21 .* EJ33) ./ detEJ;
                    invEJ(2, 2, :) = (EJ11 .* EJ33 - EJ13 .* EJ31) ./ detEJ;
                    invEJ(2, 3, :) = (EJ13 .* EJ21 - EJ11 .* EJ23) ./ detEJ;
                    invEJ(3, 1, :) = (EJ21 .* EJ32 - EJ22 .* EJ31) ./ detEJ;
                    invEJ(3, 2, :) = (EJ12 .* EJ31 - EJ11 .* EJ32) ./ detEJ;
                    invEJ(3, 3, :) = (EJ11 .* EJ22 - EJ12 .* EJ21) ./ detEJ;

                otherwise
                    EJ = obj.computeElementJacobians();
                    invEJ = zeros(size(EJ));
                    for i = 1:size(EJ, 3)
                        invEJ(:, :, i) = inv(EJ(:, :, i));
                    end
            end
        end

        function detFJ = computeAllFaceJacobianDeterminants(obj)
            % COMPUTEALLFACEJACOBIANDETERMINANTS Compute face Jacobian
            % determinants for all faces.
            %
            %   detFJ = computeAllFaceJacobianDeterminants(obj) computes
            %   Jacobian determinants for coordinate transformation on all
            %   element faces, used in boundary and interface integration.
            %
            % Inputs:
            %   obj - The Triangulation object
            %
            % Outputs:
            %   detFJ - Face Jacobian determinants (cell array by local face index)

            d = obj.nDims;
            detFJ = cell(1, d+1);

            for i = 1:d + 1
                V = obj.findElementVertices();
                U = V(:, setdiff(1:d+1, i), :);
                G = U(:, 2:end, :) - U(:, 1, :);
                G = pagemtimes(permute(G, [2, 1, 3]), G);
                gramDets = arrayfun(@(k) det(G(:, :, k)), 1:obj.nElements);
                detFJ{i} = reshape(sqrt(gramDets), [], 1);
            end
        end

        function outN = computeAllOutwardNormals(obj)
            % COMPUTEALLOUTWARDNORMALS Compute outward normal vectors for
            % all faces.
            %
            %   outN = computeAllOutwardNormals(obj) computes unit outward
            %   normal vectors for all element faces, essential for
            %   boundary conditions and flux computations.
            %
            % Inputs:
            %   obj - The Triangulation object
            %
            % Outputs:
            %   outN - Outward normal vectors (cell array by local face index)

            d = obj.nDims;
            outN = cell(1, d+1);

            for i = 1:d + 1
                V = obj.findElementVertices();
                U = V(:, setdiff(1:d+1, i), :);
                E = U(:, 2:end, :) - U(:, 1, :);

                outN{i} = zeros(d, obj.nElements);
                for j = 1:obj.nElements
                    nullVec = null(squeeze(E(:, :, j)).');
                    if size(nullVec, 2) > 1
                        nullVec = nullVec(:, 1);
                    end
                    toVertex = V(:, i, j) - U(:, 1, j);
                    if dot(toVertex, nullVec) > 0
                        nullVec = -nullVec;
                    end
                    outN{i}(:, j) = nullVec / norm(nullVec);
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
            %
            % Inputs:
            %   obj - The Triangulation object
            %
            % Outputs:
            %   obj - The Triangulation object with face connectivity

            d = obj.nDims;
            E = obj.elements.';
            T = nchoosek(1:d+1, d).';
            F = reshape(E(T(:), :), d, []);
            F = sort(F, 1).';
            [obj.faces, ~, rows] = unique(F, 'rows', 'stable');
            cols = kron((1:obj.nElements).', ones(d+1, 1));
            vals = mod(0:(length(rows) - 1), d+1) + 1;
            obj.faceToElementTable = sparse(rows, cols, vals, obj.nFaces, obj.nElements);
            obj.boundary = find(accumarray(rows, 1) == 1);
        end

        function obj = forcePositiveOrientation(obj)
            % FORCEPOSITIVEORIENTATION Ensure positive element
            % orientations.
            %
            %   obj = forcePositiveOrientation(obj) reorders element
            %   vertices to ensure positive Jacobian determinants, which is
            %   essential for proper integration and geometric
            %   computations.
            %
            % Inputs:
            %   obj - The Triangulation object
            %
            % Outputs:
            %   obj - The Triangulation object with corrected orientations

            detEJ = obj.computeElementJacobianDeterminants();
            IE = find(detEJ < 0);

            if ~isempty(IE)
                switch obj.nDims
                    case 1
                        obj.elements(IE, :) = obj.elements(IE, [2, 1]);
                    case 2
                        obj.elements(IE, :) = obj.elements(IE, [1, 3, 2]);

                    case 3
                        obj.elements(IE, :) = obj.elements(IE, [1, 3, 2, 4]);

                    otherwise
                        nVerts = obj.nVerticesPerElement;
                        obj.elements(IE, :) = obj.elements(IE, [1:nVerts - 2, nVerts, nVerts - 1]);
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
%
% Inputs:
%   V0 - Original vertices (nVertices × nDims)
%   E0 - Original elements (nElements × 2)
%
% Outputs:
%   V - Refined vertices including midpoints (nVertices + nElements × nDims)
%   E - Refined elements (2×nElements × 2)

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
% Inputs:
%   V0 - Original vertices (nVertices × 2)
%   E0 - Original triangular elements (nElements × 3)
%
% Outputs:
%   V - Refined vertices including edge midpoints
%   E - Refined triangular elements (4×nElements × 3)
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
% Inputs:
%   V0 - Original vertices (nVertices × 3)
%   E0 - Original tetrahedral elements (nElements × 4)
%
% Outputs:
%   V - Refined vertices including edge midpoints
%   E - Refined tetrahedral elements (8×nElements × 4)
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