classdef SimplexMesh < approx.mesh.Mesh
    % SIMPLEXMESH Mesh class for triangular and tetrahedral elements.
    %
    %   SimplexMesh provides a concrete implementation of the Mesh class
    %   for simplicial elements (triangles in 2D, tetrahedra in 3D, line
    %   segments in 1D). It leverages the core.geometry.Simplex class for
    %   geometric operations and coordinate transformations.
    %
    %   The class handles mesh topology generation, coordinate transformations
    %   between reference and physical coordinates, and uniform mesh refinement.
    %   Reference coordinates use standard Cartesian coordinates on the
    %   reference simplex.
    %
    % Examples:
    %   % Create 2D triangular mesh
    %   vertices2D = [0, 0; 1, 0; 0.5, 1; 1, 1];
    %   elements2D = [1, 2, 3; 2, 4, 3];
    %   triMesh = SimplexMesh(vertices2D, elements2D);
    %
    %   % Create 3D tetrahedral mesh
    %   vertices3D = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];
    %   elements3D = [1, 2, 3, 4];
    %   tetMesh = SimplexMesh(vertices3D, elements3D);
    %
    %   % Map reference Cartesian coordinates to physical coordinates
    %   xRef = [0.2, 0.3; 0.1, 0.4];  % Reference Cartesian coordinates
    %   xPhys = triMesh.collocate(xRef, 1);  % Map to element 1
    %
    %   % Refine the mesh uniformly
    %   refinedMesh = triMesh.refine(2);
    %
    % Notes:
    %   - 1D: Line segments (2 vertices per element)
    %   - 2D: Triangles (3 vertices per element)
    %   - 3D: Tetrahedra (4 vertices per element)
    %
    %   Reference coordinates use Cartesian coordinates on the reference
    %   simplex domain.
    %
    % See also:
    %   approx.mesh.Mesh, approx.mesh.Graph, core.geometry.Simplex

    properties (Access = private)
        refElement % Reference simplex geometry (core.geometry.Simplex)
    end

    methods
        function obj = SimplexMesh(vertices, elements)
            % SIMPLEXMESH Constructor for SimplexMesh.
            %
            %   obj = SimplexMesh(vertices, elements) creates a simplex
            %   mesh from vertex coordinates and element connectivity.
            %
            % Inputs:
            %   vertices - Matrix of vertex coordinates (nVertices × nDims)
            %   elements - Matrix of element connectivity (nElements × nVerticesPerElement)
            %
            % Outputs:
            %   obj - Constructed SimplexMesh object

            obj@approx.mesh.Mesh(vertices, elements);
            d = obj.nDims;
            obj.refElement = core.geometry.Simplex([zeros(d, 1), eye(d, d)]);
        end

        function graphObj = graphify(obj)
            % GRAPHIFY Generate graph from simplex elements.
            %
            %   obj = graphify(obj) creates edge connectivity by extracting
            %   all vertex pairs from simplex elements and constructs the
            %   graph object.
            %
            % Inputs:
            %   obj - The SimplexMesh object
            %
            % Outputs:
            %   graphObj - The Graph object

            E = obj.elements.';
            T = nchoosek(1:obj.nVerticesPerElement, 2).';
            edges = reshape(E(T(:), :), 2, []);
            edges = unique(sort(edges, 1).', 'rows', 'stable');
            graphObj = approx.mesh.Graph(obj.vertices, edges);
        end

        function Y = collocate(obj, X, L)
            % COLLOCATE Map reference nodes to physical nodes.
            %
            %   Y = collocate(obj, X, L) maps reference Cartesian
            %   coordinates to physical Cartesian coordinates using affine
            %   transformation for simplex elements.
            %
            % Inputs:
            %   obj - The SimplexMesh object
            %   X - Reference Cartesian coordinates (nDims × nPoints)
            %   L - Element linear indices (nElements x 1)
            %
            % Outputs:
            %   Y - Physical Cartesian coordinates (nDims × nPoints x nElements)

            if nargin < 3 || isempty(L), L = 1:obj.nElements; end

            V = obj.findElementVertices(L);
            Y = obj.refElement.transform(X, V);
        end

        function newObj = refine(obj, nLevels)
            % REFINE Perform uniform refinement of simplex mesh.
            %
            %   newMesh = refine(obj, nLevels) creates a refined mesh
            %   by subdividing each simplex element. The refinement level
            %   determines how many subdivisions are made.
            %
            %   Let k denotes the number of refinement levels. Each
            %   triangle is subdivided into 4^k triangles. Each tetrahedron
            %   is subdivided into 8^k tetrahedra.
            %
            % Inputs:
            %   obj - The SimplexMesh object
            %   nLevels - Number of refinement levels (scalar)
            %
            % Outputs:
            %   newObj - Refined SimplexMesh object

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
            newObj = approx.mesh.SimplexMesh(V, E);
        end

        function h = computeMeasure(obj)
            % COMPUTEMEASURE Compute characteristic mesh size.
            %
            %   h = computeMeshMeasure(obj) returns the minimum
            %   characteristic size among all simplex elements, typically
            %   the minimum edge length or element diameter.
            %
            % Inputs:
            %   obj - The SimplexMesh object
            %
            % Outputs:
            %   h - Mesh measure (scalar)

            d = obj.nDims;
            T = nchoosek(1:(d + 1), 2).';
            E = obj.elements(:, T(:)).';
            V = reshape(obj.vertices(E(:), :), 2, [], d);
            V = sqrt(sum(squeeze(diff(V, 1, 1)).^2, 2));
            h = min(V, [], 'all');
        end

        function detJ = computeElementJacobianDeterminants(obj)
            % COMPUTEELEMENTJACOBIANDETERMINANTS Compute element Jacobian
            % determinants.
            %
            %   detJ = computeElementJacobianDeterminants(obj) computes the
            %   Jacobian determinant for each simplex element based on the
            %   transformation from reference to physical coordinates.
            %
            % Inputs:
            %   obj - The SimplexMesh object
            %
            % Outputs:
            %   detJ - Element Jacobian determinants (vector)

            J = obj.computeElementJacobians();

            switch obj.nDims
                case 1
                    detJ = J;
                case 2
                    J11 = squeeze(J(1, 1, :));
                    J12 = squeeze(J(1, 2, :));
                    J21 = squeeze(J(2, 1, :));
                    J22 = squeeze(J(2, 2, :));
                    detJ = J11 .* J22 - J12 .* J21;
                case 3
                    J11 = squeeze(J(1, 1, :));
                    J12 = squeeze(J(1, 2, :));
                    J13 = squeeze(J(1, 3, :));
                    J21 = squeeze(J(2, 1, :));
                    J22 = squeeze(J(2, 2, :));
                    J23 = squeeze(J(2, 3, :));
                    J31 = squeeze(J(3, 1, :));
                    J32 = squeeze(J(3, 2, :));
                    J33 = squeeze(J(3, 3, :));
                    detJ = J11 .* (J22 .* J33 - J23 .* J32) - J12 .* (J21 .* J33 - J23 .* J31) + J13 .* (J21 .* J32 - J22 .* J31);
                otherwise
                    detJ = arrayfun(@(i) det(J(:, :, i)), 1:size(J, 3));
            end
        end

        function J = computeElementJacobians(obj)
            % COMPUTEELEMENTJACOBIANS Compute element Jacobian matrices.
            %
            %   J = computeElementJacobians(obj) computes the Jacobian
            %   matrix for coordinate transformation. For nonuniform grids,
            %   elements have different diagonal Jacobian matrices.
            %
            % Inputs:
            %   obj - The NonuniformGrid object
            %
            % Outputs:
            %   J - Jacobian matrices (nDims × nDims × nElements)

            V = obj.findElementVertices();
            J = permute(V(:, 2:end, :) - V(:, 1, :), [2, 1, 3]);
        end

        function invJ = computeElementInverseJacobians(obj)
            % COMPUTEELEMENTINVERSEJACOBIANS Compute element inverse
            % Jacobian matrices.
            %
            %   invJ = computeElementInverseJacobians(obj) computes the
            %   inverse Jacobian matrix for coordinate transformation.
            %
            % Inputs:
            %   obj - The SimplexMesh object
            %
            % Outputs:
            %   invJ - Inverse Jacobian matrix (nDims × nDims x nElements)

            switch obj.nDims
                case 1
                    detJ = obj.computeElementJacobianDeterminants();
                    invJ = 1 ./ detJ;

                case 2
                    J = obj.computeElementJacobians();
                    detJ = obj.computeElementJacobianDeterminants();
                    J11 = squeeze(J(1, 1, :));
                    J12 = squeeze(J(1, 2, :));
                    J21 = squeeze(J(2, 1, :));
                    J22 = squeeze(J(2, 2, :));

                    invJ = zeros(size(J));
                    invJ(1, 1, :) = J22 ./ detJ;
                    invJ(1, 2, :) = -J12 ./ detJ;
                    invJ(2, 1, :) = -J21 ./ detJ;
                    invJ(2, 2, :) = J11 ./ detJ;

                case 3
                    J = obj.computeElementJacobians();
                    detJ = obj.computeElementJacobianDeterminants();
                    J11 = squeeze(J(1, 1, :));
                    J12 = squeeze(J(1, 2, :));
                    J13 = squeeze(J(1, 3, :));
                    J21 = squeeze(J(2, 1, :));
                    J22 = squeeze(J(2, 2, :));
                    J23 = squeeze(J(2, 3, :));
                    J31 = squeeze(J(3, 1, :));
                    J32 = squeeze(J(3, 2, :));
                    J33 = squeeze(J(3, 3, :));

                    invJ = zeros(size(J));
                    invJ(1, 1, :) = (J22 .* J33 - J23 .* J32) ./ detJ;
                    invJ(1, 2, :) = (J13 .* J32 - J12 .* J33) ./ detJ;
                    invJ(1, 3, :) = (J12 .* J23 - J13 .* J22) ./ detJ;

                    invJ(2, 1, :) = (J23 .* J31 - J21 .* J33) ./ detJ;
                    invJ(2, 2, :) = (J11 .* J33 - J13 .* J31) ./ detJ;
                    invJ(2, 3, :) = (J13 .* J21 - J11 .* J23) ./ detJ;

                    invJ(3, 1, :) = (J21 .* J32 - J22 .* J31) ./ detJ;
                    invJ(3, 2, :) = (J12 .* J31 - J11 .* J32) ./ detJ;
                    invJ(3, 3, :) = (J11 .* J22 - J12 .* J21) ./ detJ;

                otherwise
                    J = obj.computeElementJacobians();
                    invJ = arrayfun(@(i) inv(J(:, :, i), 1:size(J, 3)));
            end
        end

        function detJFace = computeFaceJacobianDeterminants(obj, i)
            % COMPUTEFACEJACOBIANDETERMINANTS Compute face Jacobian
            % determinants.
            %
            %   detJFace = computeFaceJacobianDeterminants(obj, i)
            %   computes the Jacobian determinant for coordinate
            %   transformation on the specified face.
            %
            % Inputs:
            %   obj - The UniformGrid object
            %   i - Face index (1 to nDims+1)
            %
            % Outputs:
            %   detJFace - Face Jacobian determinants (scalar)

            d = obj.nDims;
            V = obj.findElementVertices();
            U = V(:, setdiff(1:d+1, i), :);
            E = U(:, 2:end, :) - U(:, 1, :);
            G = pagemtimes(permute(E, [2, 1, 3]), E);
            detJFace = sqrt(arrayfun(@(k) det(G(:, :, k)), 1:obj.nElements));
            detJFace = detJFace(:);
        end

        function normals = computeOutwardNormals(obj, i)
            % COMPUTEOUTWARDNORMALS Compute outward normal vectors.
            %
            %   normals = computeOutwardNormals(obj, i) computes
            %   the outward unit normal vector for the specified face.
            %
            % Inputs:
            %   obj - The NonuniformGrid object
            %   i - Face index (1 to nDims+1)
            %
            % Outputs:
            %   normals - Outward unit normal vector (nDims × nElements)

            d = obj.nDims;
            V = obj.findElementVertices();
            U = V(:, setdiff(1:d+1, i), :);
            E = U(:, 2:end, :) - U(:, 1, :);

            normals = zeros(d, obj.nElements);
            for j = 1:size(normals, 2)
                ni = null(squeeze(E(:, :, j)).');
                ni = ni(:, 1);
                if dot(U(:, 1, j)-V(:, i, j), ni) < 0, ni = -ni; end
                normals(:, j) = ni / norm(ni);
            end
        end
    end

    methods (Access = protected)
        function obj = setFaces(obj)
            % SETFACES Generate faces from simplex elements.
            %
            %   obj = setFaces(obj) extracts face connectivity from
            %   elements information.
            %
            % Inputs:
            %   obj - The SimplexMesh object
            %
            % Outputs:
            %   obj - The SimplexMesh object

            d = obj.nDims;
            E = obj.elements.';
            T = nchoosek(1:d+1, d).';
            F = reshape(E(T(:), :), d, []);
            [obj.faces, ~, faceIdx] = unique(sort(F, 1).', 'rows', 'stable');
            obj.boundary = find(accumarray(faceIdx, 1) == 1);

            rows = faceIdx;
            cols = kron((1:obj.nElements)', ones(d+1, 1));
            vals = true(size(rows));
            obj.faceToElementTable = sparse(rows, cols, vals, obj.nFaces, obj.nElements);
        end

        function obj = forcePositiveOrientation(obj)
            % FORCEPOSITIVEORIENTATION Fix element vertex ordering.
            %
            %   obj = forcePositiveOrientation(obj) ensures that all
            %   elements have positive Jacobian determinants using proper
            %   vertex reordering based on element type and dimension.
            %
            % Inputs:
            %   obj - The SimplexMesh object
            %
            % Outputs:
            %   obj - The SimplexMesh object with corrected element orientation

            detJ = obj.computeElementJacobianDeterminants();
            L = find(detJ < 0);
            nCorrected = length(L);

            if nCorrected > 0
                switch obj.nDims
                    case 1
                        %< Line segments: swap endpoints
                        obj.elements(L, :) = obj.elements(L, [2, 1]);

                    case 2
                        %< Triangles: swap last two vertices to reverse orientation
                        obj.elements(L, :) = obj.elements(L, [1, 3, 2]);

                    case 3
                        % Tetrahedra: use standard reordering
                        obj.elements(L, :) = obj.elements(L, [1, 3, 2, 4]);
                end
            end
        end

    end
end

function [V, E] = refineLines(V0, E0)
% REFINELINES Refine 1D line segments.

V = V0;

%< Get vertex positions of endpoints
V1 = V(E0(:, 1), :);  % [nE × d]
V2 = V(E0(:, 2), :);  % [nE × d]

%< Compute midpoints
M = (V1 + V2) / 2;  % [nE × d]

%< Append midpoints to vertex list
V = [V; M];

% Indices of new points
I = size(V0, 1) + (1:size(E0, 1))';

%< Construct new edges
E = zeros(2 * size(E0, 1), 2);
E(1:2:end, :) = [E0(:, 1), I];
E(2:2:end, :) = [I, E0(:, 2)];

end

function [V, E] = refineTriangles(V0, E0)
% REFINETRIANGLES Uniformly refines 2D triangular elements into 4 smaller triangles.
%
% Inputs:
%   V0 - Original vertices [nV × 2]
%   E0 - Original elements (triangles) [nE × 3]
%
% Outputs:
%   V  - New vertices including midpoints [nV + nEdges × 2]
%   E  - New elements [4 × nE × 3]

nV = size(V0, 1);
nE = size(E0, 1);

%< Step 1: Extract and sort all triangle edges
allEdges = [E0(:, [1 2]);
            E0(:, [2 3]);
            E0(:, [3 1])];  % [3*nE × 2]

sortedEdges = sort(allEdges, 2);
[uniqueEdges, ~, edgeID] = unique(sortedEdges, 'rows');  % [nEdges × 2]
nEdges = size(uniqueEdges, 1);

%< Step 2: Compute midpoints of unique edges
midpoints = (V0(uniqueEdges(:, 1), :) + V0(uniqueEdges(:, 2), :)) / 2;

%< Step 3: Append midpoints to vertex list
V = [V0; midpoints];

%< Step 4: Index of each edge's midpoint
I = nV + (1:nEdges).';  % Midpoint indices (1-based)
I12 = I(edgeID(1:nE));           % Edge v1–v2
I23 = I(edgeID(nE+1:2*nE));      % Edge v2–v3
I31 = I(edgeID(2*nE+1:end));     % Edge v3–v1

%< Step 5: Original triangle vertex indices
v1 = E0(:, 1);
v2 = E0(:, 2);
v3 = E0(:, 3);

%< Step 6: Define 4 sub-triangles per triangle
E = zeros(4*nE, 3);
E(4*(0:nE-1)+1, :) = [v1, I12, I31];     % Triangle 1
E(4*(0:nE-1)+2, :) = [v2, I23, I12];     % Triangle 2
E(4*(0:nE-1)+3, :) = [v3, I31, I23];     % Triangle 3
E(4*(0:nE-1)+4, :) = [I12, I23, I31];    % Center triangle

end

function [V, E] = refineTetrahedra(V0, E0)
% REFINETETRAHEDRA Uniform refinement of tetrahedral mesh (8 sub-tets per tet)
%
% Inputs:
%   V0 - [nV × 3] matrix of vertex coordinates
%   E0 - [nE × 4] matrix of tetrahedral elements (indices into V0)
%
% Outputs:
%   V  - New vertex list (including midpoints)
%   E  - New tetrahedral elements (8 × nE × 4)

nV = size(V0, 1);
nE = size(E0, 1);

%< Step 1: Extract all 6 edges from tetrahedra
edgePairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
allEdges = reshape(E0(:, edgePairs), [], 2);   % [6*nE × 2]
sortedEdges = sort(allEdges, 2);

%< Step 2: Unique edges and midpoints
[uniqueEdges, ~, edgeID] = unique(sortedEdges, 'rows');
midpoints = (V0(uniqueEdges(:,1), :) + V0(uniqueEdges(:,2), :)) / 2;

%< Step 3: Append midpoints to vertex list
V = [V0; midpoints];
I = nV + (1:size(uniqueEdges, 1));   % midpoint indices

%< Step 4: Map local edge to midpoint index for each tetrahedron
m = reshape(I(edgeID), nE, 6);  % [nE × 6]
m12 = m(:,1); m13 = m(:,2); m14 = m(:,3);
m23 = m(:,4); m24 = m(:,5); m34 = m(:,6);

v1 = E0(:,1); v2 = E0(:,2); v3 = E0(:,3); v4 = E0(:,4);

%< Step 5: Define the 8 new tetrahedra per original
E = zeros(8*nE, 4);
E(8*(0:nE-1)+1, :) = [v1, m12, m13, m14];
E(8*(0:nE-1)+2, :) = [v2, m12, m23, m24];
E(8*(0:nE-1)+3, :) = [v3, m13, m23, m34];
E(8*(0:nE-1)+4, :) = [v4, m14, m24, m34];
E(8*(0:nE-1)+5, :) = [m12, m13, m23, m14];
E(8*(0:nE-1)+6, :) = [m12, m23, m24, m14];
E(8*(0:nE-1)+7, :) = [m13, m23, m34, m14];
E(8*(0:nE-1)+8, :) = [m23, m24, m34, m14];

end
