classdef (Abstract) Mesh < handle
    % MESH Base class for computational meshes with vertices, faces, and elements.
    %
    %   Mesh provides a foundational representation for computational
    %   meshes consisting of vertices, faces, and elements. It serves as
    %   the base class for specific mesh types (SimplexMesh, etc.) and
    %   integrates seamlessly with the Graph-based mesh infrastructure.
    %
    %   In an n-dimensional mesh:
    %   - n-elements are the primary n-dimensional geometric entities
    %   - (n-1)-faces are the (n-1)-dimensional boundaries between elements
    %   - vertices are 0-dimensional geometric entities
    %   - edges are 1-dimensional connections between vertices
    %
    %   For example:
    %   - In 1D: elements are line segments, faces are points
    %   - In 2D: elements are triangles/quads, faces are edges
    %   - In 3D: elements are tetrahedra/hexahedra, faces are triangles/quads
    %
    % Notes:
    %   This is an abstract class that cannot be instantiated directly.
    %   Concrete subclasses must implement the abstract methods for
    %   generating mesh topology and coordinate transformations.
    %
    % See also:
    %   approx.mesh.Graph, approx.mesh.SimplexMesh, approx.mesh.Grid

    properties
        vertices % Matrix of vertex coordinates (nVertices × nDims)
        elements % Matrix of element definitions (nElements × nVerticesPerElement)
        faces % Matrix of face definitions (nFaces × nVerticesPerFace)
        boundary % Vector of boundary face indices (nBoundaryFaces × 1)
        faceToElementTable % Connectivity table between faces and elements
        vertexToVertexTable % Connectivity table between vertices
        faceToFaceTable % Connectivity table between faces
    end

    properties (Dependent)
        nDims % Number of spatial dimensions (integer)
        nVertices % Number of vertices in the mesh (integer)
        nElements % Number of elements in the mesh (integer)
        nVerticesPerElement % Number of vertices per element (integer)
        nFaces % Number of faces in the mesh (integer)
        nVerticesPerFace % Number of vertices per face (integer)
        nBoundaryFaces % Number of boundary faces (integer)
    end

    methods
        function obj = Mesh(vertices, elements)
            % MESH Constructor for Mesh.
            %
            %   obj = Mesh(vertices, elements) creates a mesh from vertex
            %   coordinates and element connectivity, automatically
            %   generating the graph structure and face information.
            %
            % Inputs:
            %   vertices - Matrix where each row is vertex coordinates
            %   elements - Matrix where each row is vertex linear indices
            %
            % Outputs:
            %   obj - Constructed Mesh object

            obj.vertices = vertices;
            obj.elements = elements;
            obj.forcePositiveOrientation();
            obj.setFaces();
            obj.vertexToVertexTable = [];
        end

        function n = get.nDims(obj)
            % GET.NDIMS Get the number of spatial dimensions.

            n = size(obj.vertices, 2);
        end

        function n = get.nVertices(obj)
            % GET.NVERTICES Get the number of vertices in the mesh.

            n = size(obj.vertices, 1);
        end

        function n = get.nElements(obj)
            % GET.NELEMENTS Get the number of elements in the mesh.

            n = size(obj.elements, 1);
        end

        function n = get.nVerticesPerElement(obj)
            % GET.NVERTICESPERELEMENT Get the number of vertices per element.

            n = size(obj.elements, 2);
        end

        function n = get.nFaces(obj)
            % GET.NFACES Get the number of faces in the mesh.

            n = size(obj.faces, 1);
        end

        function n = get.nVerticesPerFace(obj)
            % GET.NVERTICESPERFACE Get the number of vertices per face.

            n = size(obj.faces, 2);
        end

        function n = get.nBoundaryFaces(obj)
            % GET.NBOUNDARYFACES Get the number of boundary faces.

            n = length(obj.boundary);
        end

        function V = findElementVertices(obj, L)
            % FINDELEMENTVERTICES Find element vertices using linear
            % indices.
            %
            %   V = findElementVertices(obj, L) finds vertices using
            %   element linear indices.
            %
            % Inputs:
            %   obj - The SimplexMesh object
            %   L - Element linear indices (ne x 1 vector or [])
            %
            % Outputs:
            %   Y - Physical Cartesian coordinates (nd × nvp x ne)

            if nargin < 2 || isempty(L)
                L = 1:obj.nElements;
            end

            E = obj.elements(L, :).';
            V = obj.vertices(E(:), :);
            V = reshape(V, obj.nVerticesPerElement, [], obj.nDims);
            V = permute(V, [3, 1, 2]);
        end

        function obj = setPeriodic(obj, shift, tol)
            % SETPERIODIC Set periodic mapping within the mesh.
            %
            %   obj = setPeriodic(obj, shift) builds vertex-to-vertex table
            %   and face-to-face table to represent the periodic mapping
            %   between mesh elements.
            %
            % Inputs:
            %   shift - Periodic shifts (m x d matrix)
            %   tol   - Tolerance for floating-point comparisons
            %
            % Output:
            %   obj - The Mesh object

            if nargin < 3 || isempty(tol), tol = 1e-12; end

            obj.setVertexToVertexTable(shift, tol);
            obj.setFaceToFaceTable();
        end

        function L = findInteriorElements(obj)
            % FINDINTERIORELEMENTS Find linear indices of interior
            % elements.
            %
            %   L = findInteriorElements(obj) computes the linear indices
            %   of interior elements.
            %
            % Inputs:
            %   obj - The Grid object
            %
            % Outputs:
            %   L - Linear indices of interior elements (vector)

            mask = false(obj.nFaces, 1);
            mask(obj.boundary) = true;
            L = find(~any(obj.faceToElementTable(mask, :), 1));
        end

        function L = findBoundaryElements(obj, i)
            % FINDBOUNDARYELEMENTS Find linear indices of boundary
            % elements.
            %
            %   L = findBoundaryElements(obj) computes the linear
            %   indices of elements sharing faces with boundary.
            %
            % Inputs:
            %   obj - The Mesh object
            %
            % Outputs:
            %   L - Linear indices of boundary elements (vector)

            T = nchoosek(1:obj.nVerticesPerElement, obj.nVerticesPerFace);
            F = sort(obj.elements(:, T(i, :)), 2);
            [~, I] = ismember(F, obj.faces, 'rows');
            L = find(ismember(I, obj.boundary));
        end

        function L = findNeighborElements(obj, i, K, bc)
            % FINDNEIGHBORELEMENTS Find neighbors through i-th face.
            %
            %   L = findNeighborElements(obj, i, K, bc) finds elements that
            %   share the i-th face with elements in set K.
            %
            % Inputs:
            %   obj - The Mesh object
            %   i - Face index within element (1 to nVerticesPerElement)
            %   K - Set of element indices
            %   bc - Boundary condition parameter
            %
            % Outputs:
            %   L - Neighbor elements sharing i-th face

            if nargin < 4 || isempty(bc)
                bc = 'strict';
            end

            if isempty(K)
                L = [];
                return;
            end

            T = nchoosek(1:obj.nVerticesPerElement, obj.nVerticesPerFace);
            F = sort(obj.elements(K, T(i, :)), 2);
            [~, I] = ismember(F, obj.faces, 'rows');

            switch lower(bc)
                case 'periodic'
                    % Extend faceIdx by mapping through periodic pairs
                    FP = any(obj.faceToFaceTable(I, :), 2);
                    J = find(obj.faceToFaceTable(I(FP), :));
                    I = unique([I; J(:)]);
                case 'strict'
                    % Remove boundary faces
                    I = setdiff(I, obj.boundary);
            end

            % Use faceToElementTable to find adjacent elements
            A = obj.faceToElementTable(I, :);
            A(:, K) = 0; % Remove self from neighbor set
            L = find(any(A, 1));
        end
    end

    methods (Abstract)
        % GRAPHIFY Build the graph.
        graphObj = graphify(obj)

        % COLLOCATE Map reference nodes to physical nodes.
        Y = collocate(obj, X, L)

        % REFINE Refine the mesh.
        newObj = refine(obj, nLevels)

        % COMPUTEMEASURE Compute the mesh measure.
        h = computeMeasure(obj)

        % COMPUTEELEMENTJACOBIANDETERMINANTS Compute element Jacobian
        % determinants.
        detJ = computeElementJacobianDeterminants(obj)

        % COMPUTEELEMENTINVERSEJACOBIANS Compute element inverse Jacobian
        % matrices.
        invJ = computeElementInverseJacobians(obj)

        % COMPUTEFACEJACOBIANDETERMINANTS Compute face Jacobian
        % determinants.
        detJFace = computeFaceJacobianDeterminants(obj, faceIndex)

        % COMPUTEOUTWARDNORMALS Compute outward normal vectors.
        normals = computeOutwardNormals(obj, faceIndex)

    end

    methods (Abstract, Access = protected)
        % SETFACES Generate faces from element connectivity.
        obj = setFaces(obj)

        % FORCEPOSITIVEORIENTATION Force positive orientation.
        obj = forcePositiveOrientation(obj)
    end

    methods (Access = private)
        function obj = setVertexToVertexTable(obj, shift, tol)
            % SETVERTEXTOVERTEXTABLE Builds vertex-to-vertex table.
            %
            %   obj = setVertexToVertexTable(obj, shift, tol) builds
            %   vertex-to-vertex table.
            %
            % Inputs:
            %   obj - The Mesh object
            %   shift - Periodic shifts (m x d matrix)
            %   tol   - Tolerance for floating-point comparisons
            %
            % Output:
            %   obj - The Mesh object

            nV = obj.nVertices;
            V = obj.vertices;
            T = sparse(nV, nV);
            for k = 1:size(shift, 1)
                D = pdist2(V, V-shift(k, :));
                [I, J] = find(D < tol);
                T = T + sparse(I, J, 1, nV, nV);
                T = T + sparse(J, I, 1, nV, nV);
            end
            obj.vertexToVertexTable = spones(T);
        end

        function obj = setFaceToFaceTable(obj)
            % SETFACETOFACETABLE Builds face-to-face table.
            %
            %   obj =  setFaceToFaecTable(obj) builds face-to-face table.
            %
            % Inputs:
            %   obj - The Mesh object
            %
            % Output:
            %   obj - The Mesh object

            V2V = obj.vertexToVertexTable; %< nVertices x nVertices sparse

            %< Map face vertices through vertexToVertexTable
            %< faces: (T x nvpf) indices of vertices
            %< For each vertex, get the mapped vertex index (1 per vertex)

            %< Convert sparse matrix to cell array of mapped vertices for fast indexing
            %< Because vertexToVertexTable is binary, for each vertex there should be at most one mapped vertex per periodic image.

            %< Step 1: Extract mapping indices for each vertex in each face
            mappedFaces = zeros(size(obj.faces));
            for k = 1:obj.nVerticesPerFace
                rowIdx = obj.faces(:, k);
                %< For each vertex v = rowIdx(i), find the mapped vertex from V2V(v,:)
                %< Use find with 'first' for each row:
                %< Unfortunately MATLAB doesn't have a direct vectorized way to get indices of nonzeros in each row.
                %< But we can use a trick with sparse:

                %< Get columns of nonzero entries for rows in rowIdx
                [r, c] = find(V2V(rowIdx, :));
                %< r are row subscripts relative to rowIdx (which is vector)
                %< So mapped vertex for faces(i,k) is c(r == i)
                %< Build an array for mapped vertex for each face vertex:
                idxs = accumarray(r, c, [length(rowIdx), 1], @(x) x(1), 0);
                %< accumarray with function @(x) x(1) picks first mapping if multiple exist
                mappedFaces(:, k) = idxs;
            end

            %< If any mappedFaces are zero (no mapping), remove those faces
            validMask = all(mappedFaces > 0, 2);
            mappedFaces = mappedFaces(validMask, :);
            origFaceIdx = find(validMask);

            %< Sort vertices in mappedFaces for canonical representation
            sortedMappedFaces = sort(mappedFaces, 2);

            %< Sort original faces too (for lookup)
            sortedFaces = sort(obj.faces, 2);

            %< Now find matching mapped faces in original faces
            %< Create a combined matrix for matching:
            %< For each mapped face, find which row in sortedFaces matches exactly

            %< Use ismember with 'rows' (vectorized)
            [~, loc] = ismember(sortedMappedFaces, sortedFaces, 'rows');

            %< Build sparse connection matrix: face i connects to face loc(i)
            %< Note: only keep connections where loc > 0
            validPairs = loc > 0;
            iFaces = origFaceIdx(validPairs);
            jFaces = loc(validPairs);

            %< Build symmetric sparse matrix
            faceToFace = sparse([iFaces; jFaces], [jFaces; iFaces], 1, obj.nFaces, obj.nFaces);

            obj.faceToFaceTable = spones(faceToFace);
        end
    end
end