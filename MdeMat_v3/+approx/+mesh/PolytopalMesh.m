classdef PolytopalMesh < approx.mesh.Mesh
    % POLYTOPALMESH Base class for polytopal meshes.
    %
    %   PolytopalMesh provides a foundational representation for polytopal
    %   meshes consisting of vertices, faces, and elements with explicit
    %   connectivity information. It serves as the base class for
    %   unstructured mesh types and integrates seamlessly with the
    %   Graph-based mesh infrastructure for visualization and analysis.
    %
    %   The mesh topology follows standard finite element conventions:
    %   - In 1D: elements are line segments, faces are points (vertices)
    %   - In 2D: elements are polygons, faces are edges
    %   - In 3D: elements are polyhedra, faces are polygons
    %   - In nD: elements are n-polytopes, faces are (n-1)-polytopes
    %
    %   Connectivity is stored using explicit arrays and sparse tables:
    %   - vertices: Coordinate matrix (nVertices × nDims)
    %   - elements: Connectivity matrix (nElements × nVerticesPerElement)
    %   - faces: Connectivity matrix (nFaces × nVerticesPerFace)
    %   - faceToElementTable: Sparse connectivity (nFaces × nElements)
    %
    %   The class supports periodic boundary conditions through
    %   vertex-to-vertex and face-to-face mapping tables that identify
    %   geometrically equivalent entities across periodic boundaries.
    %
    % See also:
    %   approx.mesh.Mesh, approx.mesh.Graph, approx.mesh.PolytopalMeshPlotter

    properties (Access = public)
        faces % Face connectivity (nFaces × nVerticesPerFace)
        boundary % Boundary face indices (nBoundaryFaces × 1)
        faceToElementTable % Face-element connectivity (sparse nFaces × nElements)
        vertexToVertexTables % Periodic vertex mappings (cell of sparse nVertices × nVertices)
        faceToFaceTable % Periodic face mapping (sparse nFaces × nFaces)
    end

    properties (Dependent)
        nFaces % Number of faces in the mesh (integer)
        nVerticesPerFace % Number of vertices per face (integer)
        nBoundaryFaces % Number of boundary faces (integer)
    end

    methods
        function obj = PolytopalMesh(vertices, elements)
            % POLYTOPALMESH Constructor for PolytopalMesh.
            %
            %   obj = PolytopalMesh(vertices, elements) creates a mesh from
            %   vertex coordinates and element connectivity, automatically
            %   generating face topology and ensuring proper element
            %   orientation.
            %
            % Inputs:
            %   vertices - Vertex coordinate matrix (nVertices × nDims)
            %   elements - Element connectivity matrix (nElements × nVerticesPerElement)
            %
            % Outputs:
            %   obj - Constructed PolytopalMesh object

            obj@approx.mesh.Mesh(vertices, elements);
            obj.forcePositiveOrientation();
            obj.setFaces();
            obj.vertexToVertexTables = [];
            obj.faceToFaceTable = [];
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

        function V = findElementVertices(obj, IE)
            % FINDELEMENTVERTICES Extract vertex coordinates for specified
            % elements.
            %
            %   V = findElementVertices(obj) returns vertex coordinates for
            %   all elements in the mesh.
            %
            %   V = findElementVertices(obj, IE) returns vertex coordinates
            %   only for the specified elements.
            %
            % Inputs:
            %   obj - The PolytopalMesh object
            %   IE - Element linear indices (optional, vector)
            %
            % Outputs:
            %   V - Element vertex coordinates (nDims × nVerticesPerElement × nElements)

            if nargin < 2 || isempty(IE)
                IE = 1:obj.nElements;
            end
            E = obj.elements(IE, :).';
            V = obj.vertices(E(:), :);
            V = reshape(V, obj.nVerticesPerElement, [], obj.nDims);
            V = permute(V, [3, 1, 2]);
        end

        function obj = setPeriodic(obj, bbox, tol)
            % SETPERIODIC Establish periodic boundary mappings.
            %
            %   obj = setPeriodic(obj, bbox) builds vertex-to-vertex and
            %   face-to-face connectivity tables to represent periodic
            %   mappings between geometrically equivalent mesh entities.
            %
            %   obj = setPeriodic(obj, bbox, tol) uses the specified
            %   tolerance for identifying equivalent vertices.
            %
            % Inputs:
            %   obj - The PolytopalMesh object
            %   shift - Periodic shift vectors (nShifts × nDims)
            %   tol - Numerical tolerance(optional, default: 1e-12)
            %
            % Outputs:
            %   obj - The PolytopalMesh object

            if nargin < 3 || isempty(tol)
                tol = 1e-12;
            end

            bbox = reshape(bbox, 2, []);
            nDims = size(bbox, 2);
            shift = cell(1, nDims);
            [shift{:}] = ndgrid(-1:1);
            shift = cellfun(@(x) x(:), shift, 'Un', 0);
            shift = [shift{:}];
            shift = shift(any(shift ~= 0, 2), :);
            shift = shift .* diff(bbox, 1, 1);

            obj.setVertexToVertexTables(shift, tol);
            obj.setFaceToFaceTable();
        end

        function IE = getInternalElements(obj)
            % GETINTERIORELEMENTS Get linear indices of internal elements.
            %
            %   IE = getInternalElements(obj) computes the linear indices
            %   of elements that do not touch the mesh boundary.
            %
            % Inputs:
            %   obj - The PolytopalMesh object
            %
            % Outputs:
            %   IE - Linear indices of internal elements (vector)

            mask = false(obj.nFaces, 1);
            mask(obj.boundary) = true;
            IE = find(~any(obj.faceToElementTable(mask, :), 1));
        end

        function IE = getBoundaryElements(obj, IFE)
            % GETBOUNDARYELEMENTS Get linear indices of boundary elements.
            %
            %   IE = getBoundaryElements(obj) computes the linear indices
            %   of all elements that have at least one face on the
            %   boundary.
            %
            %   IE = getBoundaryElements(obj, IFE) computes the linear
            %   indices of elements whose IFE-th face is on the boundary.
            %
            % Inputs:
            %   obj - The PolytopalMesh object
            %   IFE - Local face index (optional, scalar)
            %
            % Outputs:
            %   IE - Linear indices of boundary elements (vector)

            if nargin < 2 || isempty(IFE)
                IE = find(any(obj.faceToElementTable(obj.boundary, :), 1));
                return;
            end

            E = obj.elements.';
            T = nchoosek(1:obj.nVerticesPerElement, obj.nVerticesPerFace).';
            F = reshape(E(T(:, IFE), :), obj.nVerticesPerFace, []);
            F = sort(F, 1).';
            [~, I] = ismember(F, obj.faces, 'rows');
            IE = find(ismember(I, obj.boundary));
        end

        function JE = findInternalElements(obj, IE, IFE)
            % FINDINTERIORELEMENTS Find internal elements across faces.
            %
            %   JE = findInternalElements(obj, IE, IFE) finds elements
            %   whose IFE-th face is not on the boundary within the
            %   specified elements @a IE.
            %
            % Inputs:
            %   obj - The Grid object
            %   IE - Linear indices of primary elements (vector)
            %   IFE - Local face index (1 to 2×nDims)
            %
            % Outputs:
            %   JE - Linear indices of internal elements (vector)

            E = obj.elements(IE, :).';
            T = nchoosek(1:obj.nVerticesPerElement, obj.nVerticesPerFace).';
            F = reshape(E(T(:, IFE), :), obj.nVerticesPerFace, []);
            F = sort(F, 1).';
            [~, I] = ismember(F, obj.faces, 'rows');
            JE = find(~ismember(I, obj.boundary));
        end

        function JE = findBoundaryElements(obj, IE, IFE)
            % FINDBOUNDARYELEMENTS Find boundary elements across faces.
            %
            %   JE = findBoundaryElements(obj, IE, IFE) finds elements
            %   whose IFE-th face is on the boundary within the specified
            %   elements @a IE.
            %
            % Inputs:
            %   obj - The Grid object
            %   IE - Linear indices of primary elements (vector)
            %   IFE - Local face index (1 to 2×nDims)
            %
            % Outputs:
            %   JE - Linear indices of boundary elements (vector)

            E = obj.elements(IE, :).';
            T = nchoosek(1:obj.nVerticesPerElement, obj.nVerticesPerFace).';
            F = reshape(E(T(:, IFE), :), obj.nVerticesPerFace, []);
            F = sort(F, 1).';
            [~, I] = ismember(F, obj.faces, 'rows');
            JE = find(ismember(I, obj.boundary));
        end

        function [JE, JFE] = findNeighborElements(obj, IE, IFE, BC)
            % FINDNEIGHBORELEMENTS Find neighboring elements across faces.
            %
            %   [JE, JFE] = findNeighborElements(obj, IE, IFE, BC) computes
            %   the linear indices of neighbor elements sharing the IFE-th
            %   face of prime elements IE, with the specified boundary
            %   condition treatment.
            %
            % Inputs:
            %   obj - The PolytopalMesh object
            %   IE - Linear indices of prime elements (vector)
            %   IFE - Local face index of prime elements (scalar)
            %   BC - Boundary condition type ('periodic' or 'dirichlet')
            %
            % Outputs:
            %   JE - Linear indices of neighbor elements (vector)
            %   JFE - Local face indices of neighbor elements (vector)

            E = obj.elements.';
            T = nchoosek(1:obj.nVerticesPerElement, obj.nVerticesPerFace).';
            F = E(T(:, IFE), IE);
            F = sort(F, 1).';
            [~, I] = ismember(F, obj.faces, 'rows');

            switch lower(BC)
                case 'periodic'
                    P = any(obj.faceToFaceTable(I, :), 2);
                    I(P) = find(obj.faceToFaceTable(I(P), :), 1, 'first');
                case 'dirichlet'
                    I = setdiff(I, obj.boundary);
            end

            R = obj.faceToElementTable(I, :);
            R(:, IE) = 0;
            [JFE, JE] = max(R, [], 2);
            JE(JFE == 0) = 0;
        end
    
        function v = computeUniqueDirection(obj)
            % COMPUTEUNIQUEDIRECTION Compute a unique direction vector.
            %
            %   v = computeUniqueDirection(obj) computes a unique unit
            %   direction vector to define the "plus" and "minus" side of a
            %   face.
            %
            % Inputs:
            %   obj - The PolytopalMesh object
            %
            % Outputs:
            %   v - Direction vector (nDims × 1)

            d = obj.nDims;
            if d == 1
                v = 1;
                return;
            end

            outN = obj.computeAllOuterNormals();
            
            while true
                v = rand(d,1) + 0.1; 
                u = v / norm(v);
            
                dots = abs(outN.' * u);  
                if all(dots < 1 - 1e-12)
                    return;
                end
            end
        end
    end

    methods (Abstract, Access = protected)
        % SETFACES Generate face connectivity from element topology.
        obj = setFaces(obj)

        % FORCEPOSITIVEORIENTATION Ensure positive element orientations.
        obj = forcePositiveOrientation(obj)
    end

    methods (Access = private)
        function obj = setVertexToVertexTables(obj, shift, tol)
            % SETVERTEXTOVERTEXTABLES Build vertex-to-vertex mapping table.
            %
            %   obj = setVertexToVertexTables(obj, shift, tol) constructs a
            %   sparse matrix identifying vertices that are geometrically
            %   equivalent under the specified periodic shifts.
            %
            % Inputs:
            %   obj - The PolytopalMesh object
            %   shift - Periodic shift vectors (nShifts × nDims)
            %   tol - Geometric tolerance for vertex matching
            %
            % Outputs:
            %   obj - The PolytopalMesh object with vertex mapping table

            n = obj.nVertices;
            V = obj.vertices;

            obj.vertexToVertexTables = cell(1, size(shift, 1));
            for k = 1:size(shift, 1)
                D = pdist2(V, V-shift(k, :));
                [I, J] = find(D < tol);
                Tk = sparse(I, J, 1, n, n);
                obj.vertexToVertexTables{k} = spones(Tk);
            end
        end

        function obj = setFaceToFaceTable(obj)
            F = obj.faces; % nF × m (m = nVerticesPerFace)
            nF = obj.nFaces;
            m = obj.nVerticesPerFace;

            %< Canonical representation of existing faces
            sortedFaces = sort(F, 2);

            faceToFace = sparse(nF, nF);

            %< Loop only over small number of shifts (8 in 2D, 26 in 3D)
            for s = 1:length(obj.vertexToVertexTables)
                Ts = obj.vertexToVertexTables{s};

                %< Map every face's k-th vertex via this shift, in bulk
                mappedFaces = zeros(nF, m);
                for k = 1:m
                    Mk = Ts(F(:, k), :);
                    [val, idx] = max(Mk, [], 2);
                    idx(~val) = 0;
                    mappedFaces(:, k) = idx;
                end

                %< Keep only faces whose all vertices map under this shift
                valid = all(mappedFaces > 0, 2);
                if ~any(valid), continue; end

                %< Find matches in the original face set
                sortedMapped = sort(mappedFaces(valid, :), 2);
                [isMatch, loc] = ismember(sortedMapped, sortedFaces, 'rows');

                if any(isMatch)
                    iFaces = find(valid);
                    iFaces = iFaces(isMatch);
                    jFaces = loc(isMatch);

                    %< Add symmetric connections
                    T = sparse(iFaces, jFaces, 1, nF, nF);
                    faceToFace = faceToFace + T + T.';
                end
            end

            obj.faceToFaceTable = spones(faceToFace);
        end
    end
end