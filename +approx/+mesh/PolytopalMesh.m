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
        Faces % Face connectivity (NFaces × NVerticesPerFace)
        Boundary % Boundary face indices (NBoundaryFaces × 1)
        FaceToElementTable % Face-element connectivity (sparse NFaces × NElements)
        VertexToVertexTables % Periodic vertex mappings (cell of sparse NVertices × NVertices)
        FaceToFaceTable % Periodic face mapping (sparse NFaces × NFaces)
    end

    properties (Dependent)
        NFaces % Number of faces in the mesh (integer)
        NVerticesPerFace % Number of vertices per face (integer)
        NBoundaryFaces % Number of boundary faces (integer)
    end

    methods
        function obj = PolytopalMesh(vertices, elements)
            % POLYTOPALMESH Constructor for PolytopalMesh.
            %
            %   obj = PolytopalMesh(vertices, elements) creates a mesh from
            %   vertex coordinates and element connectivity, automatically
            %   generating face topology and ensuring proper element
            %   orientation.

            arguments
                vertices {mustBeNumeric} = []
                elements {mustBeNumeric} = []
            end

            obj@approx.mesh.Mesh(vertices, elements);
            obj.forcePositiveOrientation();
            obj.setFaces();
            obj.VertexToVertexTables = [];
            obj.FaceToFaceTable = [];
        end

        function n = get.NFaces(obj)
            % GET.NFACES Get the number of faces in the mesh.

            n = size(obj.Faces, 1);
        end

        function n = get.NVerticesPerFace(obj)
            % GET.NVERTICESPERFACE Get the number of vertices per face.

            n = size(obj.Faces, 2);
        end

        function n = get.NBoundaryFaces(obj)
            % GET.NBOUNDARYFACES Get the number of boundary faces.

            n = length(obj.Boundary);
        end

        function V = findElementVertices(obj, EI)
            % FINDELEMENTVERTICES Extract vertex coordinates for specified
            % elements.
            %
            %   V = findElementVertices(obj) returns vertex coordinates for
            %   all elements in the mesh.
            %
            %   V = findElementVertices(obj, EI) returns vertex coordinates
            %   only for the specified elements.

            arguments
                obj approx.mesh.PolytopalMesh
                EI {mustBeInteger} = []
            end

            if isempty(EI)
                EI = 1:obj.NElements;
            end
            E = obj.Elements(EI, :).';
            V = obj.Vertices(E(:), :);
            V = reshape(V, obj.NVerticesPerElement, [], obj.NDims);
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

            arguments
                obj approx.mesh.PolytopalMesh
                bbox {mustBeNumeric}
                tol {mustBeNumeric} = 1e-12
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

        function IEI = getInternalElements(obj)
            % GETINTERIORELEMENTS Get linear indices of internal elements.
            %
            %   IEI = getInternalElements(obj) computes the linear indices
            %   of elements that do not touch the mesh boundary.

            arguments   
                obj approx.mesh.PolytopalMesh
            end

            mask = false(obj.NFaces, 1);
            mask(obj.Boundary) = true;
            IEI = find(~any(obj.FaceToElementTable(mask, :), 1));
        end

        function BEI = getBoundaryElements(obj, LFI)
            % GETBOUNDARYELEMENTS Get linear indices of boundary elements.
            %
            %   BEI = getBoundaryElements(obj) computes the linear indices
            %   of all elements that have at least one face on the
            %   boundary.
            %
            %   BEI = getBoundaryElements(obj, LFI) computes the linear
            %   indices of elements whose LFI-th face is on the boundary.

            arguments   
                obj approx.mesh.PolytopalMesh
                LFI {mustBeInteger, mustBePositive} = []
            end

            if isempty(LFI)
                BEI = find(any(obj.FaceToElementTable(obj.Boundary, :), 1));
                return;
            end

            E = obj.Elements.';
            T = nchoosek(1:obj.NVerticesPerElement, obj.NVerticesPerFace).';
            F = reshape(E(T(:, LFI), :), obj.NVerticesPerFace, []);
            F = sort(F, 1).';
            [~, I] = ismember(F, obj.Faces, 'rows');
            BEI = find(ismember(I, obj.Boundary));
        end

        function IEI = findInternalElements(obj, EI, LFI)
            % FINDINTERIORELEMENTS Find internal elements across faces.
            %
            %   IEI = findInternalElements(obj, EI, LFI) finds elements
            %   whose LFI-th face is not on the boundary within the
            %   specified elements @a EI.

            arguments   
                obj approx.mesh.PolytopalMesh
                EI {mustBeInteger} = []
                LFI {mustBeInteger, mustBePositive} = 1
            end

            E = obj.Elements(EI, :).';
            T = nchoosek(1:obj.NVerticesPerElement, obj.NVerticesPerFace).';
            F = reshape(E(T(:, LFI), :), obj.NVerticesPerFace, []);
            F = sort(F, 1).';
            [~, I] = ismember(F, obj.Faces, 'rows');
            IEI = find(~ismember(I, obj.Boundary));
        end

        function BEI = findBoundaryElements(obj, EI, LFI)
            % FINDBOUNDARYELEMENTS Find boundary elements across faces.
            %
            %   BEI = findBoundaryElements(obj, EI, LFI) finds elements
            %   whose LFI-th face is on the boundary within the specified
            %   elements @a EI.

            arguments   
                obj approx.mesh.PolytopalMesh
                EI {mustBeInteger} = []
                LFI {mustBeInteger, mustBePositive} = 1
            end

            E = obj.Elements(EI, :).';
            T = nchoosek(1:obj.NVerticesPerElement, obj.NVerticesPerFace).';
            F = reshape(E(T(:, LFI), :), obj.NVerticesPerFace, []);
            F = sort(F, 1).';
            [~, I] = ismember(F, obj.Faces, 'rows');
            BEI = find(ismember(I, obj.Boundary));
        end

        function [NEI, NLFI] = findNeighborElements(obj, EI, LFI, BC)
            % FINDNEIGHBORELEMENTS Find neighboring elements across faces.
            %
            %   [NEI, NLFI] = findNeighborElements(obj, EI, LFI, BC) computes
            %   the linear indices of neighbor elements sharing the LFI-th
            %   face of prime elements EI, with the specified boundary
            %   condition treatment.

            arguments   
                obj approx.mesh.PolytopalMesh
                EI {mustBeInteger} = []
                LFI {mustBeInteger, mustBePositive} = 1
                BC {mustBeMember(BC, {'', 'periodic', 'dirichlet'})} = ''
            end

            E = obj.Elements.';
            T = nchoosek(1:obj.NVerticesPerElement, obj.NVerticesPerFace).';
            F = E(T(:, LFI), EI);
            F = sort(F, 1).';
            [~, I] = ismember(F, obj.Faces, 'rows');

            switch lower(BC)
                case 'periodic'
                    P = any(obj.FaceToFaceTable(I, :), 2);
                    I(P) = find(obj.FaceToFaceTable(I(P), :), 1, 'first');
                case 'dirichlet'
                    I = setdiff(I, obj.Boundary);
            end

            R = obj.FaceToElementTable(I, :);
            R(:, EI) = 0;
            [NLFI, NEI] = max(R, [], 2);
            NEI(NLFI == 0) = 0;
        end
    
        function v = computeUniqueDirection(obj)
            % COMPUTEUNIQUEDIRECTION Compute a unique direction vector.
            %
            %   v = computeUniqueDirection(obj) computes a unique unit
            %   direction vector to define the "plus" and "minus" side of a
            %   face.

            arguments   
                obj approx.mesh.Mesh
            end

            nd = obj.NDims;
            if nd == 1
                v = 1;
                return;
            end

            outNormal = obj.computeAllOutwardNormals();
            
            while true
                v = rand(nd,1) + 0.1; 
                u = v / norm(v);
            
                dots = abs(outNormal.' * u);  
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

            n = obj.NVertices;
            V = obj.Vertices;

            obj.VertexToVertexTables = cell(1, size(shift, 1));
            for k = 1:size(shift, 1)
                D = pdist2(V, V-shift(k, :));
                [I, J] = find(D < tol);
                Tk = sparse(I, J, 1, n, n);
                obj.VertexToVertexTables{k} = spones(Tk);
            end
        end

        function obj = setFaceToFaceTable(obj)
            % SETFACETOFACETABLES Build face-to-face mapping table.
            
            F = obj.Faces; % nF × m (m = NVerticesPerFace)
            nF = obj.NFaces;
            m = obj.NVerticesPerFace;

            %< Canonical representation of existing faces
            sortedFaces = sort(F, 2);

            faceToFace = sparse(nF, nF);

            %< Loop only over small number of shifts (8 in 2D, 26 in 3D)
            for s = 1:length(obj.VertexToVertexTables)
                Ts = obj.VertexToVertexTables{s};

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

            obj.FaceToFaceTable = spones(faceToFace);
        end
    end
end