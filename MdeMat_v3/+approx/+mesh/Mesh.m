classdef Mesh < handle
    % MESH Base class for computational meshes.
    %
    %   Mesh provides foundational functionality for computational meshes
    %   in numerical approximation methods. It defines common interfaces
    %   for coordinate transformations, mesh refinement, and geometric
    %   computations while supporting efficient caching of expensive
    %   operations.
    %
    %   The class manages coordinate transformations between reference and
    %   physical coordinate systems, computes Jacobian matrices and their
    %   determinants for integration quadrature, and provides mesh
    %   refinement capabilities for adaptive algorithms.
    %
    %   Transformation data is cached using bitwise flags to avoid
    %   redundant computations:
    %   - 0b00001: Outward normals computed
    %   - 0b00010: Face Jacobian determinants computed
    %   - 0b00100: Element inverse Jacobians computed
    %   - 0b01000: Element Jacobian determinants computed
    %   - 0b10000: Element Jacobians computed
    %
    % See also:
    %   approx.mesh.Grid, approx.mesh.PolytopalMesh

    properties (Constant)
        EVAL_EJ = 0b100000
        EVAL_DETEJ = 0b010000
        EVAL_INVEJ = 0b001000
        EVAL_DETFJ = 0b000100
        EVAL_OUTN = 0b000010
        EVAL_BETA = 0b000001
    end

    properties (Access = public)
        vertices % Vertex coordinates (nVertices × nDims or cell array)
        elements % Element connectivity (nElements × nVerticesPerElement or cell array)
        transform % Transformation data cache (struct)
        direction % Direction vector (nDims x 1)
        status % Status mask
    end

    properties (Dependent)
        nDims % Number of spatial dimensions (integer)
        nVertices % Number of vertices in the mesh (integer)
        nElements % Number of elements in the mesh (integer)
        nVerticesPerElement % Number of vertices per element (integer)
    end

    methods
        function obj = Mesh(vertices, elements)
            % MESH Constructor for Mesh base class.
            %
            %   obj = Mesh(vertices, elements) creates a mesh from vertex
            %   coordinates and element connectivity.
            %
            % Inputs:
            %   vertices - Vertex coordinate matrix or cell array
            %   elements - Element connectivity matrix or cell array
            %
            % Outputs:
            %   obj - Constructed Mesh object

            obj.vertices = vertices;
            obj.elements = elements;
            obj.resetTransform();
        end

        function n = get.nDims(obj)
            % GET.NDIMS Get the number of spatial dimensions.

            if iscell(obj.vertices)
                n = length(obj.vertices);
            else
                n = size(obj.vertices, 2);
            end
        end

        function n = get.nVertices(obj)
            % GET.NVERTICES Get the number of vertices in the mesh.

            if iscell(obj.vertices)
                n = prod(cellfun(@length, obj.vertices));
            else
                n = size(obj.vertices, 1);
            end
        end

        function n = get.nElements(obj)
            % GET.NELEMENTS Get the number of elements in the mesh.

            if iscell(obj.elements)
                n = prod(cellfun(@(x) length(x), obj.elements));
            else
                n = size(obj.elements, 1);
            end
        end

        function n = get.nVerticesPerElement(obj)
            % GET.NVERTICESPERELEMENT Get the number of vertices per element.

            if iscell(obj.elements)
                n = 2^obj.nDims;
            else
                n = size(obj.elements, 2);
            end
        end

        function obj = resetTransform(obj)
            % RESETTRANSFORM Reset transformation data cache.
            %
            %   obj = resetTransform(obj) clears all cached transformation
            %   data and resets the computation flags. This method should
            %   be called whenever the mesh geometry changes.
            %
            % Inputs:
            %   obj - The Mesh object
            %
            % Outputs:
            %   obj - The Mesh object with cleared cache

            obj.transform = struct( ...
                'EJ', [], ... % Element Jacobian matrices
                'detEJ', [], ... % Element Jacobian determinants
                'invEJ', [], ... % Element inverse Jacobian matrices
                'detFJ', [], ... % Face Jacobian determinants
                'outN', []);
            obj.status = 0b000000; % Status mask
        end

        function EJ = computeElementJacobians(obj, IE)
            % COMPUTEELEMENTJACOBIANS Compute element Jacobian matrices.
            %
            %   EJ = computeElementJacobians(obj) computes Jacobian
            %   matrices for coordinate transformation from reference to
            %   physical coordinates for all elements.
            %
            %   EJ = computeElementJacobians(obj, IE) computes Jacobian
            %   matrices only for the specified elements.
            %
            % Inputs:
            %   obj - The Mesh object
            %   IE - Element linear indices (optional, vector)
            %
            % Outputs:
            %   EJ - Jacobian matrices (nDims × nDims × nElements)
            %
            % Notes:
            %   The Jacobian matrix $J_{ij} = \frac{\partial x_i}{\partial \xi_j}$
            %   where $x$ are physical coordinates and $\xi$ are reference
            %   coordinates. For coordinate transformations: $dx = J d\xi$.

            if ~bitand(obj.status, obj.EVAL_EJ)
                obj.transform.EJ = obj.computeAllElementJacobians();
                obj.status = bitor(obj.status, obj.EVAL_EJ);
            end

            if nargin < 2 || isempty(IE)
                EJ = obj.transform.EJ;
            else
                EJ = obj.transform.EJ(:, :, IE);
            end
        end

        function detEJ = computeElementJacobianDeterminants(obj, IE)
            % COMPUTEELEMENTJACOBIANDETERMINANTS Compute element Jacobian
            % determinants.
            %
            %   detEJ = computeElementJacobianDeterminants(obj) computes
            %   Jacobian determinants for coordinate transformation from
            %   reference to physical coordinates for all elements.
            %
            %   detEJ = computeElementJacobianDeterminants(obj, IE)
            %   computes Jacobian determinants only for the specified
            %   elements.
            %
            % Inputs:
            %   obj - The Mesh object
            %   IE - Element linear indices (optional, vector)
            %
            % Outputs:
            %   detEJ - Element Jacobian determinants (vector)
            %
            % Notes:
            %   The Jacobian determinant appears in integration formulas:
            %   $\int_{\Omega} f(x) dx = \int_{\hat{\Omega}} f(F(\xi)) |det(J)| d\xi$
            %   where $F$ is the coordinate transformation and $\hat{\Omega}$
            %   is the reference element.

            if ~bitand(obj.status, obj.EVAL_DETEJ)
                obj.transform.detEJ = obj.computeAllElementJacobianDeterminants();
                obj.status = bitor(obj.status, obj.EVAL_DETEJ);
            end

            if nargin < 2 || isempty(IE) || isscalar(obj.transform.detEJ)
                detEJ = obj.transform.detEJ;
            else
                detEJ = obj.transform.detEJ(IE);
            end
        end

        function invEJ = computeElementInverseJacobians(obj, IE)
            % COMPUTEELEMENTINVERSEJACOBIANS Compute element inverse
            % Jacobian matrices.
            %
            %   invEJ = computeElementInverseJacobians(obj) computes
            %   inverse Jacobian matrices for coordinate transformation
            %   from physical to reference coordinates for all elements.
            %
            %   invEJ = computeElementInverseJacobians(obj, IE) computes
            %   inverse Jacobian matrices only for the specified elements.
            %
            % Inputs:
            %   obj - The Mesh object
            %   IE - Element linear indices (optional, vector)
            %
            % Outputs:
            %   invEJ - Inverse Jacobian matrices (nDims × nDims × nElements)
            %
            % Notes:
            %   The inverse Jacobian matrix appears in gradient computations:
            %   $\nabla_x f = J^{-T} \nabla_\xi f$
            %   where $\nabla_x$ and $\nabla_\xi$ are gradients in physical
            %   and reference coordinates respectively.

            if ~bitand(obj.status, obj.EVAL_INVEJ)
                obj.transform.invEJ = obj.computeAllElementInverseJacobians();
                obj.status = bitor(obj.status, obj.EVAL_INVEJ);
            end

            if nargin < 2 || isempty(IE) || ismatrix(obj.transform.invEJ)
                invEJ = obj.transform.invEJ;
            else
                invEJ = obj.transform.invEJ(:, :, IE);
            end
        end

        function detFJ = computeFaceJacobianDeterminants(obj, IE, IFE)
            % COMPUTEFACEJACOBIANDETERMINANTS Compute face Jacobian
            % determinants.
            %
            %   detFJ = computeFaceJacobianDeterminants(obj, IE, IFE)
            %   computes Jacobian determinants for coordinate transformation
            %   on element faces for integration over boundaries.
            %
            % Inputs:
            %   obj - The Mesh object
            %   IE - Element linear indices (vector or [])
            %   IFE - Local face index (scalar, 1 to nFaces)
            %
            % Outputs:
            %   detFJ - Face Jacobian determinants (vector)
            %
            % Notes:
            %   Face Jacobian determinants are used in boundary integration:
            %   $\int_{\partial\Omega} f(x) dS = \int_{\partial\hat{\Omega}} f(F(\xi)) |det(J_F)| d\xi$
            %   where $\partial\Omega$ is the boundary and $J_F$ is the face Jacobian.

            if ~bitand(obj.status, obj.EVAL_DETFJ)
                obj.transform.detFJ = obj.computeAllFaceJacobianDeterminants();
                obj.status = bitor(obj.status, obj.EVAL_DETFJ);
            end

            if isempty(IE) || isscalar(obj.transform.detFJ{IFE})
                detFJ = obj.transform.detFJ{IFE};
            else
                detFJ = obj.transform.detFJ{IFE}(IE);
            end
        end

        function outN = computeOutwardNormals(obj, IE, IFE)
            % COMPUTEOUTWARDNORMALS Compute outward normal vectors.
            %
            %   outN = computeOutwardNormals(obj, IE, IFE) computes unit
            %   outward normal vectors on element faces for boundary
            %   conditions and flux computations.
            %
            % Inputs:
            %   obj - The Mesh object
            %   IE - Element linear indices (vector or [])
            %   IFE - Local face index (scalar, 1 to nFaces)
            %
            % Outputs:
            %   outN - Outward unit normal vectors (nDims × nElements)
            %
            % Notes:
            %   Outward normals are used in boundary conditions and flux
            %   computations. The direction is defined such that the normal
            %   points away from the element interior.

            if ~bitand(obj.status, obj.EVAL_OUTN)
                obj.transform.outN = obj.computeAllOutwardNormals();
                obj.status = bitor(obj.status, obj.EVAL_OUTN);
            end

            if isempty(IE) || size(obj.transform.outN{IFE}, 2) == 1
                outN = obj.transform.outN{IFE};
            else
                outN = obj.transform.outN{IFE}(:, IE);
            end
        end

        function v = computeDirection(obj)
            % COMPUTEDIRECTION Compute a direction vector.
            %
            %   v = computeDirection(obj) computes a unit direction vector
            %   to define the "plus" and "minus" side of a face.
            %
            % Inputs:
            %   obj - The Mesh object
            %
            % Outputs:
            %   v - Direction vector (nDims × 1)
            %
            % Notes:
            %   If nDims = 1, v = 1.
            %   If nDims > 1, v is a unit vector that is not parallel
            %   with any normal vectors of element boundaries.

            if ~bitand(obj.status, obj.EVAL_BETA)
                obj.direction = obj.computeUniqueDirection();
                obj.status = bitor(obj.status, obj.EVAL_BETA);
            end

            v = obj.direction;
        end

        function ratio = computeFaceScalingFactors(obj, IE, IFE)
            % COMPUTEFACESCALINGFACTORS Compute face scaling factors.
            %
            %   ratio = computeFaceScalingFactors(obj, IE, IFE) computes
            %   face scaling factors for coordinate transformation on
            %   element faces for integration over boundaries.
            %
            % Inputs:
            %   obj - The Mesh object
            %   IE - Element linear indices (vector or [])
            %   IFE - Local face index (scalar, 1 to nFaces)
            %
            % Outputs:
            %   ratio - Face scaling factors (vector)

            detEJ = obj.computeElementJacobianDeterminants(IE);
            detFJ = obj.computeFaceJacobianDeterminants(IE, IFE);
            ratio = detFJ ./ detEJ;
        end
    end

    methods (Abstract)
        % GRAPHIFY Convert mesh to graph representation.
        graphObj = graphify(obj)

        % COLLOCATE Map reference coordinates to physical coordinates.
        Y = collocate(obj, X, IE)

        % REFINE Create refined mesh with smaller elements.
        newObj = refine(obj, nLevels)

        % COMPUTEMEASURE Compute characteristic mesh size.
        h = computeMeasure(obj)

        % COMPUTEALLELEMENTJACOBIANS Compute all element Jacobian matrices.
        EJ = computeAllElementJacobians(obj)

        % COMPUTEALLELEMENTJACOBIANDETERMINANTS Compute all element
        % Jacobian determinants.
        detEJ = computeAllElementJacobianDeterminants(obj)

        % COMPUTEALLELEMENTINVERSEJACOBIANS Compute all element inverse
        % Jacobian matrices.
        invEJ = computeAllElementInverseJacobians(obj)

        % COMPUTEALLFACEJACOBIANDETERMINANTS Compute all face Jacobian
        % determinants.
        detFJ = computeAllFaceJacobianDeterminants(obj)

        % COMPUTEALLOUTWARDNORMALS Compute all outward normal vectors.
        outN = computeAllOutwardNormals(obj)

        % COMPUTEUNIQUEDIRECTION Compute a unique direction vector.
        v = computeUniqueDirection(obj)
    end
end