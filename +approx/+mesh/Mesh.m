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
        EJacMask = 0b100000
        DetEJacMask = 0b010000
        InvEJacMask = 0b001000
        DetFJacMask = 0b000100
        OutNormalMask = 0b000010
        DirectionMask = 0b000001
    end
    
    properties (Access = public)
        Vertices % Vertex coordinates (NVertices × NDims or cell array)
        Elements % Element connectivity (NElements × NVerticesPerElement or cell array)
        Transform % Transformation data cache (struct)
        Direction % Direction vector (NDims x 1)
        Status % Status mask
    end
    
    properties (Dependent)
        NDims % Number of spatial dimensions (integer)
        NVertices % Number of vertices in the mesh (integer)
        NElements % Number of elements in the mesh (integer)
        NVerticesPerElement % Number of vertices per element (integer)
    end
    
    methods
        function obj = Mesh(vertices, elements)
            % MESH Constructor for Mesh base class.
            %
            %   obj = Mesh(vertices, elements) creates a mesh from vertex
            %   coordinates and element connectivity.
            
            arguments
                vertices {mustBeNumericOrCell} = []
                elements {mustBeNumericOrCell} = []
            end
            
            obj.Vertices = vertices;
            obj.Elements = elements;
            obj.resetTransform();
        end
        
        function n = get.NDims(obj)
            % GET.NDIMS Get the number of spatial dimensions.
            
            if iscell(obj.Vertices)
                n = length(obj.Vertices);
            else
                n = size(obj.Vertices, 2);
            end
        end
        
        function n = get.NVertices(obj)
            % GET.NVERTICES Get the number of vertices in the mesh.
            
            if iscell(obj.Vertices)
                n = prod(cellfun(@length, obj.Vertices));
            else
                n = size(obj.Vertices, 1);
            end
        end
        
        function n = get.NElements(obj)
            % GET.NELEMENTS Get the number of elements in the mesh.
            
            if iscell(obj.Elements)
                n = prod(cellfun(@(x) length(x), obj.Elements));
            else
                n = size(obj.Elements, 1);
            end
        end
        
        function n = get.NVerticesPerElement(obj)
            % GET.NVERTICESPERELEMENT Get the number of vertices per element.
            
            if iscell(obj.Elements)
                n = 2^obj.NDims;
            else
                n = size(obj.Elements, 2);
            end
        end
        
        function obj = resetTransform(obj)
            % RESETTRANSFORM Reset transformation data cache.
            %
            %   obj = resetTransform(obj) clears all cached transformation
            %   data and resets the computation flags. This method should
            %   be called whenever the mesh geometry changes.

            arguments   
                obj approx.mesh.Mesh
            end
            
            obj.Transform = struct( ...
                'EJac', [], ... % Element Jacobian matrices
                'detEJac', [], ... % Element Jacobian determinants
                'invEJac', [], ... % Element inverse Jacobian matrices
                'detFJac', [], ... % Face Jacobian determinants
                'normal', []);
            obj.Status = 0b000000; % Status mask
        end
        
        function EJac = computeElementJacobians(obj, EI)
            % COMPUTEELEMENTJACOBIANS Compute element Jacobian matrices.
            %
            %   EJac = computeElementJacobians(obj) computes Jacobian
            %   matrices for coordinate transformation from reference to
            %   physical coordinates for all elements.
            %
            %   EJac = computeElementJacobians(obj, EI) computes Jacobian
            %   matrices only for the specified elements.
            %
            %   The Jacobian matrix \f$J_{ij} = \frac{\partial
            %   x_i}{\partial \xi_j}\f$ where \f$x\f$ are physical
            %   coordinates and \f$\xi\f$ are reference coordinates. For
            %   coordinate transformations: \f$dx = J d\xi\f$.

            arguments   
                obj approx.mesh.Mesh
                EI {mustBeInteger} = []
            end
            
            if ~bitand(obj.Status, obj.EJacMask)
                obj.Transform.EJac = obj.computeAllElementJacobians();
                obj.Status = bitor(obj.Status, obj.EJacMask);
            end
            
            if isempty(EI)
                EJac = obj.Transform.EJac;
            else
                EJac = obj.Transform.EJac(:, :, EI);
            end
        end
        
        function detEJac = computeElementJacobianDeterminants(obj, EI)
            % COMPUTEELEMENTJACOBIANDETERMINANTS Compute element Jacobian
            % determinants.
            %
            %   detEJac = computeElementJacobianDeterminants(obj) computes
            %   Jacobian determinants for coordinate transformation from
            %   reference to physical coordinates for all elements.
            %
            %   detEJac = computeElementJacobianDeterminants(obj, EI)
            %   computes Jacobian determinants only for the specified
            %   elements.
            %
            %   The Jacobian determinant appears in integration formulas:
            %
            %   \f[
            %       \int_{\Omega} f(x) dx = \int_{\hat{\Omega}} f(F(\xi))
            %       |det(J)| d\xi
            %   \f]
            %
            %   where \f$F\f$ is the coordinate transformation and
            %   \f$\hat{\Omega}\f$ is the reference element.

            arguments   
                obj approx.mesh.Mesh
                EI {mustBeInteger} = []
            end
            
            if ~bitand(obj.Status, obj.DetEJacMask)
                obj.Transform.detEJac = obj.computeAllElementJacobianDeterminants();
                obj.Status = bitor(obj.Status, obj.DetEJacMask);
            end
            
            if isempty(EI) || isscalar(obj.Transform.detEJac)
                detEJac = obj.Transform.detEJac;
            else
                detEJac = obj.Transform.detEJac(EI);
            end
        end
        
        function invEJac = computeElementInverseJacobians(obj, EI)
            % COMPUTEELEMENTINVERSEJacACOBIANS Compute element inverse
            % Jacobian matrices.
            %
            %   invEJac = computeElementInverseJacobians(obj) computes
            %   inverse Jacobian matrices for coordinate transformation
            %   from physical to reference coordinates for all elements.
            %
            %   invEJac = computeElementInverseJacobians(obj, EI) computes
            %   inverse Jacobian matrices only for the specified elements.
            %
            %   The inverse Jacobian matrix appears in gradient
            %   computations: \f$\nabla_x f = J^{-T} \nabla_\xi \f$ where
            %   \f$\nabla_x\f$ and \f$\nabla_\xi\f$ are gradients in
            %   physical and reference coordinates respectively.

            arguments   
                obj approx.mesh.Mesh
                EI {mustBeInteger} = []
            end
            
            if ~bitand(obj.Status, obj.InvEJacMask)
                obj.Transform.invEJac = obj.computeAllElementInverseJacobians();
                obj.Status = bitor(obj.Status, obj.InvEJacMask);
            end
            
            if isempty(EI) || ismatrix(obj.Transform.invEJac)
                invEJac = obj.Transform.invEJac;
            else
                invEJac = obj.Transform.invEJac(:, :, EI);
            end
        end
        
        function detFJac = computeFaceJacobianDeterminants(obj, EI, LFI)
            % COMPUTEFACEJacACOBIANDETERMINANTS Compute face Jacobian
            % determinants.
            %
            %   detFJac = computeFaceJacobianDeterminants(obj, EI, LFI)
            %   computes Jacobian determinants for coordinate
            %   transformation on element faces for integration over
            %   boundaries.
            %
            %   Face Jacobian determinants are used in boundary
            %   integration:
            %
            %   \f[
            %       \int_{\partial\Omega} f(x) dS =
            %       \int_{\partial\hat{\Omega}} f(F(\xi)) |det(J_F)| d\xi
            %   \f]
            %
            %   where \f$\partial\Omega\f$ is the boundary and \f$J_F\f$ is
            %   the face Jacobian.

            arguments   
                obj approx.mesh.Mesh
                EI {mustBeInteger} = []
                LFI {mustBeInteger, mustBePositive} = 1
            end
            
            if ~bitand(obj.Status, obj.DetFJacMask)
                obj.Transform.detFJac = obj.computeAllFaceJacobianDeterminants();
                obj.Status = bitor(obj.Status, obj.DetFJacMask);
            end
            
            if isempty(EI) || isscalar(obj.Transform.detFJac{LFI})
                detFJac = obj.Transform.detFJac{LFI};
            else
                detFJac = obj.Transform.detFJac{LFI}(EI);
            end
        end
        
        function normal = computeOutwardNormals(obj, EI, LFI)
            % COMPUTEOUTWARDNORMALS Compute outward normal vectors.
            %
            %   normal = computeOutwardNormals(obj, EI, LFI) computes
            %   unit outward normal vectors on element faces for boundary
            %   conditions and flux computations.
            %
            %   Outward normals are used in boundary conditions and flux
            %   computations. The direction is defined such that the normal
            %   points away from the element interior.

            arguments   
                obj approx.mesh.Mesh
                EI {mustBeInteger} = []
                LFI {mustBeInteger, mustBePositive} = 1
            end
            
            if ~bitand(obj.Status, obj.OutNormalMask)
                obj.Transform.normal = obj.computeAllOutwardNormals();
                obj.Status = bitor(obj.Status, obj.OutNormalMask);
            end
            
            if isempty(EI) || size(obj.Transform.normal{LFI}, 2) == 1
                normal = obj.Transform.normal{LFI};
            else
                normal = obj.Transform.normal{LFI}(:, EI);
            end
        end
        
        function v = computeDirection(obj)
            % COMPUTEDIRECTION Compute a direction vector.
            %
            %   v = computeDirection(obj) computes a unit direction vector
            %   to define the "plus" and "minus" side of a face.
            %
            %   If nDims = 1, v = 1.
            %   If nDims > 1, v is a unit vector that is not parallel
            %   with any normal vectors of element boundaries.

            arguments   
                obj approx.mesh.Mesh
            end
            
            if ~bitand(obj.Status, obj.DirectionMask)
                obj.Direction = obj.computeUniqueDirection();
                obj.Status = bitor(obj.Status, obj.DirectionMask);
            end
            
            v = obj.Direction;
        end
        
        function ratio = computeFaceScalingFactors(obj, EI, LFI)
            % COMPUTEFACESCALINGFACTORS Compute face scaling factors.
            %
            %   ratio = computeFaceScalingFactors(obj, EI, LFI) computes
            %   face scaling factors for coordinate transformation on
            %   element faces for integration over boundaries.

            arguments   
                obj approx.mesh.Mesh
                EI {mustBeInteger} = []
                LFI {mustBeInteger, mustBePositive} = 1
            end
            
            detEJac = obj.computeElementJacobianDeterminants(EI);
            detFJac = obj.computeFaceJacobianDeterminants(EI, LFI);
            ratio = detFJac ./ detEJac;
        end
    end
    
    methods (Abstract)
        % GRAPHIFY Convert mesh to graph representation.
        graphObj = graphify(obj)
        
        % COLLOCATE Map reference coordinates to physical coordinates.
        Y = collocate(obj, X, EI)
        
        % REFINE Create refined mesh with smaller elements.
        newObj = refine(obj, nLevels)
        
        % COMPUTEMEASURE Compute characteristic mesh size.
        h = computeMeasure(obj)
        
        % COMPUTEALLELEMENTJACOBIANS Compute all element Jacobian matrices.
        EJac = computeAllElementJacobians(obj)
        
        % COMPUTEALLELEMENTJACOBIANDETERMINANTS Compute all element
        % Jacobian determinants.
        detEJac = computeAllElementJacobianDeterminants(obj)
        
        % COMPUTEALLELEMENTINVERSEJacACOBIANS Compute all element inverse
        % Jacobian matrices.
        invEJac = computeAllElementInverseJacobians(obj)
        
        % COMPUTEALLFACEJacACOBIANDETERMINANTS Compute all face Jacobian
        % determinants.
        detFJac = computeAllFaceJacobianDeterminants(obj)
        
        % COMPUTEALLOUTWARDNORMALS Compute all outward normal vectors.
        normal = computeAllOutwardNormals(obj)
        
        % COMPUTEUNIQUEDIRECTION Compute a unique direction vector.
        v = computeUniqueDirection(obj)
    end
end


function mustBeNumericOrCell(x)
if ~iscell(x) && ~isnumeric(x)
    error('Input must be numeric or cell array.');
end
end