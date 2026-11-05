classdef FiniteElementDiscretization < handle
    % FINITEELEMENTDISCRETIZATION Finite element discretization.
    %
    %   FiniteElementDiscretization contains all components needed for
    %   finite element methods: element, mesh, space, and operator.
    %   It provides a complete framework for discontinuous Galerkin (DG)
    %   and continuous finite element methods on orthotopal meshes.
    %
    %   The class manages the complete finite element discretization
    %   including geometric mesh, function spaces, and differential
    %   operators. It supports various basis function types (Legendre,
    %   Lagrange) and integration schemes for high-order accurate
    %   computations.
    %
    % Examples:
    %   % Create FE discretization with DG orthotope element
    %   disc = FiniteElementDiscretization();
    %   disc.setDgOrthotopeElement('nDims', 2, 'order', 3, 'basisType', 1);
    %   disc.setElementOperator();
    %   disc.setDerivativeOrder(1);
    %   
    %   % Set up uniform mesh and function space
    %   disc.setUniformGrid([10, 10], [0, 1; 0, 1]);
    %   disc.setMeshSpace();
    %
    %   % Access gradient operator for stiffness matrix assembly
    %   gradOp = disc.operator.gradient;
    %   volGradData = gradOp.volumeData;   % Volume gradient terms
    %   fluxGradData = gradOp.fluxData;    % Boundary flux terms
    %   
    %   % Mesh refinement
    %   fineMesh = disc.refineMesh(2);  % Refine by factor of 2
    %
    % Notes:
    %   Components must be set in order: element, then operator, then
    %   derivative order. Mesh and space can be set independently but
    %   both are needed for mesh-based computations.
    %
    % See also:
    %   approx.element.DgElement, approx.element.DgElementOperator,
    %   approx.mesh.UniformGrid, approx.mesh.NonuniformGrid,
    %   approx.space.MeshSpace


    properties (Constant)
        AVAILABLE_BASIS_TYPES = {'modal_legendre', 'nodal_gauss_legendre', 'nodal_gauss_lobatto'}
        AVAILABLE_BASIS_PATTERNS = {'Q', 'P'}
    end

    properties
        element  % Element object for local computations
        mesh     % Mesh object defining computational domain
        space    % Function space object spanning mesh
        operator % Element operator object for differential operations
    end

    methods
        function obj = FiniteElementDiscretization()
            % FINITEELEMENTDISCRETIZATION Constructor for
            % FiniteElementDiscretization.
            %
            %   obj = FiniteElementDiscretization() creates an empty finite
            %   element discretization object. Components must be set using
            %   the appropriate setter methods in the correct order.
            %
            % Outputs:
            %   obj - Constructed FiniteElementDiscretization object

            obj.element = [];
            obj.mesh = [];
            obj.space = [];
            obj.operator = [];
        end

        function obj = setDgOrthotopeElement(obj, varargin)
            % SETDGORTHOTOPEELEMENT Set DG orthotope element with specified
            % configuration.
            %
            %   obj = setDgOrthotopeElement(obj, varargin) creates a
            %   discontinuous finite element on the unit orthotope with
            %   specified polynomial order and basis configuration. This
            %   method encapsulates the complex element construction
            %   process including geometry, basis functions, integrators,
            %   and projectors for DG methods.
            %
            % Inputs:
            %   obj - The FiniteElementDiscretization object
            %   varargin - Name-value parameter pairs:
            %<    'nDims' - Number of spatial dimensions (positive integer)
            %<    'order' - Polynomial order (positive integer, default: 1)
            %<    'basisType' - Basis function type (1, 2, or 3, default: 1)
            %<    'basisPattern' - Basis ordering pattern (1 or 2, default: 1)
            %
            % Outputs:
            %   obj - The FiniteElementDiscretization object
            %
            % Notes:
            %   Basis types:
            %     'model_legendre': Monic Legendre basis
            %     'nodal_gauss_legnedre': Lagrange-Gauss-Legendre basis
            %     'nodal_gauss_lobatto': Lagrange-Gauss-Lobatto basis
            %
            %   Basis patterns:
            %     'lx': L_infinity norm ordering (tensor product)
            %     'l1': L_1 norm ordering (total degree)

            p = inputParser;
            addParameter(p, 'nDims', [], @(x) x >= 1);
            addParameter(p, 'order', 1, @(x) x >= 1);
            addParameter(p, 'basisType', 'legendre', ...
                @(x) ismember(x, obj.AVAILABLE_BASIS_TYPES));
            addParameter(p, 'basisPattern', 1, ...
                @(x) ismember(x, obj.AVAILABLE_BASIS_PATTERNS));
            parse(p, varargin{:});

            d = p.Results.nDims;
            k = p.Results.order;
            t = p.Results.basisType;
            s = p.Results.basisPattern;

            %< Create unit orthotope geometry
            G = core.geometry.Orthotope.unit(d);

            %< Create basis functions based on type
            switch lower(t)
                case 'modal_legendre'
                    U = approx.basis.OrthogonalBasisFunction( ...
                        k, 'monic_legendre', 'unit');
                case 'nodal_gauss_legendre'
                    U = approx.basis.InterpolationBasisFunction( ...
                        k, 'lagrange', 'unit', 'gauss_legendre');
                case 'nodal_gauss_lobatto'
                    U = approx.basis.InterpolationBasisFunction( ...
                        k, 'lagrange', 'unit', 'gauss_lobatto');
            end
            U = repmat(U, 1, d);

            %< Create separable basis with specified pattern
            switch upper(s)
                case 'Q'
                    B = approx.basis.SeparableBasisFunction(U, 'lx');
                case 'P'
                    B = approx.basis.SeparableBasisFunction(U, 'l1');
            end
            B.autoLoad();

            %< Create volume integrator
            switch lower(t)
                case {'modal_legendre', 'nodal_gauss_legendre'}
                    Q = 'gauss_legendre';
                case 'nodal_gauss_lobatto'
                    Q = 'gauss_lobatto';
            end
            VI = approx.integrate.OrthotopeIntegrator(G.nDims, Q);
            VI.setPoints(B.nFactorCodims, G.lower, G.upper);

            %< Create face integrators for all boundaries
            FI = arrayfun(@(i) approx.integrate.OrthotopeFaceIntegrator( ...
                G.nDims, Q, i), 1:(2 * G.nDims));
            for i = 1:(2 * G.nDims)
                FI(i).setPoints(B.nFactorCodims, G.lower, G.upper);
            end

            %< Create projector based on basis type
            switch lower(t)
                case 'modal_legendre'
                    P = approx.project.ModalProjector(B);
                    P.setMass(VI.nodes, VI.weights);
                case {'nodal_gauss_legendre', 'nodal_gauss_lobatto'}
                    P = approx.project.NodalProjector(B);
                    P.setMass(VI.nodes);
            end

            %< Construct DG orthotope finite element
            obj.element = approx.element.DgElement(G, P, VI, FI);
        end

        function obj = setElementOperator(obj)
            % SETELEMENTOPERATOR Set element operator for differential
            % computations.
            %
            %   obj = setElementOperator(obj) creates an element operator
            %   appropriate for the current element type. The operator
            %   provides gradient and other differential operators in both
            %   volume and flux (boundary) forms required for finite
            %   element methods.
            %
            % Inputs:
            %   obj - The FiniteElementDiscretization object
            %
            % Outputs:
            %   obj - The FiniteElementDiscretization object

            core.except.assert(~isempty(obj.element), ...
                'NotPrepared', 'Element must be set before operator.');

            switch class(obj.element)
                case 'approx.element.C0Element'
                    obj.operator = approx.element.C0ElementOperator(obj.element);
                case 'approx.element.DgElement'
                    obj.operator = approx.element.DgElementOperator(obj.element);
            end
        end

        function obj = setDerivativeOrder(obj, order)
            % SETDERIVATIVEORDER Set derivative order for all components.
            %
            %   obj = setDerivativeOrder(obj, order) sets up function data
            %   for volume and flux integrals with the specified maximum
            %   derivative order. This method initializes all necessary
            %   data structures for computing bilinear forms involving
            %   derivatives up to the specified order.
            %
            % Inputs:
            %   obj - The FiniteElementDiscretization object
            %   order - Maximum derivative order (non-negative integer)
            %
            % Outputs:
            %   obj - The FiniteElementDiscretization object

            core.except.assert(~isempty(obj.element) && ~isempty(obj.operator), ...
                'NotPrepared', 'Element and operator must be set before setting derivative order.');

            %< Set volume data
            obj.element.setVolumeData(order);

            %< Set flux data for DG element
            if isa(obj.element, 'approx.element.DgElement')
                obj.element.setFluxData(max(0, order-1));
            end

            %< Set volume data
            obj.operator.setVolumeData();

            %< Set flux data for DG operator
            if isa(obj.operator, 'approx.element.DgElementOperator')
                obj.operator.setFluxData();
            end
        end

        function obj = setUniformGrid(obj, n, bbox)
            % SETUNIFORMGRID Set uniform grid mesh.
            %
            %   obj = setUniformGrid(obj, n, bbox) creates a uniform grid
            %   mesh with the specified resolution and bounding box. The
            %   mesh provides the computational domain for finite element
            %   computations.
            %
            % Inputs:
            %   obj - The FiniteElementDiscretization object
            %   n - Grid resolution vector (positive integers)
            %   bbox - Bounding box matrix [lower; upper] (2 x nDims)
            %
            % Outputs:
            %   obj - The FiniteElementDiscretization object

            obj.mesh = approx.mesh.UniformGrid(n, bbox);
        end

        function obj = setNonuniformGrid(obj, x)
            % SETNONUNIFORMGRID Set nonuniform grid mesh.
            %
            %   obj = setNonuniformGrid(obj, x) creates a nonuniform grid
            %   mesh with the specified coordinate vectors. This allows for
            %   adaptive mesh refinement and custom domain discretizations.
            %
            % Inputs:
            %   obj - The FiniteElementDiscretization object
            %   x - Cell array of coordinate vectors, one per dimension
            %
            % Outputs:
            %   obj - The FiniteElementDiscretization object

            obj.mesh = approx.mesh.NonuniformGrid(x);
        end

        function obj = setMeshSpace(obj)
            % SETMESHSPACE Set mesh-based function space.
            %
            %   obj = setMeshSpace(obj) creates a mesh-based function space
            %   that spans the entire mesh domain. The space provides the
            %   global degrees of freedom for finite element computations.
            %
            % Inputs:
            %   obj - The FiniteElementDiscretization object
            %
            % Outputs:
            %   obj - The FiniteElementDiscretization object

            core.except.assert(~isempty(obj.element) && ~isempty(obj.mesh), ...
                'NotPrepared', 'Element and mesh must be set before space.');

            obj.space = approx.space.MeshSpace(obj.element, obj.mesh);
        end

        function newObj = refineMesh(obj, varargin)
            % REFINEMESH Create refined mesh discretization.
            %
            %   newObj = refineMesh(obj, varargin) creates a new finite
            %   element discretization with refined mesh while preserving
            %   the element and operator configurations. This enables
            %   h-refinement studies and adaptive mesh refinement.
            %
            % Inputs:
            %   obj - The FiniteElementDiscretization object
            %   varargin - Refinement parameters:
            %<    For uniform grid: k - Refinement factor per dimension
            %<    For nonuniform grid: k - Number of subdivisions per element
            %
            % Outputs:
            %   newObj - New FiniteElementDiscretization with refined mesh

            newObj = approx.discretization.FiniteElementDiscretization();
            newObj.element = obj.element;
            newObj.mesh = obj.mesh.refine(varargin{:});
            newObj.space = approx.space.MeshSpace(obj.element, newObj.mesh);
            newObj.operator = obj.operator;
        end
    end
end