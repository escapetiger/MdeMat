classdef SemiLagrangianFiniteElementDiscretization < approx.discretization.FiniteElementDiscretization
    % SEMILAGRANGIANFINITEELEMENTDISCRETIZATION Semi-Lagrangian finite
    % element discretization.
    %
    %   SemiLagrangianFiniteElementDiscretization extends the standard
    %   finite element discretization with semi-Lagrangian capabilities for
    %   solving advection-dominated problems. It adds element clipping and
    %   semi-Lagrangian operators to handle characteristic tracing and flux
    %   computations across element boundaries.
    %
    %   The semi-Lagrangian method traces characteristics backward in time
    %   and requires specialized operators to handle the interaction
    %   between the departure points and the finite element mesh structure.
    %
    % Examples:
    %   % Create semi-Lagrangian DG discretization
    %   disc = SemiLagrangianFiniteElementDiscretization();
    %   disc.setDgOrthotopeElement('nDims', 2, 'order', 2);
    %   disc.setElementOperator();
    %   disc.setElementClipper();
    %   disc.setSemiLagrangianElementOperator();
    %   
    %   % Set up mesh and advance solution
    %   disc.setUniformGrid([20, 20], [0, 1; 0, 1]);
    %   disc.setMeshSpace();
    %   disc.setDerivativeOrder(1);
    %   
    %   % Access semi-Lagrangian operator for time stepping
    %   slOp = disc.slOperator;
    %   clipper = disc.clipper;
    %   
    %   % Mesh refinement preserves all components
    %   fineMesh = disc.refineMesh(2);
    %
    % Notes:
    %   Element and clipper must be set before semi-Lagrangian operator.
    %   The clipper handles geometric intersections for characteristic
    %   foot computations.
    %
    % See also:
    %   approx.discretization.FiniteElementDiscretization,
    %   approx.element.DgElementClipper, approx.element.C0ElementClipper,
    %   approx.element.DgSemiLagrangianElementOperator
    
    properties
        clipper    % Element clipper for geometric intersections
        slOperator % Semi-Lagrangian element operator for advection
    end
    
    methods
        function obj = setElementClipper(obj, varargin)
            % SETELEMENTCLIPPER Set element clipper for geometric
            % operations.
            %
            %   obj = setElementClipper(obj, varargin) creates an element
            %   clipper appropriate for the current element type. The
            %   clipper handles geometric intersections needed for
            %   semi-Lagrangian characteristic tracing and flux
            %   computations.
            %
            % Inputs:
            %   obj - The SemiLagrangianFiniteElementDiscretization object
            %   varargin - Additional parameters passed to clipper constructor
            %
            % Outputs:
            %   obj - The SemiLagrangianFiniteElementDiscretization object

            core.except.assert(~isempty(obj.element), ...
                'NotPrepared', 'Element must be set before clipper.');
                
            switch class(obj.element)
                case 'approx.element.C0Element'
                    obj.clipper = approx.element.C0ElementClipper(obj.element, varargin{:});
                case 'approx.element.DgElement'
                    obj.clipper = approx.element.DgElementClipper(obj.element, varargin{:});
            end
        end

        function obj = setSemiLagrangianElementOperator(obj)
            % SETSEMILAGRANGIANELEMENTOPERATOR Set semi-Lagrangian element
            % operator.
            %
            %   obj = setSemiLagrangianElementOperator(obj) creates a
            %   semi-Lagrangian element operator that provides specialized
            %   operations for advection problems including characteristic
            %   tracing and flux reconstruction across element boundaries.
            %
            % Inputs:
            %   obj - The SemiLagrangianFiniteElementDiscretization object
            %
            % Outputs:
            %   obj - The SemiLagrangianFiniteElementDiscretization object

            core.except.assert(~isempty(obj.element) && ~isempty(obj.clipper), ...
                'NotPrepared', 'Element and clipper must be set before semi-Lagrangian operator.');

            switch class(obj.element)
                case 'approx.element.C0Element'
                    obj.slOperator = approx.element.SemiLagrangianC0ElementOperator(obj.element, obj.clipper);
                case 'approx.element.DgElement'
                    obj.slOperator = approx.element.SemiLagrangianDgElementOperator(obj.element, obj.clipper);
            end
        end

        function newObj = refineMesh(obj, varargin)
            % REFINEMESH Create refined semi-Lagrangian discretization.
            %
            %   newObj = refineMesh(obj, varargin) creates a new
            %   semi-Lagrangian finite element discretization with refined
            %   mesh while preserving all component configurations
            %   including element, operator, clipper, and semi-Lagrangian
            %   operator.
            %
            % Inputs:
            %   obj - The SemiLagrangianFiniteElementDiscretization object
            %   varargin - Refinement parameters:
            %<    For uniform grid: k - Refinement factor per dimension
            %<    For nonuniform grid: k - Number of subdivisions per element
            %
            % Outputs:
            %   newObj - New SemiLagrangianFiniteElementDiscretization with refined mesh

            newObj = approx.discretization.SemiLagrangianFiniteElementDiscretization();
            newObj.element = obj.element;
            newObj.mesh = obj.mesh.refine(varargin{:});
            newObj.space = approx.space.MeshSpace(obj.element, newObj.mesh);
            newObj.operator = obj.operator;
            newObj.clipper = obj.clipper;
            newObj.slOperator = obj.slOperator;
        end
    end
end