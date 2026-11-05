classdef PoissonScheme < physics.SteadyScheme
    % POISSONSCHEME Base class for all schemes solving poisson equations.
    %
    %   PoissonScheme implements a numerical scheme for solving poisson
    %   equations. This base class defines the basic validation and parse
    %   logic of configuration struct.
    %
    % See Also:
    %   physics.SteadyScheme

    properties
        solver % Linear solver
        elliptic % Elliptic assembly
        source % Source assembly
        A % LHS matrix
        b % RHS vector
    end

    methods
        function obj = PoissonScheme(config)
            % LOCALSCHEME Constructor for LocalScheme.
            %
            %   obj = LocalScheme(config) creates a local scheme with the
            %   specified configuration structure. Inherits ODE integration
            %   capabilities from OdeScheme.
            %
            % Inputs:
            %   config - Configuration structure.
            %
            % Outputs:
            %   obj - Constructed LocalScheme object

            obj@physics.SteadyScheme(config);
            obj.solver = core.linalg.LinearSolver();
            obj.elliptic = [];
        end

        function state = initialize(obj, state)
            % INITIALIZE Initialize the local scheme.
            %
            %   state = initialize(obj, state) sets up the local scheme by
            %   assembling the local diffusion operator, fitting initial
            %   conditions, and preparing visualization components.
            %
            % Inputs:
            %   obj - The LocalScheme object
            %   state - Initial solution state structure
            %
            % Outputs:
            %   state - Initialized state ready for time integration

            % Add dofs
            state.addDof('U', []);
            state.addCoefficient('force', state.xDisc.mean([], obj.config.force));

            % Elliptic operator
            obj.elliptic = approx.assembly.FiniteElementEllipticAssembly( ...
                state.xDisc, obj.config.bcType, obj.config.xPenaltyType);

            % Source operator
            obj.source = approx.assembly.FiniteElementSourceAssembly(state.xDisc);

            %< Set LHS and RHS terms
            obj.setLhsTerm(state);
            obj.setRhsTerm(state);
        end

        function state = preStep(~, state)
            % PRESTEP Preprocess operations before time step.
            %
            %   state = preStep(obj, state) performs preprocessing
            %   including time step size computation based on CFL condition
            %   and coefficient updates for multi-step methods.
            %
            % Inputs:
            %   obj - The UpwindScheme object
            %   state - Current solution state
            %
            % Outputs:
            %   state - Preprocessed solution state
        end

        function state = step(obj, state)
            % STEP Perform one time step.
            %
            %   state = step(obj, state) advances the solution by one time
            %   step using the upwind scheme. The method updates the
            %   solution history and applies the diffusion operator with
            %   appropriate boundary conditions.
            %
            % Inputs:
            %   obj - The LocalScheme object
            %   state - Current solution state
            %
            % Outputs:
            %   state - Updated solution state after one time step

            switch obj.config.bcType
                case 'periodic'
                    [obj.A, obj.b] = obj.elliptic.forceCompletion(obj.A, obj.b);
                    U = obj.solver.solve(obj.A, obj.b);
                    state.dofs.U = U(1:end-1);
                case 'dirichlet'
                    state.dofs.U = obj.solver.solve(obj.A, obj.b);
            end
        end

        function state = postStep(~, state)
            % POSTSTEP Postprocess operations after time step.
            %
            %   state = postStep(obj, state) performs postprocessing
            %   including visualization updates and timeline advancement
            %   for the next integration step.
            %
            % Inputs:
            %   obj - The LocalScheme object
            %   state - Current solution state
            %
            % Outputs:
            %   state - Postprocessed solution state

        end
    
        function state = finalize(obj, state)
            % FINALIZE Finalize simulation with scheme and state reporting.
            %
            %   state = finalize(obj, state) performs final reporting
            %   including scheme configuration details and state
            %   information summary.
            %
            % Inputs:
            %   obj - The LocalScheme object
            %   state - Final solution state
            %
            % Outputs:
            %   state - Final processed solution state

            %< Report state information
            state.report();

            %< Report boundary condition type
            if isempty(obj.config.bc)
                fprintf('[R] Boundary Conditions: Periodic\n');
            else
                fprintf('[R] Boundary Conditions: Dirichlet\n');
            end
        end
    end

    methods (Access = protected)
        function obj = parse(obj)
            % PARSE Check and parse the contained configuration struct.
            %
            %   obj = parse(obj) checks and parses the contained
            %   configuration struct.
            %
            % Inputs:
            %   obj - The DiffusionScheme object
            %
            % Outputs:
            %   obj - The DiffusionScheme object

            parse@physics.SteadyScheme(obj);

            %< Validate required fields
            requiredFields = {'bc'};
            [tf, mf] = core.validate.hasFields(obj.config, requiredFields);
            core.except.assert(tf, 'InvalidInput', ...
                'Required configuration %s is missing.', mf);

            %< Determine boundary condition type
            if isempty(obj.config.bc)
                obj.config.bcType = 'periodic';
            else
                obj.config.bcType = 'dirichlet';
            end
        end    
    
        function obj = setLhsTerm(obj, state)
            % SETLINEARTERM Assemble linear operator sparse matrix.

            field = eye(state.xDisc.nDims);
            obj.A = -obj.elliptic.assembleMatrix(field);
        end

        function obj = setRhsTerm(obj, state)
            % SETSOURCETERM Set up source term.

            obj.b = 0;
            obj.addForce(state);
            obj.addBc(state);
        end

        function obj = addForce(obj, state)
            G = obj.source.assemblePointwiseVector(obj.config.force);
            obj.b = obj.b + G;
        end

        function obj = addBc(obj, state)
            if strcmpi(obj.config.bcType, 'periodic'), return; end
            func = obj.config.bc;
            field = eye(state.xDisc.nDims);
            G = obj.elliptic.assembleVector(field, func);
            obj.b = obj.b + G;
        end
    end
end