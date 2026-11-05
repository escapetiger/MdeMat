classdef UpwindScheme < physics.advection.AdvectionScheme
    % UPWINDSCHEME Upwind scheme.
    %
    %   UpwindScheme implements an upwind scheme for solving advection
    %   equations. The method uses upwind numerical fluxes to handle the
    %   hyperbolic nature of advection problems, providing stability and
    %   sharp resolution of solution features.
    %
    % Notes:
    %   The upwind flux selection provides natural stabilization for
    %   hyperbolic problems but may require smaller time steps compared
    %   to semi-Lagrangian methods for stability.
    %
    % See Also:
    %   physics.OdeScheme, physics.advection.SemiLagrangianScheme

    properties
        assembly %< Upwind assembly
    end

    methods
        function obj = UpwindScheme(config)
            % UWDGSCHEME Constructor for UpwindScheme.
            %
            %   obj = UpwindScheme(config) creates an upwind scheme with
            %   the specified configuration structure. Inherits ODE
            %   integration capabilities from OdeScheme.
            %
            % Inputs:
            %   config - Configuration structure.
            %
            % Outputs:
            %   obj - Constructed UpwindScheme object

            obj@physics.advection.AdvectionScheme(config);
            obj.assembly = [];
        end

        function state = initialize(obj, state)
            % INITIALIZE Initialize the upwind scheme.
            %
            %   state = initialize(obj, state) sets up the upwind scheme by
            %   assembling the advection operator, fitting initial
            %   conditions, and preparing visualization components.
            %
            % Inputs:
            %   obj - The UpwindScheme object
            %   state - Initial solution state structure
            %
            % Outputs:
            %   state - Initialized state ready for time integration

            %< Fit initial condition
            state.addDof('U', state.xDisc.fit([], obj.config.ic));

            %< Create upwind operator
            obj.assembly = approx.assembly.FiniteElementUpwindAssembly(state.xDisc, obj.config.bcType);

            %< Set linear and source terms
            obj.setLinearTerm(state);
            obj.setSourceTerm(state);
        end

        function state = preStep(obj, state)
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

            %< Compute time step based on CFL condition
            h = state.xDisc.mesh.computeMeasure();
            obj.tDisc.setTimeStep(h, obj.config.cfl, []);
        end

        function state = step(obj, state)
            % STEP Perform one time step.
            %
            %   state = step(obj, state) advances the solution by one time
            %   step using the upwind scheme. The method updates the
            %   solution history and applies the advection operator with
            %   appropriate boundary conditions.
            %
            % Inputs:
            %   obj - The UpwindScheme object
            %   state - Current solution state
            %
            % Outputs:
            %   state - Updated solution state after one time step

            %< Update ODE integrator with current solution
            obj.tDisc.update(state.dofs.U);

            %< Perform ODE integration step
            state.dofs.U = obj.tDisc.step(obj.L, obj.S, []);
        end

        function state = postStep(obj, state)
            % POSTSTEP Postprocess operations after time step.
            %
            %   state = postStep(obj, state) performs postprocessing
            %   including visualization updates and timeline advancement
            %   for the next integration step.
            %
            % Inputs:
            %   obj - The UpwindScheme object
            %   state - Current solution state
            %
            % Outputs:
            %   state - Postprocessed solution state

            %< Advance timeline for next step
            obj.tDisc.timeline.advance();
        end

        function state = finalize(obj, state)
            % FINALIZE Finalize simulation with scheme and state reporting.
            %
            %   state = finalize(obj, state) performs final reporting
            %   including scheme configuration details and state
            %   information summary.
            %
            % Inputs:
            %   obj - The UpwindScheme object
            %   state - Final solution state
            %
            % Outputs:
            %   state - Final processed solution state

            %< Report state information
            state.report();

            fprintf('[R] ODE Integrator: %s\n', class(obj.tDisc));
            fprintf('[R] Final Time: %.6g\n', obj.config.tFinal);
            fprintf('[R] CFL Number: %.3f\n', obj.config.cfl);

            %< Report boundary condition type
            if isempty(obj.config.bc)
                fprintf('[R] Boundary Conditions: Periodic\n');
            else
                fprintf('[R] Boundary Conditions: Dirichlet\n');
            end
        end
    
    end

    methods (Access = private) %< Rule 246: helper methods
        function obj = setLinearTerm(obj, state) %#ok<INUSD>
            % SETLINEARTERM Assemble linear operator sparse matrix.

            field = obj.config.advection;
            obj.L = -obj.assembly.assembleMatrix(field);
        end

        function obj = setSourceTerm(obj, ~)
            % SETSOURCETERM Set up source term for boundary conditions.

            switch obj.config.bcType
                case 'periodic'
                    obj.S = [];
                case 'dirichlet'
                    obj.S = @(U, t) obj.computeSourceTerm(U, t);
            end
        end

        function S = computeSourceTerm(obj, ~, t)
            % COMPUTESOURCETERM Compute source term for boundary conditions.

            func = obj.config.bc;
            field = obj.config.advection;
            direction = obj.config.advection;
            S = -obj.assembly.assembleVector(field, direction, func, t);
        end
    end
end