classdef LocalScheme < physics.diffusion.DiffusionScheme
    % LOCALSCHEME Local scheme. 

    properties
        assembly % Diffusion assembly
    end

    methods
        function obj = LocalScheme(config)
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

            obj@physics.diffusion.DiffusionScheme(config);
            obj.assembly = [];
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

            % Initial condition
            state.addDof('U', state.xDisc.fit([], obj.config.ic));

            % Diffusion operator
            obj.assembly = approx.assembly.FiniteElementEllipticAssembly( ...
                state.xDisc, obj.config.bcType, obj.config.xPenaltyType);

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
            %   solution history and applies the diffusion operator with
            %   appropriate boundary conditions.
            %
            % Inputs:
            %   obj - The LocalScheme object
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
            %   obj - The LocalScheme object
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
            %   obj - The LocalScheme object
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

            field = obj.config.diffusion;
            obj.L = obj.assembly.assembleMatrix(field);
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
            field = obj.config.diffusion;
            S = obj.assembly.assembleVector(field, func, t);
        end
    end
end
