classdef SemiLagrangianScheme < physics.advection.AdvectionScheme
    % SEMILAGRANGIANSCHEME Semi-Lagrangian scheme.
    %
    %   SemiLagrangianScheme implements a Semi-Lagrangian (SL) method for
    %   solving advection equations. The method combines the stability
    %   advantages of Lagrangian methods with the flexibility of spatial
    %   discretizations.
    %
    %   The SL method traces characteristics backward in time and uses
    %   high-order interpolation to evaluate solutions at departure points,
    %   making it particularly effective for transport-dominated problems
    %   with large time steps and minimal numerical diffusion.
    %
    % Notes:
    %   The semi-Lagrangian approach allows for larger time steps compared
    %   to traditional Eulerian methods, making it computationally
    %   efficient for advection-dominated flows.
    %
    % See Also:
    %   physics.OdeScheme, physics.advection.UpwindScheme

    properties
        assembly % Semi-Lagrangian assembly
    end

    methods
        function obj = SemiLagrangianScheme(config)
            % SEMILAGRANGIANSCHEME Constructor for SemiLagrangianScheme.
            %
            %   obj = SemiLagrangianScheme(config) creates a SL scheme with
            %   the specified configuration structure. Inherits ODE
            %   integration capabilities from OdeScheme.
            %
            % Inputs:
            %   config - Configuration structure.
            %
            % Outputs:
            %   obj - Constructed SemiLagrangianScheme object

            obj@physics.advection.AdvectionScheme(config);
            obj.assembly = [];
        end

        function state = initialize(obj, state)
            % INITIALIZE Initialize the SL scheme.
            %
            %   state = initialize(obj, state) sets up the semi-Lagrangian
            %   discontinuous Galerkin scheme by configuring finite element
            %   data, assembling the advection operator, fitting initial
            %   conditions, and preparing visualization components.
            %
            % Inputs:
            %   obj - The SemiLagrangianScheme object
            %   state - Initial solution state structure
            %
            % Outputs:
            %   state - Initialized state ready for time integration

            %< Fit initial condition
            state.addDof('U', state.xDisc.fit([], obj.config.ic));

            %< Create Semi-Lagrangian assembly
            obj.assembly = approx.assembly.FiniteElementSemiLagrangianAssembly(state.xDisc, obj.config.bcType);

            %< Set mass terms
            obj.setMassTerm(state);
        end

        function state = preStep(obj, state)
            % PRESTEP Preprocess operations before time step.
            %
            %   state = preStep(obj, state) performs preprocessing
            %   including time step size computation based on CFL condition
            %   and coefficient updates for multi-step methods.
            %
            % Inputs:
            %   obj - The SemiLagrangianScheme object
            %   state - Current solution state
            %
            % Outputs:
            %   state - Preprocessed solution state

            %< Compute time step based on CFL condition
            h = state.xDisc.mesh.computeMeasure();
            obj.tDisc.setTimeStep(h, obj.config.cfl, []);
        end


        function state = step(obj, state)
            % STEP Perform one semi-Lagrangian time step.
            %
            %   state = step(obj, state) advances the solution by one time
            %   step using the semi-Lagrangian transport operator. The
            %   method updates the solution history and applies the
            %   characteristic-based transport with appropriate boundary
            %   condition.
            %
            % Inputs:
            %   obj - The SemiLagrangianScheme object
            %   state - Current solution state
            %
            % Outputs:
            %   state - Updated solution state after one time step

            %< Update ODE integrator with current solution
            obj.tDisc.update(state.dofs.U);

           %< Perform ODE integration step
            state.dofs.U = obj.tDisc.step([], [], obj.M);
        end

        function state = postStep(obj, state)
            % POSTSTEP Postprocess operations after time step.
            %
            %   state = postStep(obj, state) performs postprocessing
            %   including visualization updates and timeline advancement
            %   for the next integration step.
            %
            % Inputs:
            %   obj - The SemiLagrangianScheme object
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
            %   obj - The SemiLagrangianScheme object
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
        function obj = setMassTerm(obj, ~)
            % SETMASSTERM Set up mass term.

            obj.M = @(tBegin, tEnd) obj.computeM(tBegin, tEnd);
        end

        function [V, M] = computeM(obj, tBegin, tEnd)
            % COMPUTEMASSTERM Compute mass term.

            field = obj.config.advection;
            
            obj.assembly.reset();

            if nargout >= 2
                M = obj.assembly.assembleMatrix(field, tBegin, tEnd);
            end

            if nargout >= 1
                V = [];
                if strcmpi(obj.config.bcType, 'dirichlet')
                    if abs(tBegin-tEnd) < 1e-8, return; end
                    V = obj.assembly.assembleVector(field, tBegin, tEnd, obj.config.bc);
                end
            end
        end
    end
end