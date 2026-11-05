classdef UwdgScheme < physics.Scheme
    % UWDGSCHEME Upwind Discontinuous Galerkin scheme.
    %
    %   UwdgScheme implements an Upwind Discontinuous Galerkin (UWDG)
    %   method for solving advection equations. The method uses upwind
    %   numerical fluxes to handle the hyperbolic nature of advection
    %   problems, providing stability and sharp resolution of solution
    %   features.
    %
    %   The UWDG method discretizes the spatial domain using discontinuous
    %   Galerkin finite elements and handles inter-element communication
    %   through upwind fluxes, making it particularly suitable for
    %   transport-dominated problems with sharp gradients or
    %   discontinuities.
    %
    % Examples:
    %   % Create UWDG scheme with configuration
    %   config = struct('ic', @(x) sin(x), 'advection', [1], 'cfl', 0.5, ...
    %                   'bc', [], 'verbose', true);
    %   scheme = UwdgScheme(config, timeIntegrator, visualizer, analyzer);
    %   
    %   % Full simulation workflow
    %   state = scheme.initialize(initialState);
    %   finalState = scheme.run(state);
    %   scheme.finalize(finalState);
    %
    % Notes:
    %   The upwind flux selection provides natural stabilization for
    %   hyperbolic problems but may require smaller time steps compared
    %   to semi-Lagrangian methods for stability.
    %
    % See Also:
    %   physics.Scheme, physics.advection.SldgScheme,
    %   fem.assembly.UpwindGridAssembly

    properties
        assembly % Transport assembly for upwind operations
        operator % Transport operator matrix
    end

    methods
        function state = initialize(obj, state)
            % INITIALIZE Initialize the UWDG scheme.
            %
            %   state = initialize(obj, state) sets up the upwind
            %   discontinuous Galerkin scheme by configuring finite element
            %   data, assembling the advection operator, projecting initial
            %   conditions, and preparing visualization components.
            %
            % Inputs:
            %   obj - The UwdgScheme object
            %   state - Initial solution state structure
            %
            % Outputs:
            %   state - Initialized state ready for time integration

            %< Reset timer
            obj.timer.reset();

            %< Local data
            state.disc.x.fe.setVolumeData(1);
            state.disc.x.fe.setFluxData(0);
            state.disc.x.op.setVolumeData();
            state.disc.x.op.setFluxData();

            %< Initial condition
            state.dofs.U = state.disc.x.space.project(obj.config.ic);

            %< Advection operator
            obj.assembly = fem.assembly.UpwindGridAssembly( ...
                state.disc.x.fe, state.disc.x.mesh, state.disc.x.op, ...
                ~isempty(obj.config.bc));
            obj.operator = obj.assembly.divergence(obj.config.advection);

            %< Reset time discretization
            obj.tDisc.reset();

            %< Reset visualizer
            obj.visualizer.reset();
            obj.visualizer.addDataset('DG', { ...
                'Color', 'b', ...
                'Marker', 'o', ...
                'LineStyle', 'none'});
            if obj.visualizer.hasExact
                obj.visualizer.addDataset('REF', ...
                    {'Color', 'r', ...
                    'Marker', 'none', ...
                    'LineStyle', '-'});
            end

            %< Visualize
            if obj.visualizer.isEnabled
                obj.timer.start(2, '[M] Visualization.');
                obj.visualizer.plot(state.disc.x.space, state.dofs, obj.tDisc);
                obj.timer.stop(2);
            end
        end

        function state = preStep(obj, state)
            % PRESTEP Preprocess operations before time step.
            %
            %   state = preStep(obj, state) performs preprocessing
            %   including time step size computation based on CFL condition
            %   and coefficient updates for multi-step methods.
            %
            % Inputs:
            %   obj - The UwdgScheme object
            %   state - Current solution state
            %
            % Outputs:
            %   state - Preprocessed solution state

            h = state.disc.x.mesh.measure;
            obj.tDisc.timeline.setTimeStep(h, obj.config.cfl);
            if isa(obj.tDisc, 'approx.odeint.BdfIntegrator')
                obj.tDisc.setCoefficients();
            end
        end

        function state = step(obj, state)
            % STEP Perform one upwind DG time step.
            %
            %   state = step(obj, state) advances the solution by one time
            %   step using the upwind discontinuous Galerkin method. The
            %   method updates the solution history and applies the linear
            %   advection operator with appropriate boundary conditions.
            %
            % Inputs:
            %   obj - The UwdgScheme object
            %   state - Current solution state
            %
            % Outputs:
            %   state - Updated solution state after one time step

            if ~isempty(obj.config.bc)
                c = obj.config.advection;
                S = @(t) -obj.assembly.divergenceBc(c, c, obj.config.bc, t);
            else
                S = [];
            end

            obj.tDisc.update(state.dofs.U);
            L = -obj.operator;
            state.dofs.U = obj.tDisc.step(L, S, []);
        end

        function state = postStep(obj, state)
            % POSTSTEP Postprocess operations after time step.
            %
            %   state = postStep(obj, state) performs postprocessing
            %   including visualization updates and timeline advancement
            %   for the next integration step.
            %
            % Inputs:
            %   obj - The UwdgScheme object
            %   state - Current solution state
            %
            % Outputs:
            %   state - Postprocessed solution state

            if obj.visualizer.isEnabled
                obj.timer.start(2, '[M] Visualization.');
                obj.visualizer.plot(state.disc.x.space, state.dofs, obj.tDisc);
                obj.timer.stop(2);
            end

            obj.tDisc.timeline.advance();
        end

        function state = finalize(obj, varargin)
            % FINALIZE Finalize simulation with error analysis.
            %
            %   state = finalize(obj, state) performs final error analysis
            %   for absolute error computation against exact solutions.
            %
            %   [fine, coarse] = finalize(obj, fine, coarse) performs
            %   relative error analysis comparing fine and coarse grid
            %   solutions for convergence studies.
            %
            % Inputs:
            %   obj - The UwdgScheme object
            %   varargin - Variable arguments containing state(s) for analysis
            %
            % Outputs:
            %   state - Final processed solution state

            if obj.analyzer.isEnabled
                obj.timer.start(3, '[M] Error Evaluation.');
                if length(varargin) == 1
                    state = varargin{1};
                    obj.analyzer.setLevel(state.disc.x.mesh);
                    obj.analyzer.absolute(state.disc.x.space, state.dofs, ...
                        obj.tDisc.timeline.now);
                else
                    [state, coarse] = varargin{:};
                    obj.analyzer.setLevel(state.disc.x.mesh);
                    obj.analyzer.relative(state.disc.x.space, state.dofs, ...
                        coarse.disc.x.space, coarse.dofs);
                end
                obj.timer.stop(3);
            end

            obj.timer.report();
        end
    end
end