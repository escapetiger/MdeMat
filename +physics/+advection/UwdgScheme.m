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
    %   approx.assembly.UpwindGridAssembly

    properties
        operator % Upwind advection operator
        A % Advection term 
        L % Linear term
        S % Source term
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

            %< Initial condition
            state.dofs.U = state.xDisc.space.project(obj.config.ic);

            %< Advection operator
            if isempty(obj.config.bc)
                bcType = 'periodic';
            else
                bcType = 'dirichlet';
            end
            obj.operator = approx.assembly.ConvectionOperator( ...
                state.xDisc.space, state.xDisc.operator, bcType);

            %< Linear term
            obj.setLinearTerm(state);

            %< Source term
            obj.setSourceTerm(state);

            %< Reset time discretization
            obj.tDisc.reset();

            %< Reset visualizer
            if ~isempty(obj.visualizer) && obj.visualizer.isEnabled
                if ~isempty(obj.config.exact)
                    obj.visualizer.addExact('U', obj.config.exact);
                end
                obj.visualizer.reset();

                obj.timer.start(2, '[M] Visualization.');
                obj.visualizer.plot(state.xDisc.space, state.dofs, obj.tDisc);
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

            h = state.xDisc.mesh.measure;
            obj.tDisc.setTimeStep(h, obj.config.cfl, []);
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

            obj.tDisc.update(state.dofs.U);
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
            %   obj - The UwdgScheme object
            %   state - Current solution state
            %
            % Outputs:
            %   state - Postprocessed solution state

            if ~isempty(obj.visualizer) && obj.visualizer.isEnabled
                obj.timer.start(2, '[M] Visualization.');
                obj.visualizer.plot(state.xDisc.space, state.dofs, obj.tDisc);
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

            if length(varargin) == 1
                state = varargin{1};
                coarse = [];
            else
                [state, coarse] = varargin{:};
            end

            if ~isempty(obj.analyzer) && obj.analyzer.isEnabled
                obj.timer.start(3, '[M] Error Evaluation.');
                if isempty(coarse)
                    obj.analyzer.addExact('U', obj.config.exact);
                    obj.analyzer.setLevel(state.xDisc.mesh);
                    obj.analyzer.absolute(state.xDisc.space, state.dofs, ...
                        obj.tDisc.timeline.now);
                else
                    obj.analyzer.setLevel(state.xDisc.mesh);
                    obj.analyzer.relative(state.xDisc.space, state.dofs, ...
                        coarse.xDisc.space, coarse.dofs);
                end
                obj.timer.stop(3);
            end

            obj.timer.report();
        end
    end

    methods
        function obj = setLinearTerm(obj, state)
            % Assemble linear operator sparse matrix

            obj.A = obj.operator.linear(obj.config.advection);

            obj.L = -obj.A;
        end

        function obj = setSourceTerm(obj, ~)
            if ~isempty(obj.config.bc)
                obj.S = @(U, t) computeS(U, t);
            else
                obj.S = [];
            end

            function S = computeS(~, t)
                c = obj.config.advection;
                S = -obj.operator.linearBc(c, c, obj.config.bc, t);
            end
        end
    end
end