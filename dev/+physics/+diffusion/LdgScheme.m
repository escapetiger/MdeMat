classdef LdgScheme < physics.Scheme
    % LDGSCHEME Local DG scheme for diffusion problems. 

    properties
        operator % Diffusion operator
        L % Linear term
        S % Source term
    end

    methods
        function state = initialize(obj, state)
            % INITIALIZE Initialize the scheme.

            % Reset timer
            obj.timer.reset();

            % Initial condition
            state.dofs.U = state.xDisc.space.project(obj.config.ic);

            % Diffusion operator
            if isempty(obj.config.bc)
                bcType = 'periodic';
            else
                bcType = 'dirichlet';
            end
            obj.operator = approx.assembly.DiffusionOperator( ...
                state.xDisc.space, state.xDisc.operator, bcType, ...
                obj.config.xAuxFluxType, ...
                obj.config.xPrmFluxType, ...
                obj.config.xAuxBoundaryJumpType);

            % Linear term
            obj.L = obj.operator.linear(obj.config.diffusion);

            % Source term
            if ~isempty(obj.config.bc)
                c = obj.config.diffusion;
                if strcmpi(obj.config.xAuxBoundaryJumpType, 'implicit')
                    g = [];
                else
                    g = @(X) state.xDisc.space.evaluate([], X, state.dofs.U);
                end
                obj.S = @(~, t) obj.operator.linearBc(c, obj.config.bc, g, t);
            else
                obj.S = [];
            end

            % Reset time discretization
            obj.tDisc.reset();

            % Reset visualizer
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
            % PRESTEP Preprocess step.

            % Set adaptive time step based on CFL condition.
            h = state.xDisc.mesh.measure;
            if strcmpi(obj.config.xAuxBoundaryJumpType, 'explicit')
                obj.tDisc.timeline.setTimeStep(h, obj.config.cfl, 2);
            else
                obj.tDisc.timeline.setTimeStep(h, obj.config.cfl);
            end
            if isa(obj.tDisc, 'approx.odeint.BdfIntegrator')
                obj.tDisc.setCoefficients();
            end
        end

        function state = step(obj, state)
            % STEP Perform one time step integration.
            %
            % Syntax:
            %   step(obj, state)
            %
            % Inputs:
            %   state - Solution state to advance

            obj.tDisc.update(state.dofs.U);
            state.dofs.U = obj.tDisc.step(obj.L, obj.S, []);
        end

        function state = postStep(obj, state)
            % POSTSTEP Postprocess step.

            % Update visualization if enabled
            if ~isempty(obj.visualizer) && obj.visualizer.isEnabled
                obj.timer.start(2, '[M] Visualization.');
                obj.visualizer.plot(state.xDisc.space, state.dofs, obj.tDisc);
                obj.timer.stop(2);
            end

            obj.tDisc.timeline.advance();
        end
    
        function state = finalize(obj, varargin)
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
end
