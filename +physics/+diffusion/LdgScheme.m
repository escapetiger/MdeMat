classdef LdgScheme < physics.Scheme
    % LDGSCHEME

    properties 
        assembly % Diffusion assembly
        operator % Diffusion operator
    end

    methods
        function state = initialize(obj, state)
            % INITIALIZE Initialize the scheme.

            % Reset timer
            obj.timer.reset();

            % Local data
            state.disc.x.fe.setVolumeData(1);
            state.disc.x.fe.setFluxData(0);
            state.disc.x.op.setVolumeData();
            state.disc.x.op.setFluxData();

            % Initial condition
            state.dofs.U = state.disc.x.space.project(obj.config.ic);

            % Diffusion operator
            obj.assembly = fem.assembly.DiffusionGridAssembly( ...
                state.disc.x.fe, state.disc.x.mesh, state.disc.x.op, ...
                ~isempty(obj.config.bc), obj.config.xDiffusionFlux(1), ...
                obj.config.xDiffusionFlux(2), obj.config.xDiffusionBoundaryJump);
            obj.operator = obj.assembly.hessian(obj.config.diffusion);

            % Reset time discretization
            obj.tDisc.reset();

            % Reset visualizer
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

            % Update visualization if enabled
            if obj.visualizer.isEnabled
                obj.timer.start(2, '[M] Visualization.');
                obj.visualizer.plot(state.disc.x.space, state.dofs, obj.tDisc);
                obj.timer.stop(2);
            end
        end

        function state = preStep(obj, state)
            % PRESTEP Preprocess step.

            % Set adaptive time step based on CFL condition.
            h = state.disc.x.mesh.measure;
            obj.tDisc.timeline.setTimeStep(h, obj.config.cfl);
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

            % Update boundary condition
            if ~isempty(obj.config.bc)
                c = obj.config.diffusion;
                if obj.config.xDiffusionBoundaryJump == 1
                    g = [];
                else
                    g = @(X) state.disc.x.space.evaluate(X, state.dofs.U);
                end
                S = @(t) obj.assembly.hessianBc(c, obj.config.bc, g, t);
            else
                S = [];
            end

            obj.tDisc.update(state.dofs.U);
            L = obj.operator;
            state.dofs.U = obj.tDisc.step(L, S, []);
        end

        function state = postStep(obj, state)
            % POSTSTEP Postprocess step.

            % Update visualization if enabled
            if obj.visualizer.isEnabled
                obj.timer.start(2, '[M] Visualization.');
                obj.visualizer.plot(state.disc.x.space, state.dofs, obj.tDisc);
                obj.timer.stop(2);
            end

            obj.tDisc.timeline.advance();
        end
    
        function state = finalize(obj, varargin)
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
