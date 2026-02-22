classdef SldgScheme < physics.Scheme
    % SLDGSCHEME Semi-Lagrangian Discontinuous Galerkin scheme.
    %
    %   SldgScheme implements a Semi-Lagrangian Discontinuous Galerkin
    %   (SLDG) method for solving advection equations. The method combines
    %   the stability advantages of Lagrangian methods with the flexibility
    %   of discontinuous Galerkin spatial discretizations.
    %
    %   The SLDG method traces characteristics backward in time and uses
    %   high-order interpolation to evaluate solutions at departure points,
    %   making it particularly effective for transport-dominated problems
    %   with large time steps and minimal numerical diffusion.
    %
    % Examples:
    %   % Create SLDG scheme with configuration
    %   config = struct('ic', @(x) sin(x), 'cfl', 0.8, 'bc', [], 'verbose', true);
    %   scheme = SldgScheme(config, timeIntegrator, visualizer, analyzer);
    %   
    %   % Full simulation workflow
    %   state = scheme.initialize(initialState);
    %   finalState = scheme.run(state);
    %   scheme.finalize(finalState);
    %
    % Notes:
    %   The semi-Lagrangian approach allows for larger time steps compared
    %   to traditional Eulerian methods, making it computationally efficient
    %   for advection-dominated flows.
    %
    % See Also:
    %   physics.Scheme, physics.advection.UwdgScheme,
    %   fem.assembly.SemiLagrangianUniformGridAssembly

    properties
        assembly % Transport assembly for semi-Lagrangian operations
        operator % Transport operator function handle
    end

    methods
        function state = initialize(obj, state)
            % INITIALIZE Initialize the SLDG scheme.
            %
            %   state = initialize(obj, state) sets up the semi-Lagrangian
            %   discontinuous Galerkin scheme by configuring finite element
            %   data, assembling the transport operator, projecting initial
            %   conditions, and preparing visualization components.
            %
            % Inputs:
            %   obj - The SldgScheme object
            %   state - Initial solution state structure
            %
            % Outputs:
            %   state - Initialized state ready for time integration

            % Reset timer
            obj.timer.reset();

            % Local data
            state.disc.x.fe.setVolumeData(1);
            state.disc.x.fe.setFluxData(0);
            state.disc.x.op.setVolumeData();
            state.disc.x.op.setFluxData();

            % Semi-Lagrangian advection operator
            obj.assembly = fem.assembly.SemiLagrangianUniformGridAssembly( ...
                state.disc.x.fe, state.disc.x.mesh, state.disc.x.slop, ...
                ~isempty(obj.config.bc));
            obj.operator = @(tBegin, tEnd) obj.transportFn(state, tBegin, tEnd);

            % Initial condition
            state.dofs.U = state.disc.x.space.project(obj.config.ic);

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
            %   obj - The SldgScheme object
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
            % STEP Perform one semi-Lagrangian time step.
            %
            %   state = step(obj, state) advances the solution by one time
            %   step using the semi-Lagrangian transport operator. The
            %   method updates the solution history and applies the
            %   characteristic-based transport.
            %
            % Inputs:
            %   obj - The SldgScheme object
            %   state - Current solution state
            %
            % Outputs:
            %   state - Updated solution state after one time step

            obj.tDisc.update(state.dofs.U);
            M = obj.operator;
            state.dofs.U = obj.tDisc.step([], [], M);
        end

        function state = postStep(obj, state)
            % POSTSTEP Postprocess operations after time step.
            %
            %   state = postStep(obj, state) performs postprocessing
            %   including visualization updates and timeline advancement
            %   for the next integration step.
            %
            % Inputs:
            %   obj - The SldgScheme object
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
            %   [state, coarse] = finalize(obj, state, coarse) performs
            %   relative error analysis comparing fine and coarse grid
            %   solutions for convergence studies.
            %
            % Inputs:
            %   obj - The SldgScheme object
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

        function varargout = transportFn(obj, state, tBegin, tEnd)
            % TRANSPORTFN Compute semi-Lagrangian transport operator.
            %
            %   M = transportFn(obj, state, tBegin, tEnd) computes the
            %   transport operator for the semi-Lagrangian method by
            %   setting up characteristic line parameters and assembling
            %   the transport matrix.
            %
            % Inputs:
            %   obj - The SldgScheme object
            %   state - Current solution state
            %   tBegin - Beginning time for transport step
            %   tEnd - End time for transport step
            %
            % Outputs:
            %   varargout - Transport operator matrix and related data

            ht = tBegin - tEnd;
            hx = cell2mat(state.disc.x.mesh.spacings);
            state.disc.x.clp.setVolumePieceData(ht, hx);
            state.disc.x.slop.setVolumePieceData();
            [varargout{1:nargout}] = obj.assembly.transport(tBegin, tEnd, obj.config.bc);
        end
    end
end