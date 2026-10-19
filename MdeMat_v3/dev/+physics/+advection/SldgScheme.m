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
    % Notes:
    %   The semi-Lagrangian approach allows for larger time steps compared
    %   to traditional Eulerian methods, making it computationally efficient
    %   for advection-dominated flows.
    %
    % See Also:
    %   physics.Scheme, physics.advection.UwdgScheme,
    %   approx.assembly.SemiLagrangianUniformGridAssembly

    properties
        operator % Semi-Lagrangian transport operator
        M % Dynamic mass
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

            % Semi-Lagrangian advection operator
            if isempty(obj.config.bc)
                bcType = 'periodic';
            else
                bcType = 'dirichlet';
            end
            obj.operator = approx.assembly.SemiLagrangianConvectionOperator( ...
                state.xDisc.space, state.xDisc.slOperator, bcType);
            
            % Mass term
            obj.setMassTerm(state);

            % Initial condition
            state.dofs.U = state.xDisc.space.project(obj.config.ic);

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

            h = state.xDisc.mesh.measure;
            obj.tDisc.setTimeStep(h, obj.config.cfl, []);
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
            %   obj - The SldgScheme object
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
        function obj = setMassTerm(obj, state)
            obj.M = @(tBegin, tEnd) computeM(tBegin, tEnd);
        
            function varargout = computeM(tBegin, tEnd)
                ht = tBegin - tEnd;
                hx = cell2mat(state.xDisc.mesh.spacings);
                state.xDisc.clipper.setVolumePieceData(ht, hx);
                state.xDisc.slOperator.setVolumePieceData();
                [varargout{1:nargout}] = obj.operator.linear(tBegin, tEnd, obj.config.bc);
            end
        end
    end
end