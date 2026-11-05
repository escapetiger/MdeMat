classdef SeparableOdeIntegrator < core.linalg.Separable
    % SEPARABLEODELINTEGRATOR Base class for separable ODE integrators.
    %
    %   SeparableOdeIntegrator extends the Separable class to provide
    %   specialized functionality for ODE integration methods that
    %   decompose the evolution operator into separable factors. Each
    %   factor corresponds to a sub-integrator that handles one part
    %   of the split equation.
    %
    %   The mathematical foundation is based on decomposing:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = \sum_{i=1}^{n} (L_i u + F_i(u) + S_i)
    %   \f]
    %
    %   into n separable sub-problems, each solved by an appropriate
    %   factor integrator stored in the separable representation.
    %
    % Examples:
    %   % Create factor integrators
    %   integrators = {approx.odeint.BeIntegrator(1.0), ...
    %                        approx.odeint.FeIntegrator(1.0)};
    %   
    %   % Create separable integrator (abstract - use concrete subclass)
    %   integrator = approx.odeint.LieTrotterIntegrator(integrators);
    %   
    %   % Access factor integrators using inherited methods
    %   diffusionIntegrator = integrator.getFactor(1);
    %   reactionIntegrator = integrator.getFactor(2);
    %
    % See also:
    %   core.linalg.Separable, approx.odeint.LieTrotterIntegrator,
    %   approx.odeint.StrangIntegrator

    properties (Constant)
        IMPLICIT_TYPES = {'approx.odeint.BeIntegrator', ...
                         'approx.odeint.DirkIntegrator', ...
                         'approx.odeint.BdfIntegrator'}
    end

    properties
        timeline % Timeline object shared across all factor integrators
    end

    properties (Dependent)
        U0
    end

    methods
        function obj = SeparableOdeIntegrator(varargin)
            % SEPARABLEODELINTEGRATOR Constructor for SeparableOdeIntegrator.
            %
            %   obj = SeparableOdeIntegrator() creates a separable ODE
            %   integrator with empty factors.
            %
            %   obj = SeparableOdeIntegrator(integrators) creates a
            %   separable ODE integrator with the specified factor integrators.
            %
            % Inputs:
            %   integrators - Cell array of integrator objects (optional)
            %
            % Outputs:
            %   obj - Constructed SeparableOdeIntegrator object

            obj@core.linalg.Separable(varargin{:});
            
            if nargin == 0 || obj.nFactors == 0
                obj.timeline = [];
            else
                final = obj.getFactor(1).timeline.final;
                obj.timeline = approx.mesh.DynamicTimeline(final, 1);
            end
        end

        function U0 = get.U0(obj)
            U0 = obj.getFactor(1).U0;
        end

        function obj = setTimeStep(obj, varargin)
            % SETTIMESTEP Set time step for the main timeline and all factors.
            %
            %   obj = setTimeStep(obj, h, C) sets time step as dt = C*h
            %   for dynamic timelines.
            %
            %   obj = setTimeStep(obj, h, C, p) sets time step as dt = C*h^p
            %   for dynamic timelines.
            %
            % Inputs:
            %   obj - The SeparableOdeIntegrator object
            %   varargin - Arguments passed to timeline.setTimeStep
            %
            % Outputs:
            %   obj - The SeparableOdeIntegrator object

            if ~isempty(obj.timeline)
                obj.timeline.setTimeStep(varargin{:});
                
                for i = 1:obj.nFactors
                    factorTimeline = obj.getFactor(i).timeline;
                    if isa(factorTimeline, 'approx.mesh.DynamicTimeline')
                        factorTimeline.setTimeStep(varargin{:});
                    end
                end
            end
        end

        function obj = reset(obj, i)
            % RESET Reset factor integrators to initial state.
            %
            %   obj = reset(obj) resets all factor integrators to
            %   their initial states.
            %
            %   obj = reset(obj, i) resets the specified factor
            %   integrator to its initial state.
            %
            % Inputs:
            %   obj - The SeparableOdeIntegrator object
            %   i - Factor index (positive integer, optional)
            %
            % Outputs:
            %   obj - The SeparableOdeIntegrator object

            if nargin < 2
                i = [];
            end

            if isempty(i)
                for i = 1:obj.nFactors
                    obj.getFactor(i).reset();
                end
                if ~isempty(obj.timeline)
                    obj.timeline.reset();
                end
            else
                obj.getFactor(i).reset();
            end
        end

        function obj = update(obj, U, i)
            % UPDATE Update solution history for factor integrators.
            %
            %   obj = update(obj, U) updates the solution history for all
            %   factor integrators with the new solution.
            %   
            %   obj = update(obj, U, i) updates the solution
            %   history for the specified factor integrator with the new
            %   solution.
            %
            % Inputs:
            %   obj - The SeparableOdeIntegrator object
            %   U - The new solution to add to history
            %   i - Factor index (positive integer, optional)
            %
            % Outputs:
            %   obj - The SeparableOdeIntegrator object

            if nargin < 3
                i = [];
            end

            if ~isempty(i)
                obj.validateFactorIndex(i);
            end

            if isempty(i)
                for i = 1:obj.nFactors
                    obj.getFactor(i).update(U);
                end
            else
                obj.getFactor(i).update(U);
            end
        end
        
        function obj = advance(obj, i)
            % ADVANCE Advance timeline for factor integrators.
            %
            %   obj = advance(obj) advances the timeline for the main
            %   integrator and all factor integrators. This should be
            %   called once per complete splitting step.
            %
            %   obj = advance(obj, i) advances the timeline for the
            %   specified factor integrator only. This is used internally
            %   during sub-steps.
            %
            % Inputs:
            %   obj - The SeparableOdeIntegrator object
            %   i - Factor index (positive integer, optional)
            %
            % Outputs:
            %   obj - The SeparableOdeIntegrator object

            if nargin < 2
                i = [];
            end

            if isempty(i)
                if ~isempty(obj.timeline)
                    obj.timeline.advance();
                end
                for i = 1:obj.nFactors
                    obj.getFactor(i).advance();
                end
            else
                obj.getFactor(i).advance();
            end
        end

        function U = applyFactorStep(obj, i, h, L, F, S, M)
            % APPLYFACTORSTEP Apply a single factor step.
            %
            %   U = applyFactorStep(obj, i, h, L, F, S, M)
            %   applies one step of the specified factor integrator with
            %   the given step size and operators. The factor's timeline
            %   is temporarily modified and then restored.
            %
            % Inputs:
            %   obj - The SeparableOdeIntegrator object
            %   i - Index of the factor integrator to apply
            %   h - Step size for this factor step
            %   L - Linear operator for this factor
            %   F - Nonlinear operator for this factor
            %   S - Source term for this factor
            %   M - Mass operator (shared across factors)
            %
            % Outputs:
            %   U - Updated solution vector

            integrator = obj.getFactor(i);

            h0 = integrator.timeline.dt;
            
            if isa(integrator.timeline, 'approx.mesh.DynamicTimeline')
                integrator.timeline.h(1) = h;
                integrator.timeline.next = integrator.timeline.now + h;
            end

            if any(cellfun(@(cls) isa(obj.getFactor(i), cls), obj.IMPLICIT_TYPES))
                U = integrator.step(L, S, M);
            else
                U = integrator.step(L, F, S, M);
            end
            
            if isa(integrator.timeline, 'approx.mesh.DynamicTimeline')
                integrator.timeline.h(1) = h0;
                integrator.timeline.next = integrator.timeline.now + h0;
            end
        end
    end

    methods (Abstract)
        % STEP Advance one time step using the specific splitting method.
        U = step(obj, L, F, S, M)
    end
end