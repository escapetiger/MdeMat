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
    %                  approx.odeint.FeIntegrator(1.0)};
    %   
    %   % Create separable integrator (abstract - use concrete subclass)
    %   obj = approx.odeint.LieTrotterIntegrator(1.0, integrators);
    %   
    %   % Access factor integrators using inherited methods
    %   diffusionIntegrator = obj.getFactor(1);
    %   reactionIntegrator = obj.getFactor(2);
    %
    % See also:
    %   core.linalg.Separable, approx.odeint.LieTrotterIntegrator,
    %   approx.odeint.StrangIntegrator

    properties (Constant)
        IMPLICIT_TYPES = {'approx.odeint.BeIntegrator', ...
            'approx.odeint.DirkIntegrator', ...
            'approx.odeint.BdfIntegrator'};
    end

    methods
        function obj = reset(obj, k)
            % RESET Reset all factor integrators to initial state.
            %
            %   obj = reset(obj) resets all factor integrators to
            %   their initial states.
            %
            %   obj = reset(obj, k) reset the specefied factor
            %   integrator to its initial state.
            %
            % Inputs:
            %   obj - The SeparableOdeIntegrator object
            %   k - Factor index (positive integer)
            %
            % Outputs:
            %   obj - The SeparableOdeIntegrator object

            if nargin < 2 || isempty(k)
                k = [];
            end

            if isempty(k)
                for i = 1:obj.nFactors
                    obj.getFactor(i).reset();
                end
            else
                obj.getFactor(k).reset();
            end
        end

        function obj = update(obj, U, k)
            % UPDATE Update solution history for all factor integrators.
            %
            %   obj = update(obj, U) updates the solution history for all
            %   factor integrators with the new solution.
            %   
            %   obj = reset(obj, U, k) updates the solution history for
            %   the specified factor integrator with the new solution.
            %
            % Inputs:
            %   obj - The SeparableOdeIntegrator object
            %   U - The new solution to add to history
            %   k - Factor index (positive integer)
            %
            % Outputs:
            %   obj - The SeparableOdeIntegrator object

            if nargin < 2 || isempty(k)
                k = [];
            end

            if isempty(k)
                for i = 1:obj.nFactors
                    obj.getFactor(i).update(U);
                end
            else
                obj.getFactor(k).update(U);
            end
        end

        function obj = advance(obj)
            % ADVANCE Advance all factor integrator timelines.
            %
            %   obj = advance(obj) advances the timeline for all factor
            %   integrators.
            %
            % Inputs:
            %   obj - The SeparableOdeIntegrator object
            %
            % Outputs:
            %   obj - The SeparableOdeIntegrator object

            obj.timeline.advance();

            for i = 1:obj.nFactors
                obj.getFactor(i).advance();
            end
        end

        function U = applyFactorStep(obj, k, h, L, F, S, M)
            % APPLYFACTORSTEP Apply a single factor step.
            %
            %   U = applyFactorStep(obj, k, h, L, F, S, M) applies one
            %   step of the specified factor integrator with the given step
            %   size and operators.
            %
            % Inputs:
            %   obj - The SeparableOdeIntegrator object
            %   k - Index of the factor integrator to apply
            %   f - Step size for this factor step
            %   L - Linear operator for this factor
            %   F - Nonlinear operator for this factor
            %   S - Source term for this factor
            %   M - Mass operator (shared across factors)
            %
            % Outputs:
            %   U - Updated solution vector

            integrator = obj.getFactor(k);

            h0 = integrator.timeline.dt;
            integrator.timeline.setTimeStep(h, 1, 1);

            if any(strcmp(class(obj), obj.IMPLICIT_TYPES))
                U = integrator.step(L, S, M);
            else
                U = integrator.step(L, F, S, M);
            end
            
            integrator.timeline.setTimeStep(h0, 1, 1);
        end
    end

    methods (Abstract)
        % STEP Advance one time step using the specific splitting method.
        U = step(obj, L, F, S, M)
    end
end