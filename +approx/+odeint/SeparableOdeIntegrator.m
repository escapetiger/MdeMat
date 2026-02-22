classdef SeparableOdeIntegrator < core.linalg.Separable
    % SEPARABLEODEINTEGRATOR Base class for separable ODE integrators.
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
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = \sum_{i=1}^{n} (L_i u +
    %     F_i(u) + S_i)
    %   \f]
    %
    %   into n separable sub-problems, each solved by an appropriate
    %   factor integrator stored in the separable representation.
    %
    % See also:
    %   core.linalg.Separable, approx.odeint.LieTrotterIntegrator,
    %   approx.odeint.StrangIntegrator
    
    properties (Constant)
        ImplicitTypes = {'approx.odeint.BeIntegrator', ...
            'approx.odeint.DirkIntegrator', ...
            'approx.odeint.BdfIntegrator'}
    end
    
    properties
        Timeline % Timeline object shared across all factor integrators
    end
    
    properties (Dependent)
        U0 % Solution history from the first factor integrator
    end
    
    methods
        function obj = SeparableOdeIntegrator(options)
            % SEPARABLEODEINTEGRATOR Constructor for
            % SeparableOdeIntegrator.
            %
            %   obj = SeparableOdeIntegrator() creates a separable ODE
            %   integrator with empty factors.
            %
            %   obj = SeparableOdeIntegrator(factors=factors) creates a
            %   separable ODE integrator with the specified factor @a
            %   integrators.

            arguments
                options.factors {core.except.mustBeCellOrObject} = {}
            end
            
            obj@core.linalg.Separable(factors=options.factors);
            
            if obj.NFactors == 0
                obj.Timeline = [];
            else
                final = obj.getFactor(1).Timeline.Final;
                obj.Timeline = approx.mesh.DynamicTimeline(final, nSteps=1);
            end
        end
        
        function U0 = get.U0(obj)
            U0 = obj.getFactor(1).U0;
        end
        
        function obj = setTimeStep(obj, options)
            % SETTIMESTEP Set time step for the main Timeline and all
            % factors.
            %
            %   obj = setTimeStep(obj, h=h, C=C) sets time step as dt = C*h
            %   for dynamic Timelines.
            %
            %   obj = setTimeStep(obj, h=h, C=C, p=p) sets time step as
            %   StepSize = C*h^p for dynamic Timelines.
            %
            %   obj = setTimeStep(obj, dt=dt) directly sets the time step.

            arguments
                obj approx.odeint.SeparableOdeIntegrator
                options.dt {mustBeNonnegative} = []
                options.h {mustBePositive} = []
                options.C {mustBePositive} = 1
                options.p {mustBePositive} = 1
            end
            
            if ~isempty(obj.Timeline)
                if ~isempty(options.dt)
                    obj.Timeline.setTimeStep(dt=options.dt);
                    for i = 1:obj.NFactors
                        obj.getFactor(i).setTimeStep(dt=options.dt);
                    end
                else
                    obj.Timeline.setTimeStep(h=options.h, C=options.C, p=options.p);
                    for i = 1:obj.NFactors
                        obj.getFactor(i).setTimeStep(h=options.h, C=options.C, p=options.p);
                    end
                end
            end
        end
        
        function obj = reset(obj, options)
            % RESET Reset factor integrators to initial state.
            %
            %   obj = reset(obj) resets all factor integrators to their
            %   initial states.
            %
            %   obj = reset(obj, i=i) resets the specified factor
            %   integrator to its initial state.

            arguments
                obj
                options.i {mustBePositive, mustBeInteger} = []
                options.linOpts = struct()
            end
            
            if isempty(options.i)
                for i = 1:obj.NFactors
                    obj.getFactor(i).reset(linOpts=options.linOpts);
                end
                if ~isempty(obj.Timeline)
                    obj.Timeline.reset();
                end
            else
                obj.getFactor(options.i).reset(linOpts=options.linOpts);
            end
        end
        
        function obj = update(obj, U, options)
            % UPDATE Update solution history for factor integrators.
            %
            %   obj = update(obj, U) updates the solution history for all
            %   factor integrators with the new solution @a U.
            %
            %   obj = update(obj, U, i=i) updates the solution history for
            %   the specified factor integrator with the new solution @a U.

            arguments
                obj approx.odeint.SeparableOdeIntegrator
                U {mustBeNumeric}
                options.i {mustBePositive, mustBeInteger} = []
            end
            
            if ~isempty(options.i)
                obj.validateFactorIndex(options.i);
            end
            
            if isempty(options.i)
                for i = 1:obj.NFactors
                    obj.getFactor(i).update(U);
                end
            else
                obj.getFactor(options.i).update(U);
            end
        end
        
        function obj = add(obj, U, options)
            % ADD Add new solution for factor integrators.
            %
            %   obj = add(obj, U) adds the newest solution @a U for all
            %   factor integrators.
            %
            %   obj = add(obj, U, i=i) adds the newest solution @a U for
            %   the specified factor integrator.

            arguments
                obj approx.odeint.SeparableOdeIntegrator
                U {mustBeNumeric}
                options.i {mustBePositive, mustBeInteger} = []
            end
            
            if ~isempty(options.i)
                obj.validateFactorIndex(options.i);
            end
            
            if isempty(options.i)
                for i = 1:obj.NFactors
                    obj.getFactor(i).add(U);
                end
            else
                obj.getFactor(options.i).add(U);
            end
        end
        
        function obj = advance(obj, options)
            % ADVANCE Advance Timeline for factor integrators.
            %
            %   obj = advance(obj) advances the Timeline for the main
            %   integrator and all factor integrators. This should be
            %   called once per complete splitting step.
            %
            %   obj = advance(obj, i=i) advances the Timeline for the
            %   specified factor integrator only. This is used internally
            %   during sub-steps.

            arguments
                obj approx.odeint.SeparableOdeIntegrator
                options.i {mustBePositive, mustBeInteger} = []
            end
            
            if isempty(options.i)
                if ~isempty(obj.Timeline)
                    obj.Timeline.advance();
                end
                for i = 1:obj.NFactors
                    obj.getFactor(i).advance();
                end
            else
                obj.getFactor(options.i).advance();
            end
        end
        
        function U = factorStep(obj, i, options)
            % FACTORSTEP Apply a single factor step.
            %
            %   U = factorStep(obj, i, L=L, F=F, S=S, M=M) performs the i-th regular 
            %   factor time step with the given operators.
            %
            %   U = factorStep(obj, i, L=L, S=S, M=M, P=P, A=A, B=B, C=C) performs
            %   the i-th augmented factor time step with the given operators.
            %
            %   U = factorStep(obj, i, dt=dt, fn=fn) performs the i-th custom factor 
            %   time step with the given operators.

            arguments
                obj approx.odeint.SeparableOdeIntegrator
                i {mustBePositive, mustBeInteger}
                options.L = [] % Linear operator
                options.F = [] % Nonlinear operator
                options.S = [] % Source term
                options.M = [] % Mass operator
                options.P = [] % Coupling operator
                options.A = [] % Constraint matrix
                options.B = [] % Constraint coupling matrix
                options.C = [] % Constraint source term
                options.dt = [] % Time step override
                options.fn = [] % Custom step function (function handle)
            end

            %< Use custom step function if provided
            if ~isempty(options.fn)
                U = options.fn(obj, options);
                return;
            end
            
            integrator = obj.getFactor(i);

            L = options.L;
            F = options.F;
            S = options.S;
            M = options.M;
            P = options.P;
            A = options.A;
            B = options.B;
            C = options.C;
            dt = options.dt;
            fn = options.fn;
            isAugmented = ~isempty(P) || ~isempty(A) || ~isempty(B) || ~isempty(C);

            %< Call step with new interface
            if any(cellfun(@(cls) isa(obj.getFactor(i), cls), obj.ImplicitTypes))
                if isAugmented
                    U = integrator.step(L=L, S=S, M=M, P=P, A=A, B=B, C=C, dt=dt, fn=fn);
                else
                    U = integrator.step(L=L, S=S, M=M, dt=dt, fn=fn);
                end
            else
                if isAugmented
                    U = integrator.step(L=L, F=F, S=S, M=M, P=P, A=A, B=B, C=C, dt=dt, fn=fn);
                else
                    U = integrator.step(L=L, F=F, S=S, M=M, dt=dt, fn=fn);
                end
            end
        end
    end

    methods (Static)
        function xi = getEntry(i, x)
            % GETENTRY Utility to extract the i-th entry from a cell array or return
            % the input directly if it's not a cell array.

            if isempty(x)
                xi = [];
            elseif iscell(x)
                xi = x{i};
            else
                xi = x;
            end
        end
    end

    methods (Abstract)
        % STEP Advance one time step using the specific splitting method.
        U = step(obj, options)
    end
end