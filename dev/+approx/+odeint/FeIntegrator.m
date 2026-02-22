classdef FeIntegrator < approx.odeint.OdeIntegrator
    % FEINTEGRATOR Forward Euler integrator.
    %
    %   FeIntegrator implements the Forward Euler method, a first-order
    %   explicit time integration scheme for solving ordinary differential
    %   equations of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + F(u) + S(t),
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator,
    %   \f$F(u)\f$ is a nonlinear operator, and \f$S(t)\f$ is a source
    %   term. Any of these terms can be empty. For mass term, it means it
    %   reduces to an identity. For linear, nonlinear and source term, it
    %   means they vanish from the equation. This is the simplest explicit
    %   method, using only one evaluation per time step.
    %
    %   Forward Euler has a stability region of \f$|1 + z| \leq 1\f$, which
    %   restricts the time step for stability. For eigenvalue
    %   \f$\lambda\f$, the stability condition is \f$|1 + \lambda \Delta t|
    %   \leq 1\f$.
    %
    %   Forward Euler is conditionally stable and may require very small
    %   time steps for stiff problems. Use implicit methods for better
    %   stability properties with stiff systems.
    %
    % See also:
    %   approx.odeint.OdeIntegrator, approx.odeint.BeIntegrator,
    %   approx.odeint.ExrkIntegrator
    
    properties (Constant)
        Order = 1 % Accuracy order of the Forward Euler method
    end
    
    methods
        function obj = FeIntegrator(final)
            % FEINTEGRATOR Constructor for FeIntegrator.
            %
            %   obj = FeIntegrator(final) creates a Forward Euler integrator
            %   with the specified final time.
            
            arguments
                final {mustBeNumeric, mustBeNonnegative}
            end
            
            obj@approx.odeint.OdeIntegrator(1, 1, final);
        end

        function obj = reset(obj, options)
            % RESET Reset the FE integrator to initial state.
            %
            %   obj = reset(obj) calls the parent reset method and resets
            %   the underlying linear solver.

            arguments
                obj approx.odeint.FeIntegrator
                options.linOpts struct = struct() % Linear solver options
            end

            reset@approx.odeint.OdeIntegrator(obj);

            args = namedargs2cell(options.linOpts);
            obj.Solver = core.linalg.LinearSolver(args{:});
        end
        
        function obj = setLhs(obj, options)
            % SETLhs Set the left-hand side matrix for the Forward Euler
            % system.
            %
            %   obj = setLhs(obj, M=M) sets the left-hand side matrix to the
            %   mass operator.

            arguments
                obj approx.odeint.FeIntegrator
                options.M = [] % Mass operator
                options.dt = [] % Time step override
                options.t = [] % Time values for time-dependent operators
            end

            %< Evaluate mass operator for LHS
            Me = obj.evalMLhs(options.M, options.t);

            obj.Solver.setLhs(Me);
            obj.Solver.precompute();
        end
        
        function obj = setRhs(obj, options)
            % SETRhs Set the right-hand side vector for the Forward Euler
            % system.
            %
            %   obj = setRhs(obj, F=F, S=S, U=U, M=M) constructs the right-hand
            %   side vector using the explicit Forward Euler formula: M*U0
            %   + dt*(L*U0 + F + S).

            arguments
                obj approx.odeint.FeIntegrator
                options.F = [] % Nonlinear operator
                options.S = [] % Source term
                options.U = [] % Solution values (defaults to obj.U0{1})
                options.M = [] % Mass operator for RHS
                options.dt = [] % Time step override
                options.t = [] % Time values
            end

            % Use custom time step if provided
            if ~isempty(options.dt)
                dt = options.dt;
            else
                dt = obj.Timeline.StepSize;
            end

            % Use provided solution or current history
            if ~isempty(options.U)
                U0 = options.U;
            else
                U0 = obj.U0{1};
            end

            %< Evaluate RHS mass operator
            Me = obj.evalMRhs([], options.M, options.t);

            %< Add mass term contribution from previous solution
            if isempty(Me)
                Rhs = U0;
            elseif iscell(Me)
                Rhs = Me{2} * U0;
                if ~isempty(Me{1})
                    Rhs = Rhs + Me{1};
                end
            else
                Rhs = Me * U0;
            end

            %< Add nonlinear term contribution: dt*F(U0)
            if ~isempty(options.F)
                Fe = obj.evalOp(options.F, U0, options.t);
                Rhs = Rhs + dt * Fe;
            end

            %< Add source term contribution: dt*S(t0)
            if ~isempty(options.S)
                Se = obj.evalOp(options.S, U0, options.t);
                Rhs = Rhs + dt * Se;
            end

            obj.Solver.setRhs(Rhs);
        end
        
        function U = stage(obj, options)
            % STAGE Compute the single stage of Forward Euler method.
            %
            %   U = stage(obj, L=L, F=F, S=S, M=M) computes the solution at the
            %   next time step using the explicit Forward Euler formula.
            %   For Forward Euler, there is only one stage per time step.

            arguments
                obj approx.odeint.FeIntegrator
                options.L = [] % Linear operator
                options.F = [] % Nonlinear operator
                options.S = [] % Source term
                options.M = [] % Mass operator
                options.dt = [] % Time step override
                options.fn = [] % Custom step function
            end

            % Use custom step function if provided
            if ~isempty(options.fn)
                U = options.fn(obj, options);
                return;
            end

            % Use custom time step if provided
            if ~isempty(options.dt)
                dt = options.dt;
            else
                dt = obj.Timeline.StepSize;
            end

            t0 = obj.Timeline.Now;

            %< Update Lhs matrix if step size changed or custom dt provided
            if obj.Timeline.HasStepSizeChanged || ~isempty(options.dt)
                obj.setLhs(M=options.M, dt=dt, t=t0+dt);
            end

            %< Construct Rhs vector with explicit formula
            obj.setRhs(F=options.F, S=options.S, M=options.M, dt=dt, t=[t0+dt, t0]);

            %< Solve the linear system
            U = obj.Solver.solve();
        end
        
        function U = step(obj, options)
            % STEP Advance one time step using Forward Euler method.
            %
            %   U = step(obj, L=L, F=F, S=S, M=M) performs a complete time step
            %   using the explicit Forward Euler method. For Forward Euler,
            %   this is equivalent to computing a single stage.

            arguments
                obj approx.odeint.FeIntegrator
                options.L = [] % Linear operator
                options.F = [] % Nonlinear operator
                options.S = [] % Source term
                options.M = [] % Mass operator
                options.dt = [] % Time step override
                options.fn = [] % Custom step function
            end

            U = obj.stage(L=options.L, F=options.F, S=options.S, ...
                M=options.M, dt=options.dt, fn=options.fn);
        end
    end
end