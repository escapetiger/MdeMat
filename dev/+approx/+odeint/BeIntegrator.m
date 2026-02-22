classdef BeIntegrator < approx.odeint.OdeIntegrator
    % BEINTEGRATOR Backward Euler integrator.
    %
    %   BeIntegrator implements the Backward Euler method, a first-order
    %   implicit time integration scheme for solving ordinary differential
    %   equations of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + F(u) + S,
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator,
    %   \f$F(u)\f$ is a nonlinear operator, and \f$S\f$ is a source term.
    %   Any of these terms can be empty. For mass operator, it means it
    %   reduces to an identity. For linear, nonlinear and source terms,
    %   it means they vanish from the equation. The method is A-stable
    %   and suitable for stiff problems.
    %
    %   The Backward Euler method is the simplest implicit method,
    %   requiring only one function evaluation per time step but involving
    %   the solution of a linear system at each step.
    %
    % See also:
    %   approx.odeint.OdeIntegrator, approx.odeint.BdfIntegrator,
    %   approx.odeint.FeIntegrator

    properties (Constant)
        Order = 1 % Accuracy order of the Backward Euler method
    end

    methods
        function obj = BeIntegrator(final)
            % BEINTEGRATOR Constructor for BeIntegrator.
            %
            %   obj = BeIntegrator(final) creates a Backward Euler
            %   integrator with the specified final time.

            arguments
                final {mustBeNumeric, mustBeNonnegative}
            end

            obj@approx.odeint.OdeIntegrator(1, 1, final);
        end

        function obj = reset(obj, options)
            % RESET Reset the BE integrator to initial state.
            %
            %   obj = reset(obj) calls the parent reset method and resets
            %   the underlying linear solver.

            arguments
                obj approx.odeint.BeIntegrator
                options.linOpts struct = struct() % Linear solver options
            end

            reset@approx.odeint.OdeIntegrator(obj);

            args = namedargs2cell(options.linOpts);
            obj.Solver = core.linalg.LinearSolver(args{:});
        end

        function obj = setLhs(obj, options)
            % SETLhs Set the left-hand side matrix for the Backward Euler
            % system.
            %
            %   obj = setLhs(obj, L=L, M=M) constructs the left-hand side
            %   matrix as M - dt*L for the implicit Backward Euler system.
            %
            %   obj = setLhs(obj, L=L, M=M, P=P, A=A, B=B) constructs the
            %   left-hand side matrix for augmented Backward Euler systems:
            %
            %   \f[
            %     \begin{bmatrix}
            %       M - dt \cdot L & -dt \cdot P \\
            %       A & B
            %     \end{bmatrix}
            %   \f]

            arguments
                obj approx.odeint.BeIntegrator
                options.L = [] % Linear operator
                options.M = [] % Mass operator
                options.P = [] % Coupling operator (triggers augmented mode)
                options.A = [] % Constraint matrix (triggers augmented mode)
                options.B = [] % Constraint coupling matrix
                options.dt = [] % Time step override
                options.t = [] % Time values for time-dependent operators
            end

            % Use custom time step if provided
            if ~isempty(options.dt)
                dt = options.dt;
            else
                dt = obj.Timeline.StepSize;
            end

            % Detect augmented mode
            isAugmented = ~isempty(options.P) || ~isempty(options.A) || ~isempty(options.B);

            %< Evaluate mass operator for LHS
            Me = obj.evalMLhs(options.M, options.t);

            %< Evaluate linear operator if provided
            if ~isempty(options.L)
                Le = obj.evalOp(options.L, [], options.t);
            else
                Le = [];
            end

            if isAugmented
                %< Construct augmented block matrix
                Lhs = cell(2, 2);

                %< Block (1,1): M - dt*L
                Lhs{1, 1} = Me;
                if ~isempty(Le)
                    Lhs{1, 1} = Lhs{1, 1} - dt * Le;
                end

                %< Block (1,2): -dt*P
                if ~isempty(options.P)
                    Pe = obj.evalOp(options.P, [], options.t);
                    Lhs{1, 2} = -dt * Pe;
                else
                    Lhs{1, 2} = [];
                end

                %< Block (2,1): A
                Lhs{2, 1} = options.A;

                %< Block (2,2): B
                Lhs{2, 2} = options.B;

                Lhs = core.linalg.block(Lhs);
            else
                %< Regular: M - dt*L
                Lhs = Me;
                if ~isempty(Le)
                    Lhs = Lhs - dt * Le;
                end
            end

            %< Set LHS in solver and precompute factorization
            obj.Solver.setLhs(Lhs);
            obj.Solver.precompute();
        end

        function obj = setRhs(obj, options)
            % SETRhs Set the right-hand side vector for the Backward Euler
            % system.
            %
            %   obj = setRhs(obj, F=F, S=S, U=U) constructs the right-hand side
            %   vector using the mass operator, nonlinear term, and source term.
            %
            %   obj = setRhs(obj, F=F, S=S, U=U, M=M, C=C) constructs the
            %   right-hand side vector for augmented Backward Euler systems:
            %
            %   \f[
            %     \begin{bmatrix}
            %       \text{regular RHS} \\
            %       C
            %     \end{bmatrix}
            %   \f]

            arguments
                obj approx.odeint.BeIntegrator
                options.F = [] % Nonlinear operator
                options.S = [] % Source term
                options.U = [] % Solution values (defaults to obj.U0{1})
                options.M = [] % Mass operator for RHS
                options.C = [] % Constraint source term (triggers augmented mode)
                options.dt = [] % Time step override
                options.t = [] % Time values
            end

            % Use custom time step if provided
            if ~isempty(options.dt)
                dt = options.dt;
            else
                dt = obj.Timeline.StepSize;
            end

            % Detect augmented mode
            isAugmented = ~isempty(options.C);

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
                Rhs1 = U0;
            elseif iscell(Me)
                if isAugmented
                    Rhs1 = Me{2} * U0(1:size(Me{2}, 2));
                    if ~isempty(Me{1})
                        Rhs1 = Rhs1 + Me{1}(1:size(Me{2}, 2));
                    end
                else
                    Rhs1 = Me{2} * U0;
                    if ~isempty(Me{1})
                        Rhs1 = Rhs1 + Me{1};
                    end
                end
            else
                if isAugmented
                    Rhs1 = Me * U0(1:size(Me, 2));
                else
                    Rhs1 = Me * U0;
                end
            end

            %< Add nonlinear term contribution
            if ~isempty(options.F)
                Fe = obj.evalOp(options.F, U0, options.t);
                Rhs1 = Rhs1 + dt * Fe;
            end

            %< Add source term contribution
            if ~isempty(options.S)
                Se = obj.evalOp(options.S, U0, options.t);
                Rhs1 = Rhs1 + dt * Se;
            end

            if isAugmented
                %< Construct augmented RHS vector
                Rhs = cell(2, 1);
                Rhs{1} = Rhs1;
                Rhs{2} = options.C;
                Rhs = core.linalg.block(Rhs);
            else
                %< Regular RHS
                Rhs = Rhs1;
            end

            obj.Solver.setRhs(Rhs);
        end

        function U = stage(obj, options)
            % STAGE Compute the single stage of the Backward Euler method.
            %
            %   U = stage(obj, L=L, F=F, S=S, M=M) computes the solution at the
            %   current time step by solving the implicit Backward Euler
            %   system. For Backward Euler, there is only one stage per
            %   time step.
            %
            %   U = stage(obj, L=L, F=F, S=S, M=M, P=P, A=A, B=B, C=C) computes
            %   the solution for augmented Backward Euler systems with constraints.

            arguments
                obj approx.odeint.BeIntegrator
                options.L = [] % Linear operator
                options.F = [] % Nonlinear operator
                options.S = [] % Source term
                options.M = [] % Mass operator
                options.P = [] % Coupling operator (triggers augmented mode)
                options.A = [] % Constraint matrix (triggers augmented mode)
                options.B = [] % Constraint coupling matrix
                options.C = [] % Constraint source term
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
            t1 = t0 + dt;

            %< Update Lhs matrix if step size changed or if custom dt provided
            if obj.Timeline.HasStepSizeChanged || ~isempty(options.dt)
                obj.setLhs(L=options.L, M=options.M, P=options.P, A=options.A, B=options.B, dt=dt, t=t1);
            end

            %< Construct Rhs vector
            obj.setRhs(F=options.F, S=options.S, M=options.M, C=options.C, dt=dt, t=[t1, t0]);

            %< Solve the linear system
            U = obj.Solver.solve();
        end

        function U = step(obj, options)
            % STEP Advance one time step using the Backward Euler method.
            %
            %   U = step(obj, L=L, F=F, S=S, M=M) performs a complete time
            %   step by solving the implicit Backward Euler system. For
            %   Backward Euler, this is equivalent to computing a single
            %   stage.
            %
            %   U = step(obj, L=L, F=F, S=S, M=M, P=P, A=A, B=B, C=C)
            %   performs a complete time step for augmented Backward Euler
            %   systems.

            arguments
                obj approx.odeint.BeIntegrator
                options.L = [] % Linear operator
                options.F = [] % Nonlinear operator
                options.S = [] % Source term
                options.M = [] % Mass operator
                options.P = [] % Coupling operator (triggers augmented mode)
                options.A = [] % Constraint matrix (triggers augmented mode)
                options.B = [] % Constraint coupling matrix
                options.C = [] % Constraint source term
                options.dt = [] % Time step override
                options.fn = [] % Custom step function
            end

            U = obj.stage(L=options.L, F=options.F, S=options.S, M=options.M, ...
                P=options.P, A=options.A, B=options.B, C=options.C, ...
                dt=options.dt, fn=options.fn);
        end
    end
end
