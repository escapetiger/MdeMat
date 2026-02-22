classdef DirkIntegrator < approx.odeint.OdeIntegrator
    % DIRKINTEGRATOR Diagonally Implicit Runge-Kutta integrator.
    %
    %   DirkIntegrator implements Diagonally Implicit Runge-Kutta (DIRK)
    %   methods for solving ordinary differential equations of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + F(u) + S,
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator,
    %   \f$F(u)\f$ is a nonlinear operator, and \f$S\f$ is a source term.
    %   Any of these terms can be empty. For mass operator, it means it
    %   reduces to an identity. For linear, nonlinear and source terms,
    %   it means they vanish from the equation.
    %
    %   DIRK methods are characterized by a lower triangular Butcher
    %   tableau with identical diagonal elements, allowing each stage to be
    %   solved independently. This makes them more efficient than fully
    %   implicit Runge-Kutta methods while maintaining good stability
    %   properties.
    %
    % See also:
    %   approx.odeint.OdeIntegrator, approx.odeint.BeIntegrator,
    %   approx.odeint.BdfIntegrator

    properties
        A % Butcher tableau coefficient matrix (lower triangular)
        b % Stage weights vector for final update
        c % Stage time coefficients vector
        K % Cell array of derivative approximations for each stage
    end

    methods
        function obj = DirkIntegrator(nStages, final)
            % DIRKINTEGRATOR Constructor for DirkIntegrator.
            %
            %   obj = DirkIntegrator(nStages, final) creates a DIRK integrator
            %   with the specified number of stages and final time.

            arguments
                nStages {mustBePositive, mustBeInteger}
                final {mustBeNumeric, mustBeNonnegative}
            end

            obj@approx.odeint.OdeIntegrator(1, nStages, final);
        end

        function obj = reset(obj, options)
            % RESET Reset the DIRK integrator to initial state.
            %
            %   obj = reset(obj) clears the derivative approximation
            %   storage and calls the parent reset method.

            arguments
                obj approx.odeint.DirkIntegrator
                options.linOpts struct = struct() % Linear solver options
            end

            reset@approx.odeint.OdeIntegrator(obj);
            obj.K = cell(obj.NStages, 1);

            args = namedargs2cell(options.linOpts);
            obj.Solver = cell(1, obj.NStages+1);
            for i = 1 : obj.NStages+1
                obj.Solver{i} = core.linalg.LinearSolver(args{:});
            end
        end

        function obj = setLhs(obj, i, options)
            % SETLhs Set the left-hand side matrix for stage i.
            %
            %   obj = setLhs(obj, i=i, L=L, M=M) constructs the left-hand side
            %   matrix as M - dt*A(i,i)*L for stage i of the DIRK system.
            %
            %   obj = setLhs(obj, i=i, L=L, M=M, P=P, A=A, B=B) constructs the
            %   left-hand side matrix for augmented DIRK systems:
            %
            %   \f[
            %     \begin{bmatrix}
            %       M - dt \cdot A(i,i) \cdot L & -dt \cdot A(i,i) \cdot P \\
            %       A & B
            %     \end{bmatrix}
            %   \f]

            arguments
                obj approx.odeint.DirkIntegrator
                i {mustBePositive, mustBeInteger} % Stage index
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

                %< Block (1,1): M - dt*A(i,i)*L
                Lhs{1, 1} = Me;
                if ~isempty(Le)
                    Lhs{1, 1} = Lhs{1, 1} - dt * obj.A(i, i) * Le;
                end

                %< Block (1,2): -dt*A(i,i)*P
                if ~isempty(options.P)
                    Pe = obj.evalOp(options.P, [], options.t);
                    Lhs{1, 2} = -dt * obj.A(i, i) * Pe;
                else
                    Lhs{1, 2} = [];
                end

                %< Block (2,1): A
                Lhs{2, 1} = options.A;

                %< Block (2,2): B
                Lhs{2, 2} = options.B;

                Lhs = core.linalg.block(Lhs);
            else
                %< Regular: M - dt*A(i,i)*L
                Lhs = Me;
                if ~isempty(Le)
                    Lhs = Lhs - dt * obj.A(i, i) * Le;
                end
            end

            obj.Solver{i}.setLhs(Lhs);
            obj.Solver{i}.precompute();
        end

        function obj = setRhs(obj, i, options)
            % SETRhs Set the right-hand side vector for stage i.
            %
            %   obj = setRhs(obj, i=i, F=F, S=S, U=U, M=M) constructs the
            %   right-hand side vector using mass operator, accumulated
            %   previous stages, nonlinear term, and source term.
            %
            %   obj = setRhs(obj, i=i, F=F, S=S, U=U, M=M, C=C) constructs
            %   the right-hand side vector for augmented DIRK systems with
            %   constraints:
            %
            %   \f[
            %     \begin{bmatrix}
            %       \text{regular RHS} \\
            %       C
            %     \end{bmatrix}
            %   \f]

            arguments
                obj approx.odeint.DirkIntegrator
                i {mustBePositive, mustBeInteger} % Stage index
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
            Me = obj.evalMRhs(i, options.M, options.t);

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

            %< Add accumulated term from previous stages
            for j = 1:(i - 1)
                Rhs1 = Rhs1 + dt * obj.A(i, j) * obj.K{j};
            end

            %< Add nonlinear term contribution
            if ~isempty(options.F)
                Fe = obj.evalOp(options.F, U0, options.t);
                Rhs1 = Rhs1 + dt * obj.A(i, i) * Fe;
            end

            %< Add source term contribution with diagonal coefficient
            if ~isempty(options.S)
                Se = obj.evalOp(options.S, U0, options.t);
                Rhs1 = Rhs1 + dt * obj.A(i, i) * Se;
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

            obj.Solver{i}.setRhs(Rhs);
        end

        function U = stage(obj, i, options)
            % STAGE Compute stage i of the DIRK method.
            %
            %   U = stage(obj, i=i, L=L, F=F, S=S, M=M) computes the solution at stage i
            %   by solving the implicit system and computing the derivative
            %   approximation for use in subsequent stages.
            %
            %   U = stage(obj, i=i, L=L, F=F, S=S, M=M, P=P, A=A, B=B, C=C) computes
            %   the solution for augmented DIRK systems with constraints.

            arguments
                obj approx.odeint.DirkIntegrator
                i {mustBePositive, mustBeInteger} % Stage number
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
                U = options.fn(obj, i, options);
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
                obj.setLhs(i, L=options.L, M=options.M, P=options.P, A=options.A, B=options.B, dt=dt, t=t0+obj.c(i)*dt);
            end

            %< Construct Rhs vector
            obj.setRhs(i, F=options.F, S=options.S, M=options.M, C=options.C, dt=dt, t=[t0+obj.c(i)*dt, t0]);

            %< Solve stage system
            U = obj.Solver{i}.solve();

            %< Evaluate operators for derivative computation
            Le = obj.evalOp(options.L, U, [t0+obj.c(i)*dt, t0]);
            Se = obj.evalOp(options.S, U, [t0+obj.c(i)*dt, t0]);
            Pe = obj.evalOp(options.P, U, [t0+obj.c(i)*dt, t0]);

            %< Compute derivative approximation for this stage
            obj.K{i} = 0;
            if ~isempty(Le)
                if ~isempty(options.P) || ~isempty(options.A) || ~isempty(options.B) || ~isempty(options.C)
                    obj.K{i} = obj.K{i} + Le * U(1:size(Le, 2));
                else
                    obj.K{i} = obj.K{i} + Le * U;
                end
            end
            if ~isempty(Pe)
                obj.K{i} = obj.K{i} + Pe * U(end-size(Pe,2)+1:end);
            end
            if ~isempty(Se)
                obj.K{i} = obj.K{i} + Se;
            end
        end

        function U = step(obj, options)
            % STEP Advance one time step using the DIRK method.
            %
            %   U = step(obj, L=L, F=F, S=S, M=M) performs a complete time step by
            %   computing all stages and combining them according to the
            %   Butcher tableau weights.
            %
            %   U = step(obj, L=L, F=F, S=S, M=M, P=P, A=A, B=B, C=C) performs a
            %   complete time step for augmented DIRK systems with constraints.

            arguments
                obj approx.odeint.DirkIntegrator
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

            % Detect augmented mode
            isAugmented = ~isempty(options.P) || ~isempty(options.A) || ~isempty(options.B) || ~isempty(options.C);

            t0 = obj.Timeline.Now;
            U0 = obj.U0{1};
            s = obj.NStages;

            %< Compute all stage approximations
            for i = 1:s
                [~] = obj.stage(i, L=options.L, F=options.F, S=options.S, M=options.M, P=options.P, A=options.A, B=options.B, C=options.C, dt=dt);
            end

            %< Update Lhs matrix for final solve if step size changed or custom dt provided
            if obj.Timeline.HasStepSizeChanged || ~isempty(options.dt)
                %< Evaluate Lhs mass term
                Me = obj.evalMLhs(options.M, t0+dt);

                if isAugmented
                    %< Construct augmented block matrix for final solve
                    Lhs = cell(2, 2);
                    Lhs{1, 1} = Me;
                    Lhs{1, 2} = [];
                    Lhs{2, 1} = options.A;
                    Lhs{2, 2} = options.B;
                    Lhs = core.linalg.block(Lhs);
                else
                    %< Regular mass matrix
                    Lhs = Me;
                end

                %< Set LHS in final solver and precompute factorization
                obj.Solver{s+1}.setLhs(Lhs);
                obj.Solver{s+1}.precompute();
            end

            %< Evaluate mass term for final solve
            Me = obj.evalMRhs(s+1, options.M, [t0+dt, t0]);

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

            %< Add weighted stage derivatives
            for i = 1:s
                Rhs1 = Rhs1 + dt * obj.b(i) * obj.K{i};
            end

            if isAugmented
                %< Construct augmented RHS vector for final solve
                Rhs = cell(2, 1);
                Rhs{1} = Rhs1;
                Ce = obj.evalOp(options.C, U0, [t0+dt, t0]);
                Rhs{2} = Ce;
                Rhs = core.linalg.block(Rhs);
            else
                %< Regular RHS
                Rhs = Rhs1;
            end

            obj.Solver{s+1}.setRhs(Rhs);

            %< Solve for final step approximation using stateful API
            U = obj.Solver{s+1}.solve();
        end
    end
end