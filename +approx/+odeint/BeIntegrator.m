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

            if ~isempty(options.fn)
                U = options.fn(obj, options);
                return;
            end

            if ~isempty(options.dt)
                dt = options.dt;
            else
                dt = obj.Timeline.StepSize;
            end

            U = obj.stage(options.L, options.F, options.S, options.M, ...
                options.P, options.A, options.B, options.C, dt);
        end
    end

    methods (Access = protected)
        function obj = setLhs(obj, L, M, P, A, B, dt, t, U)
            % SETLhs Set the left-hand side matrix for the Backward Euler system.

            if ~isempty(U)
                U0 = U;
            else
                U0 = obj.U0{1};
            end

            isAugmented = ~isempty(P) || ~isempty(A) || ~isempty(B);

            %< Evaluate mass operator for LHS
            Me = obj.evalMLhs(M, t);

            %< Evaluate linear operator if provided
            if ~isempty(L)
                Le = obj.evalOp(L, U0, t);
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
                if ~isempty(P)
                    Pe = obj.evalOp(P, [], t);
                    Lhs{1, 2} = -dt * Pe;
                else
                    Lhs{1, 2} = [];
                end

                %< Block (2,1): A
                Lhs{2, 1} = A;

                %< Block (2,2): B
                Lhs{2, 2} = B;

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

        function obj = setRhs(obj, F, S, M, C, dt, t, U)
            % SETRhs Set the right-hand side vector for the Backward Euler system.

            isAugmented = ~isempty(C);

            if ~isempty(U)
                U0 = U;
            else
                U0 = obj.U0{1};
            end

            %< Evaluate RHS mass operator
            Me = obj.evalMRhs([], M, t);

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
            if ~isempty(F)
                Fe = obj.evalOp(F, U0, t);
                Rhs1 = Rhs1 + dt * Fe;
            end

            %< Add source term contribution
            if ~isempty(S)
                Se = obj.evalOp(S, U0, t);
                Rhs1 = Rhs1 + dt * Se;
            end

            if isAugmented
                %< Construct augmented RHS vector
                Rhs = cell(2, 1);
                Rhs{1} = Rhs1;
                Rhs{2} = C;
                Rhs = core.linalg.block(Rhs);
            else
                %< Regular RHS
                Rhs = Rhs1;
            end

            obj.Solver.setRhs(Rhs);
        end

        function U = stage(obj, L, F, S, M, P, A, B, C, dt)
            % STAGE Compute the single stage of the Backward Euler method.

            U0 = obj.U0{1};
            t0 = obj.Timeline.Now;
            t1 = t0 + dt;

            %< Update Lhs matrix if step size changed
            if obj.Timeline.HasStepSizeChanged
                obj.setLhs(L, M, P, A, B, dt, t1, U0);
            end

            %< Construct Rhs vector
            obj.setRhs(F, S, M, C, dt, [t1, t0], U0);

            %< Solve the linear system
            U = obj.Solver.solve();
        end
    end
end
