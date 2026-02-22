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
        AugmentSolver
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
            obj.AugmentSolver = core.linalg.LinearSolver(args{:});
        end

        function U = step(obj, options)
            % STEP Advance one time step using the DIRK method.
            %
            %   U = step(obj, L=L, F=F, S=S, M=M) performs a complete time
            %   step by computing all stages and combining them according
            %   to the Butcher tableau weights.
            %
            %   U = step(obj, L=L, F=F, S=S, M=M, P=P, A=A, B=B, C=C)
            %   performs a complete time step for augmented DIRK systems
            %   with constraints.

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

            if ~isempty(options.fn)
                U = options.fn(obj, options);
                return;
            end

            if ~isempty(options.dt)
                dt = options.dt;
            else
                dt = obj.Timeline.StepSize;
            end

            isAugmented = ~isempty(options.P) || ~isempty(options.A) || ~isempty(options.B) || ~isempty(options.C);

            t0 = obj.Timeline.Now;
            U0 = obj.U0{1};
            s = obj.NStages;

            %< Compute all stage approximations
            U1 = U0;
            for i = 1:s
                U1 = obj.stage(i, options.L, options.F, options.S, options.M, ...
                    options.P, options.A, options.B, options.C, dt, U1);
            end

            %< Update Lhs matrix for final solve if step size changed or custom dt provided
            if obj.Timeline.HasStepSizeChanged || ~isempty(options.dt)
                obj.setLhs(s+1, [], options.M, [], [], [], dt, t0+dt, []);
            end

            %< Evaluate mass term for final solve
            Me = obj.evalMRhs(s+1, options.M, [t0+dt, t0]);

            %< Add mass term contribution from previous solution
            if isAugmented
                if isempty(Me)
                    Rhs = U0(1:size(Me, 2));
                elseif iscell(Me)
                    Rhs = Me{2} * U0(1:size(Me, 2));
                    if ~isempty(Me{1})
                        Rhs = Rhs + Me{1};
                    end
                else
                    Rhs = Me * U0(1:size(Me, 2));
                end
            else
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
            end

            %< Add weighted stage derivatives
            for i = 1:s
                Rhs = Rhs + dt * obj.b(i) * obj.K{i};
            end

            obj.Solver{s+1}.setRhs(Rhs);

            %< Solve for final step approximation
            if isAugmented
                U = obj.Solver{s+1}.solve();
                obj.AugmentSolver.setLhs(options.B);
                obj.AugmentSolver.setRhs(-options.A * U + options.C);
                Ua = obj.AugmentSolver.solve();
                U = [U; Ua];
            else
                U = obj.Solver{s+1}.solve();
            end
        end
    end

    methods (Access = protected)
        function obj = setLhs(obj, i, L, M, P, A, B, dt, t, U)
            % SETLhs Set the left-hand side matrix for stage i.

            U1 = U;
            isAugmented = ~isempty(P) || ~isempty(A) || ~isempty(B);

            %< Evaluate mass operator for LHS
            Me = obj.evalMLhs(M, t);

            %< Evaluate linear operator if provided
            if ~isempty(L)
                Le = obj.evalOp(L, U1, t);
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
                if ~isempty(P)
                    Pe = obj.evalOp(P, U1, t);
                    Lhs{1, 2} = -dt * obj.A(i, i) * Pe;
                else
                    Lhs{1, 2} = [];
                end

                %< Block (2,1): A
                Lhs{2, 1} = A;

                %< Block (2,2): B
                Lhs{2, 2} = B;

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

        function obj = setRhs(obj, i, F, S, M, C, dt, t, U)
            % SETRhs Set the right-hand side vector for stage i.

            isAugmented = ~isempty(C);

            U0 = obj.U0{1};
            U1 = U;

            %< Evaluate RHS mass operator
            Me = obj.evalMRhs(i, M, t);

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
            if ~isempty(F)
                Fe = obj.evalOp(F, U1, t);
                Rhs1 = Rhs1 + dt * obj.A(i, i) * Fe;
            end

            %< Add source term contribution with diagonal coefficient
            if ~isempty(S)
                Se = obj.evalOp(S, U1, t);
                Rhs1 = Rhs1 + dt * obj.A(i, i) * Se;
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

            obj.Solver{i}.setRhs(Rhs);
        end

        function U = stage(obj, i, L, F, S, M, P, A, B, C, dt, U)
            % STAGE Compute stage i of the DIRK method.

            if ~isempty(U)
                U1 = U;
            else
                U1 = obj.U0{1};
            end

            t0 = obj.Timeline.Now;

            %< Update Lhs matrix if step size changed or custom dt provided
            if obj.Timeline.HasStepSizeChanged
                obj.setLhs(i, L, M, P, A, B, dt, t0+obj.c(i)*dt, U1);
            end

            %< Construct Rhs vector
            obj.setRhs(i, F, S, M, C, dt, [t0+obj.c(i)*dt, t0], U1);

            %< Solve stage system
            U = obj.Solver{i}.solve();

            %< Evaluate operators for derivative computation
            Le = obj.evalOp(L, U1, [t0+obj.c(i)*dt, t0]);
            Se = obj.evalOp(S, U1, [t0+obj.c(i)*dt, t0]);
            Pe = obj.evalOp(P, U1, [t0+obj.c(i)*dt, t0]);

            %< Compute derivative approximation for this stage
            obj.K{i} = 0;
            if ~isempty(Le)
                if ~isempty(P) || ~isempty(A) || ~isempty(B) || ~isempty(C)
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
    end
end
