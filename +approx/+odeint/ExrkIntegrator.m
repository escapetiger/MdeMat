classdef ExrkIntegrator < approx.odeint.OdeIntegrator
    % EXRKINTEGRATOR Explicit Runge-Kutta integrator.
    %
    %   ExrkIntegrator implements explicit Runge-Kutta methods for solving
    %   ordinary differential equations of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + F(u) + S,
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator,
    %   \f$F(u)\f$ is a nonlinear operator, and \f$S\f$ is a source term.
    %   Any of these terms can be empty (indicating they vanish from the
    %   equation). When M is empty, it defaults to identity. All non-empty
    %   terms are treated explicitly, making this suitable for non-stiff
    %   problems where stability constraints are manageable.
    %
    %   Explicit RK methods have excellent stability properties for
    %   non-stiff problems but may require small time steps for stability
    %   when applied to stiff systems.
    %
    %   Augmented systems with constraints:
    %
    %   \f[
    %     \begin{bmatrix}
    %       M \\
    %       A
    %     \end{bmatrix}
    %     \frac{\mathrm{d}}{\mathrm{d}t}
    %     \begin{bmatrix}
    %       u \\
    %       \lambda
    %     \end{bmatrix}
    %     =
    %     \begin{bmatrix}
    %       Lu + P\lambda + F(u) + S \\
    %       C
    %     \end{bmatrix}
    %   \f]
    %
    %   where P, A, B, C represent constraint coupling terms.
    %
    % See also:
    %   approx.odeint.OdeIntegrator, approx.odeint.DirkIntegrator,
    %   approx.odeint.FeIntegrator

    properties
        A % Butcher tableau coefficient matrix (strictly lower triangular)
        b % Stage weights vector for final update
        c % Stage time coefficients vector
        K % Cell array of derivative approximations for each stage
        AugmentSolver % Linear solver for augmented constraint solve
    end

    methods
        function obj = ExrkIntegrator(nStages, final)
            % EXRKINTEGRATOR Constructor for ExrkIntegrator.
            %
            %   obj = ExrkIntegrator(nStages, final) creates an explicit
            %   Runge-Kutta integrator with the specified number of stages
            %   and final time.

            arguments
                nStages {mustBePositive, mustBeInteger}
                final {mustBeNumeric, mustBeNonnegative}
            end

            obj@approx.odeint.OdeIntegrator(1, nStages, final);
        end

        function obj = reset(obj, options)
            % RESET Reset the EXRK integrator to initial state.
            %
            %   obj = reset(obj) clears the derivative approximation
            %   storage and calls the parent reset method.

            arguments
                obj approx.odeint.ExrkIntegrator
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
            % STEP Advance one time step using the explicit RK method.
            %
            %   U = step(obj, L=L, F=F, S=S, M=M) performs a complete time step by
            %   evaluating all stages explicitly and combining them according
            %   to the Butcher tableau weights. Empty terms contribute zero
            %   to all stages.
            %
            %   U = step(obj, L=L, F=F, S=S, M=M, P=P, A=A, B=B, C=C) performs
            %   a complete time step for augmented systems with constraints.

            arguments
                obj approx.odeint.ExrkIntegrator
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

            t0 = obj.Timeline.Now;
            U0 = obj.U0{1};
            s = obj.NStages;

            %< Compute all stage approximations
            for i = 1:s
                obj.stage(i, options.L, options.F, options.S, options.M, ...
                    options.P, options.A, options.B, options.C, dt, []);
            end

            isAugmented = ~isempty(options.P) || ~isempty(options.A) || ~isempty(options.B) || ~isempty(options.C);

            %< Update Lhs matrix for final solve if step size changed or custom dt provided
            if obj.Timeline.HasStepSizeChanged || ~isempty(options.dt)
                %< Evaluate Lhs mass term
                Me = obj.evalMLhs(options.M, t0+dt);

                %< Construct Lhs matrix
                obj.setLhs(s+1, Me, dt, t0+dt);
            end

            %< Evaluate mass term for final solve
            Me = obj.evalMRhs(s+1, options.M, [t0+dt, t0]);

            %< Add mass term contribution from previous solution
            if isAugmented
                if isempty(Me)
                    Rhs = U0(1:size(Me, 2));
                elseif iscell(Me)
                    Rhs = Me{2} * U0(1:size(Me{2}, 2));
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

            %< Solve for final step approximation using stateful API
            if isAugmented
                U = obj.Solver{s+1}.solve();
                obj.AugmentSolver.setLhs(options.B);
                obj.AugmentSolver.setRhs(-options.A * U + options.C);
                Ua = obj.AugmentSolver.solve();
                U = [U; Ua];
            else
                U = obj.Solver{s + 1}.solve();
            end
        end
    end

    methods (Access = protected)
        function obj = setLhs(obj, i, M, dt, t)
            % SETLhs Set the mass matrix for explicit solve at stage i.

            %< Evaluate mass operator for LHS
            Me = obj.evalMLhs(M, t);

            obj.Solver{i}.setLhs(Me);
            obj.Solver{i}.precompute();
        end

        function obj = setRhs(obj, i, M, C, dt, t)
            % SETRhs Set the right-hand side vector for stage i.

            U0 = obj.U0{1};

            %< Evaluate RHS mass operator
            Me = obj.evalMRhs(i, M, t);

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

            %< Add accumulated term from previous stages
            for j = 1:(i - 1)
                Rhs = Rhs + dt * obj.A(i, j) * obj.K{j};
            end

            obj.Solver{i}.setRhs(Rhs);
        end

        function U = stage(obj, i, L, F, S, M, P, A, B, C, dt, U)
            % STAGE Compute stage i of the explicit RK method.

            if ~isempty(U)
                U0 = U;
            else
                U0 = obj.U0{1};
            end

            t0 = obj.Timeline.Now;

            %< Evaluate operators at U0 (previous time step solution)
            Le = obj.evalOp(L, U0, [t0, t0]);
            Se = obj.evalOp(S, U0, [t0+obj.c(i)*dt, t0]);
            Pe = obj.evalOp(P, U0, [t0, t0]);

            %< Update Lhs matrix if step size changed
            if obj.Timeline.HasStepSizeChanged
                obj.setLhs(i, M, dt, t0+obj.c(i)*dt);
            end

            %< Construct Rhs vector
            obj.setRhs(i, M, C, dt, [t0+obj.c(i)*dt, t0]);

            %< Solve for stage approximation using stateful API
            U = obj.Solver{i}.solve();

            %< Evaluate nonlinear term at stage approximation
            Fe = obj.evalOp(F, U, [t0+obj.c(i)*dt, t0]);

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
            if ~isempty(Fe)
                obj.K{i} = obj.K{i} + Fe;
            end
            if ~isempty(Se)
                obj.K{i} = obj.K{i} + Se;
            end
        end
    end
end
