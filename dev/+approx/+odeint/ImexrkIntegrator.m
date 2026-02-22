classdef ImexrkIntegrator < approx.odeint.OdeIntegrator
    % IMEXRKINTEGRATOR Implicit-Explicit Runge-Kutta integrator.
    %
    %   ImexrkIntegrator implements IMEX Runge-Kutta methods for solving
    %   ordinary differential equations of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + F(u) + S(t),
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator
    %   (treated implicitly), \f$F(u)\f$ is a nonlinear operator (treated
    %   explicitly), and \f$S(t)\f$ is a source term (treated explicitly).
    %   Any of these terms can be empty (indicating they vanish from the
    %   equation). This approach allows efficient handling of problems with
    %   both stiff linear parts and nonlinear parts that benefit from
    %   explicit treatment.
    %
    %   IMEX methods combine the stability of implicit methods for stiff
    %   linear terms with the efficiency of explicit methods for nonlinear
    %   terms, making them ideal for many PDE applications.
    %
    %   When L is empty, the method reduces to an explicit RK method.
    %   When F and S are empty, the method reduces to a DIRK method.
    %   Empty terms contribute zero to the corresponding evaluations.
    %
    % See also:
    %   approx.odeint.OdeIntegrator, approx.odeint.DirkIntegrator,
    %   approx.odeint.ExrkIntegrator

    properties
        AI % Butcher tableau for implicit part (lower triangular matrix)
        bI % Stage weights for implicit part (vector)
        AE % Butcher tableau for explicit part (strictly lower triangular matrix)
        bE % Stage weights for explicit part (vector)
        c % Stage time coefficients (vector)
        KI % Cell array of implicit derivative approximations
        KE % Cell array of explicit derivative approximations
    end

    methods
        function obj = ImexrkIntegrator(nStages, final)
            % IMEXRKINTEGRATOR Constructor for ImexrkIntegrator.
            %
            %   obj = ImexrkIntegrator(nStages, final) creates an IMEX
            %   Runge-Kutta integrator with the specified number of stages
            %   and final time.

            arguments
                nStages {mustBePositive, mustBeInteger}
                final {mustBePositive}
            end

            obj@approx.odeint.OdeIntegrator(1, nStages, final);
        end

        function obj = reset(obj, options)
            % RESET Reset the IMEX RK integrator to initial state.
            %
            %   obj = reset(obj) clears the derivative approximation
            %   storage for both implicit and explicit parts and calls the
            %   parent reset method.

            arguments
                obj approx.odeint.ImexrkIntegrator
                options.linOpts struct = struct() % Linear solver options
            end

            reset@approx.odeint.OdeIntegrator(obj);
            obj.KI = cell(obj.NStages, 1);
            obj.KE = cell(obj.NStages, 1);
            args = namedargs2cell(options.linOpts);
            obj.Solver = cell(1, obj.NStages+1);
            for i = 1 : obj.NStages+1
                obj.Solver{i} = core.linalg.LinearSolver(args{:});
            end
        end

        function obj = setLhs(obj, i, options)
            % SETLHS Set the left-hand side matrix for stage i.
            %
            %   obj = setLhs(obj, i, L=L, M=M) constructs the left-hand side
            %   matrix as M - dt*AI(i,i)*L where AI(i,i) is the diagonal
            %   element of the implicit Butcher tableau for stage i. If L
            %   is empty, LHS = M (identity operation).

            arguments
                obj approx.odeint.ImexrkIntegrator
                i {mustBePositive, mustBeInteger} % Stage index
                options.L = [] % Linear operator
                options.M = [] % Mass operator
                options.dt = [] % Time step override
                options.t = [] % Time values for time-dependent operators
            end

            % Use custom time step if provided
            if ~isempty(options.dt)
                dt = options.dt;
            else
                dt = obj.Timeline.StepSize;
            end

            %< Evaluate mass operator for LHS
            Me = obj.evalMLhs(options.M, options.t);

            %< Start with mass term
            Lhs = Me;

            %< Add linear term contribution with implicit coefficient
            if ~isempty(options.L)
                Le = obj.evalOp(options.L, [], options.t);
                Lhs = Lhs - dt * obj.AI(i, i) * Le;
            end

            obj.Solver{i}.setLhs(Lhs);
            obj.Solver{i}.precompute();
        end

        function obj = setRhs(obj, i, options)
            % SETRhs Set the right-hand side vector for stage i.
            %
            %   obj = setRhs(obj, i, M=M) constructs the right-hand side
            %   vector incorporating both implicit and explicit terms along
            %   with the mass operator. Empty terms contribute zero to the
            %   evaluation.

            arguments
                obj approx.odeint.ImexrkIntegrator
                i {mustBePositive, mustBeInteger} % Stage index
                options.M = [] % Mass operator
                options.dt = [] % Time step override
                options.t = [] % Time values
            end

            % Use custom time step if provided
            if ~isempty(options.dt)
                dt = options.dt;
            else
                dt = obj.Timeline.StepSize;
            end

            U0 = obj.U0{1};

            %< Evaluate RHS mass operator
            Me = obj.evalMRhs(i, options.M, options.t);

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

            %< Accumulate terms from previous implicit stages
            for j = 1:(i - 1)
                Rhs = Rhs + dt * obj.AI(i, j) * obj.KI{j};
            end

            %< Accumulate terms from explicit stages
            for j = 1:i
                Rhs = Rhs + dt * obj.AE(i + 1, j) * obj.KE{j};
            end

            obj.Solver{i}.setRhs(Rhs);
        end

        function U = stage(obj, i, options)
            % STAGE Compute stage i of the IMEX RK method.
            %
            %   U = stage(obj, i, L=L, F=F, S=S, M=M) solves the implicit system
            %   for the derivative approximation at stage i and evaluates
            %   the explicit terms. Empty terms are handled gracefully.

            arguments
                obj approx.odeint.ImexrkIntegrator
                i {mustBePositive, mustBeInteger} % Stage number
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
            U0 = obj.U0{1};

            % -------------------------------------------------------------
            % IMPLICIT PART
            % -------------------------------------------------------------
            %< Evaluate linear term at stage time
            Le = obj.evalOp(options.L, U0, [t0+obj.c(i)*dt, t0]);

            %< Update LHS matrix if step size changed or custom dt provided
            if obj.Timeline.HasStepSizeChanged || ~isempty(options.dt)
                %< Construct LHS matrix
                obj.setLhs(i, L=options.L, M=options.M, dt=dt, t=t0+obj.c(i)*dt);
            end

            %< Construct Rhs vector
            obj.setRhs(i, M=options.M, dt=dt, t=[t0+obj.c(i)*dt, t0]);

            %< Solve for implicit stage approximation
            U = obj.Solver{i}.solve();

            %< Compute implicit derivative approximation
            obj.KI{i} = 0;
            if ~isempty(Le)
                obj.KI{i} = obj.KI{i} + Le * U;
            end

            % -------------------------------------------------------------
            % EXPLICIT PART
            % -------------------------------------------------------------

            %< Evaluate source term at next stage time
            Se = obj.evalOp(options.S, U, [t0+obj.c(i+1)*dt, t0]);

            %< Evaluate nonlinear term at current stage approximation
            Fe = obj.evalOp(options.F, U, [t0+obj.c(i)*dt, t0]);

            %< Compute explicit derivative approximation
            obj.KE{i + 1} = 0;

            %< Add nonlinear term contribution
            if ~isempty(Fe)
                obj.KE{i + 1} = obj.KE{i + 1} + Fe;
            end

            %< Add source term contribution
            if ~isempty(Se)
                obj.KE{i + 1} = obj.KE{i + 1} + Se;
            end
        end

        function U = step(obj, options)
            % STEP Advance one time step using the IMEX RK method.
            %
            %   U = step(obj, L=L, F=F, S=S, M=M) performs a complete time step by
            %   solving implicit stages for the linear term and explicit
            %   evaluations for the nonlinear term, then combining
            %   according to the Butcher tableau weights. Empty terms are
            %   handled gracefully by contributing zero.

            arguments
                obj approx.odeint.ImexrkIntegrator
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
            U0 = obj.U0{1};
            r = obj.NStages;
            s = r - 1;

            %< Initialize first stage for explicit part
            obj.KE{1} = 0;

            %< Evaluate source term at initial time
            Se = obj.evalOp(options.S, U0, [t0+obj.c(1)*dt, t0]);

            %< Evaluate nonlinear term at initial solution
            Fe = obj.evalOp(options.F, U0, [t0, t0]);

            %< Initialize explicit stage derivative approximation
            obj.KE{1} = 0;

            %< Add nonlinear term contribution
            if ~isempty(Fe)
                obj.KE{1} = obj.KE{1} + Fe;
            end

            %< Add source term contribution
            if ~isempty(Se)
                obj.KE{1} = obj.KE{1} + Se;
            end

            %< Compute all subsequent stages
            for i = 1:s
                [~] = obj.stage(i, L=options.L, F=options.F, S=options.S, M=options.M, dt=dt);
            end

            %< Update LHS matrix for final solve if step size changed
            if obj.Timeline.HasStepSizeChanged
                %< Construct LHS matrix
                obj.setLhs(s+1, L=[], M=options.M, dt=dt, t=t0+dt);
            end

            %< Evaluate mass term for final solve
            Me = obj.evalMRhs(s+1, options.M, [t0 + dt, t0]);

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

            %< Add weighted explicit derivative terms
            for i = 1:r
                Rhs = Rhs + dt * obj.bE(i) * obj.KE{i};
            end

            %< Add weighted implicit derivative terms
            for i = 1:s
                Rhs = Rhs + dt * obj.bI(i) * obj.KI{i};
            end

            obj.Solver{s+1}.setRhs(Rhs);

            %< Solve for final step approximation
            U = obj.Solver{s+1}.solve();
        end
    end
end