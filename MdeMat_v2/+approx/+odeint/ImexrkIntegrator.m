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
    % Notes:
    %   When L is empty, the method reduces to an explicit RK method.
    %   When F and S are empty, the method reduces to a DIRK method.
    %   Empty terms contribute zero to the corresponding evaluations.
    %
    % Examples:
    %   % Basic IMEX problem
    %   integrator = approx.odeint.Ars222Integrator(1.0);
    %   U_new = integrator.step(L, F, S, M);
    %   integrator.update(U_new);
    %   integrator.advance();
    %
    %   % Pure implicit diffusion
    %   U_new = integrator.step(L, [], [], []);
    %
    %   % Mixed implicit-explicit (L implicit, F explicit)
    %   U_new = integrator.step(L, F, [], []);
    %
    %   % Implicit linear with explicit forcing
    %   U_new = integrator.step(L, [], S, []);
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
            %
            % Inputs:
            %   nStages - Number of stages for the IMEX method (positive integer)
            %   final - Final time for integration (positive scalar)
            %
            % Outputs:
            %   obj - Constructed ImexrkIntegrator object

            obj@approx.odeint.OdeIntegrator(1, nStages, final);
        end

        function obj = reset(obj)
            % RESET Reset the IMEX RK integrator to initial state.
            %
            %   obj = reset(obj) clears the derivative approximation
            %   storage for both implicit and explicit parts and calls the
            %   parent reset method.
            %
            % Inputs:
            %   obj - The ImexrkIntegrator object
            %
            % Outputs:
            %   obj - The ImexrkIntegrator object

            reset@approx.odeint.OdeIntegrator(obj);
            obj.KI = cell(obj.nStages, 1);
            obj.KE = cell(obj.nStages, 1);
        end

        function obj = setLhs(obj, i, L, M)
            % SETLHS Set the left-hand side matrix for stage i.
            %
            %   obj = setLhs(obj, i, L, M) constructs the left-hand side
            %   matrix as M - dt*AI(i,i)*L where AI(i,i) is the diagonal
            %   element of the implicit Butcher tableau for stage i. If L
            %   is empty, LHS = M (identity operation).
            %
            % Inputs:
            %   obj - The ImexrkIntegrator object
            %   i - Stage index (positive integer)
            %   L - Linear operator (matrix or empty)
            %   M - Mass operator (matrix or empty)
            %
            % Outputs:
            %   obj - The ImexrkIntegrator object

            dt = obj.timeline.dt;

            obj.lhs{i} = M;

            %< Add linear term contribution with implicit coefficient
            if ~isempty(L)
                obj.lhs{i} = obj.lhs{i} - dt * obj.AI(i, i) * L;
            end
        end

        function obj = setRhs(obj, i, M)
            % SETRHS Set the right-hand side vector for stage i.
            %
            %   obj = setRhs(obj, i, M) constructs the right-hand side
            %   vector incorporating both implicit and explicit terms along
            %   with the mass operator. Empty terms contribute zero to the
            %   evaluation.
            %
            % Inputs:
            %   obj - The ImexrkIntegrator object
            %   i - Stage index (positive integer)
            %   M - Mass operator (matrix, function, or empty)
            %
            % Outputs:
            %   obj - The ImexrkIntegrator object

            dt = obj.timeline.dt;
            U0 = obj.U0{1};

            %< Add mass term contribution from previous solution
            if isempty(M)
                obj.rhs{i} = U0;
            elseif iscell(M)
                obj.rhs{i} = M{2} * U0;
                if ~isempty(M{1})
                    obj.rhs{i} = obj.rhs{i} + M{1};
                end
            else
                obj.rhs{i} = M * U0;
            end

            %< Accumulate terms from previous implicit stages
            for j = 1:(i - 1)
                obj.rhs{i} = obj.rhs{i} + dt * obj.AI(i, j) * obj.KI{j};
            end

            %< Accumulate terms from explicit stages
            for j = 1:i
                obj.rhs{i} = obj.rhs{i} + dt * obj.AE(i + 1, j) * obj.KE{j};
            end
        end

        function U = stage(obj, i, L, F, S, M)
            % STAGE Compute stage i of the IMEX RK method.
            %
            %   U = stage(obj, i, L, F, S, M) solves the implicit system for
            %   the derivative approximation at stage i and evaluates the
            %   explicit terms. Empty terms are handled gracefully.
            %
            % Inputs:
            %   obj - The ImexrkIntegrator object
            %   i - Stage index (positive integer)
            %   L - Linear operator (matrix, function, or empty)
            %   F - Nonlinear operator (function or empty)
            %   S - Source term (function, vector, or empty)
            %   M - Mass operator (matrix, function, or empty)
            %
            % Outputs:
            %   U - Solution approximation at stage i

            t0 = obj.timeline.now;
            dt = obj.timeline.dt;
            U0 = obj.U0{1};

            % -------------------------------------------------------------
            % IMPLICIT PART
            % -------------------------------------------------------------
            %< Evaluate linear term at stage time
            Le = obj.evaluateOperator(L, U0, t0+obj.c(i)*dt);

            %< Update LHS matrix if step size changed
            if obj.timeline.hasStepSizeChanged
                %< Evaluate LHS mass term
                Me = obj.evaluateLhsMass(M, t0+obj.c(i)*dt);

                %< Construct LHS matrix
                obj.setLhs(i, Le, Me);
            end

            %< Evaluate RHS mass term at previous time level
            Me = obj.evaluateRhsMass(i, M, t0+obj.c(i)*dt);

            %< Construct RHS vector
            obj.setRhs(i, Me);

            %< Solve for implicit stage approximation
            U = obj.solver.solve(obj.lhs{i}, obj.rhs{i});

            %< Compute implicit derivative approximation
            obj.KI{i} = 0;
            if ~isempty(Le)
                obj.KI{i} = obj.KI{i} + Le * U;
            end

            % -------------------------------------------------------------
            % EXPLICIT PART
            % -------------------------------------------------------------

            %< Evaluate source term at next stage time
            Se = obj.evaluateOperator(S, U, t0+obj.c(i+1)*dt);

            %< Evaluate nonlinear term at current stage approximation
            Fe = obj.evaluateOperator(F, U, t0+obj.c(i)*dt);

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

        function U = step(obj, L, F, S, M)
            % STEP Advance one time step using the IMEX RK method.
            %
            %   U = step(obj, L, F, S, M) performs a complete time step by
            %   solving implicit stages for the linear term and explicit
            %   evaluations for the nonlinear term, then combining according
            %   to the Butcher tableau weights. Empty terms are handled
            %   gracefully by contributing zero.
            %
            % Inputs:
            %   obj - The ImexrkIntegrator object
            %   L - Linear operator (matrix, function, or empty)
            %   F - Nonlinear operator (function or empty)
            %   S - Source term (function, vector, or empty)
            %   M - Mass operator (matrix, function, or empty)
            %
            % Outputs:
            %   U - Solution at the next time step

            t0 = obj.timeline.now;
            dt = obj.timeline.dt;
            U0 = obj.U0{1};
            r = obj.nStages;
            s = r - 1;

            %< Initialize first stage for explicit part
            obj.KE{1} = 0;

            %< Evaluate source term at initial time
            Se = obj.evaluateOperator(S, U0, t0+obj.c(1)*dt);

            %< Evaluate nonlinear term at initial solution
            Fe = obj.evaluateOperator(F, U0, t0);

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
                [~] = obj.stage(i, L, F, S, M);
            end

            %< Update LHS matrix for final solve if step size changed
            if obj.timeline.hasStepSizeChanged
                %< Evaluate LHS mass term
                Me = obj.evaluateLhsMass(M, t0+dt);

                %< Construct LHS matrix
                obj.setLhs(s+1, [], Me);
            end

            %< Evaluate mass term for final solve
            Me = obj.evaluateRhsMass(s+1, M, [t0 + dt, t0]);

            %< Add mass term contribution from previous solution
            if isempty(Me)
                obj.rhs{s+1} = U0;
            elseif iscell(Me)
                obj.rhs{s+1} = Me{2} * U0;
                if ~isempty(Me{1})
                    obj.rhs{s+1} = obj.rhs{s+1} + Me{1};
                end
            else
                obj.rhs{s+1} = Me * U0;
            end

            %< Add weighted explicit derivative terms
            for i = 1:r
                obj.rhs{s+1} = obj.rhs{s+1} + dt * obj.bE(i) * obj.KE{i};
            end

            %< Add weighted implicit derivative terms
            for i = 1:s
                obj.rhs{s+1} = obj.rhs{s+1} + dt * obj.bI(i) * obj.KI{i};
            end

            %< Solve for final step approximation
            U = obj.solver.solve(obj.lhs{s+1}, obj.rhs{s+1});
        end
    end
end