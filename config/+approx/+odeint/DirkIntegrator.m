classdef DirkIntegrator < approx.odeint.OdeIntegrator
    % DIRKINTEGRATOR Diagonally Implicit Runge-Kutta integrator.
    %
    %   DirkIntegrator implements Diagonally Implicit Runge-Kutta (DIRK)
    %   methods for solving ordinary differential equations of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + S,
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator, and
    %   \f$S\f$ is a source term. Any of these terms can be empty. For mass
    %   operator, it means it reduces to an identity. For linear and source
    %   term, it means they vanish from the equation.
    %
    %   DIRK methods are characterized by a lower triangular Butcher
    %   tableau with identical diagonal elements, allowing each stage to be
    %   solved independently. This makes them more efficient than fully
    %   implicit Runge-Kutta methods while maintaining good stability
    %   properties.
    %
    % Examples:
    %   % Basic usage
    %   integrator = approx.odeint.Sdirk2Integrator(1.0);
    %   U_new = integrator.step(L, S, M);
    %   integrator.update(U_new);
    %   integrator.advance();
    %
    %   % Pure diffusion (only linear term)
    %   U_new = integrator.step(L, [], []);
    %
    %   % Pure forcing (only source term)
    %   U_new = integrator.step([], S, []);
    %
    %   % Identity evolution (no dynamics)
    %   U_new = integrator.step([], [], []);
    %
    %   % With mass matrix
    %   U_new = integrator.step(L, S, M);
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
            %
            % Inputs:
            %   nStages - Number of stages for the DIRK method (positive integer)
            %   final - Final time for integration (positive scalar)
            %
            % Outputs:
            %   obj - Constructed DirkIntegrator object

            obj@approx.odeint.OdeIntegrator(1, nStages, final);
        end

        function obj = reset(obj)
            % RESET Reset the DIRK integrator to initial state.
            %
            %   obj = reset(obj) clears the derivative approximation
            %   storage and calls the parent reset method.
            %
            % Inputs:
            %   obj - The DirkIntegrator object
            %
            % Outputs:
            %   obj - The DirkIntegrator object

            reset@approx.odeint.OdeIntegrator(obj);
            obj.K = cell(obj.nStages, 1);
        end

        function obj = setLhs(obj, i, L, M)
            % SETLHS Set the left-hand side matrix for stage i.
            %
            %   obj = setLhs(obj, i, L, M) constructs the left-hand side
            %   matrix as M - dt*A(i,i)*L for stage i of the DIRK system.
            %
            % Inputs:
            %   obj - The DirkIntegrator object
            %   i - Stage index (positive integer)
            %   L - Linear operator (matrix or empty)
            %   M - Mass operator (matrix or empty)
            %
            % Outputs:
            %   obj - The DirkIntegrator object

            dt = obj.timeline.dt;

            %< Add mass term
            obj.lhs{i} = M;

            %< Add linear term contribution with diagonal coefficient
            if ~isempty(L)
                obj.lhs{i} = obj.lhs{i} - dt * obj.A(i, i) * L;
            end
        end

        function obj = setRhs(obj, i, S, M)
            % SETRHS Set the right-hand side vector for stage i.
            %
            %   obj = setRhs(obj, i, S, M) constructs the right-hand side
            %   vector using mass operator, accumulated previous stages,
            %   and source term.
            %
            % Inputs:
            %   obj - The DirkIntegrator object
            %   i - Stage index (positive integer)
            %   S - Source term (vector, function, or empty)
            %   M - Mass operator (matrix, function, or empty)
            %
            % Outputs:
            %   obj - The DirkIntegrator object

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

            %< Add accumulated term from previous stages
            for j = 1:(i - 1)
                obj.rhs{i} = obj.rhs{i} + dt * obj.A(i, j) * obj.K{j};
            end

            %< Add source term contribution with diagonal coefficient
            if ~isempty(S)
                obj.rhs{i} = obj.rhs{i} + dt * obj.A(i, i) * S;
            end
        end

        function U = stage(obj, i, L, S, M)
            % STAGE Compute stage i of the DIRK method.
            %
            %   U = stage(obj, i, L, S, M) computes the solution at stage i
            %   by solving the implicit system and computing the derivative
            %   approximation for use in subsequent stages.
            %
            % Inputs:
            %   obj - The DirkIntegrator object
            %   i - Stage index (positive integer)
            %   L - Linear operator (matrix, function, or empty)
            %   S - Source term (vector, function, or empty)
            %   M - Mass operator (matrix, function, or empty)
            %
            % Outputs:
            %   U - Solution approximation at stage i

            t0 = obj.timeline.now;
            dt = obj.timeline.dt;
            U0 = obj.U0{1};

            %< Evaluate linear term at stage time
            Le = obj.evaluateOperator(L, U0, [t0+obj.c(i)*dt, t0]);

            %< Evaluate source term at stage time
            Se = obj.evaluateOperator(S, U0, [t0+obj.c(i)*dt, t0]);

            %< Update LHS matrix if step size changed
            if obj.timeline.hasStepSizeChanged
                %< Evaluate LHS mass term
                Me = obj.evaluateLhsMass(M, t0+obj.c(i)*dt);

                %< Construct LHS matrix
                obj.setLhs(i, Le, Me);
            end

            %< Evaluate RHS mass term at previous time level
            Me = obj.evaluateRhsMass(i, M, [t0 + obj.c(i)*dt, t0]);

            %< Construct RHS vector
            obj.setRhs(i, Se, Me);

            %< Solve linear system for stage approximation
            U = obj.solver.solve(obj.lhs{i}, obj.rhs{i});

            %< Compute derivative approximation for this stage
            obj.K{i} = 0;
            if ~isempty(Le)
                obj.K{i} = obj.K{i} + Le * U;
            end
            if ~isempty(Se)
                obj.K{i} = obj.K{i} + Se;
            end
        end

        function U = step(obj, L, S, M)
            % STEP Advance one time step using the DIRK method.
            %
            %   U = step(obj, L, S, M) performs a complete time step by
            %   computing all stages and combining them according to the
            %   Butcher tableau weights.
            %
            % Inputs:
            %   obj - The DirkIntegrator object
            %   L - Linear operator (matrix, function, or empty)
            %   S - Source term (vector, function, or empty)
            %   M - Mass operator (matrix, function, or empty)
            %
            % Outputs:
            %   U - Solution at the next time step

            t0 = obj.timeline.now;
            dt = obj.timeline.dt;
            U0 = obj.U0{1};
            s = obj.nStages;

            %< Compute all stage approximations
            for i = 1:s
                [~] = obj.stage(i, L, S, M);
            end

            %< Update LHS matrix for final solve if step size changed
            if obj.timeline.hasStepSizeChanged
                %< Evaluate LHS mass term
                Me = obj.evaluateLhsMass(M, t0+dt);
                obj.lhs{s+1} = Me;
            end

            %< Evaluate mass term for final solve
            Me = obj.evaluateRhsMass(s+1, M, [t0+dt, t0]);

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

            %< Add weighted stage derivatives
            for i = 1:s
                obj.rhs{s+1} = obj.rhs{s+1} + dt * obj.b(i) * obj.K{i};
            end

            %< Solve for final step approximation
            U = obj.solver.solve(obj.lhs{s+1}, obj.rhs{s+1});
        end
    end
end