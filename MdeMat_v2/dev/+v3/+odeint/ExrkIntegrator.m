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
    % Examples:
    %   % Basic usage with all terms
    %   integrator = approx.odeint.Exrk4Integrator(1.0);
    %   U_new = integrator.step(L, F, S, M);
    %   integrator.update(U_new);
    %   integrator.advance();
    %
    %   % Pure diffusion (only linear term)
    %   U_new = integrator.step(L, [], [], []);
    %
    %   % Pure nonlinear dynamics
    %   U_new = integrator.step([], F, [], []);
    %
    %   % Linear system with forcing
    %   U_new = integrator.step(L, [], S, []);
    %
    %   % With mass matrix
    %   U_new = integrator.step(L, F, S, M);
    %
    % See also:
    %   approx.odeint.OdeIntegrator, approx.odeint.DirkIntegrator,
    %   approx.odeint.FeIntegrator

    properties
        A % Butcher tableau coefficient matrix (strictly lower triangular)
        b % Stage weights vector for final update
        c % Stage time coefficients vector
        K % Cell array of derivative approximations for each stage
    end

    methods
        function obj = ExrkIntegrator(nStages, final)
            % EXRKINTEGRATOR Constructor for ExrkIntegrator.
            %
            %   obj = ExrkIntegrator(nStages, final) creates an explicit
            %   Runge-Kutta integrator with the specified number of stages
            %   and final time.
            %
            % Inputs:
            %   nStages - Number of stages for the explicit RK method (positive integer)
            %   final - Final time for integration (positive scalar)
            %
            % Outputs:
            %   obj - Constructed ExrkIntegrator object

            obj@approx.odeint.OdeIntegrator(1, nStages, final);
        end

        function obj = reset(obj)
            % RESET Reset the explicit RK integrator to initial state.
            %
            %   obj = reset(obj) clears the derivative approximation
            %   storage and calls the parent reset method.
            %
            % Inputs:
            %   obj - The ExrkIntegrator object
            %
            % Outputs:
            %   obj - The ExrkIntegrator object

            reset@approx.odeint.OdeIntegrator(obj);
            obj.K = cell(obj.nStages, 1);
        end

        function obj = setLhs(obj, i, M)
            % SETLHS Set the mass matrix for explicit solve at stage i.
            %
            %   obj = setLhs(obj, i, M) sets the left-hand side matrix to
            %   the mass operator for stage i. For explicit methods, this
            %   is used to solve M*K = RHS for the derivative
            %   approximation.
            %
            % Inputs:
            %   obj - The ExrkIntegrator object
            %   i - Stage index (positive integer)
            %   M - Mass operator (matrix or empty)
            %
            % Outputs:
            %   obj - The ExrkIntegrator object

            obj.lhs{i} = M;
        end

        function obj = setRhs(obj, i, M)
            % SETRHS Set the right-hand side vector for stage i.
            %
            %   obj = setRhs(obj, i, M) constructs the right-hand side
            %   vector using mass operator and accumulated previous stages.
            %
            % Inputs:
            %   obj - The ExrkIntegrator object
            %   i - Stage index (positive integer)
            %   M - Mass operator (matrix, function, or empty)
            %
            % Outputs:
            %   obj - The ExrkIntegrator object

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
        end

        function U = stage(obj, i, L, F, S, M)
            % STAGE Compute stage i of the explicit RK method.
            %
            %   U = stage(obj, i, L, F, S, M) computes the solution at stage i
            %   by evaluating all terms explicitly and computing the derivative
            %   approximation for use in subsequent stages.
            %
            % Inputs:
            %   obj - The ExrkIntegrator object
            %   i - Stage index (positive integer)
            %   L - Linear operator (matrix, function, or empty)
            %   F - Nonlinear operator (function or empty)
            %   S - Source term (vector, function, or empty)
            %   M - Mass operator (matrix, function, or empty)
            %
            % Outputs:
            %   U - Solution approximation at stage i

            t0 = obj.timeline.now;
            dt = obj.timeline.dt;
            U0 = obj.U0{1};

            %< Evaluate linear term at initial time
            if isempty(L)
                Le = [];
            elseif isa(L, 'function_handle')
                Le = L(t0);
            else
                Le = L;
            end

            %< Evaluate source term at stage time
            if isempty(S)
                Se = [];
            elseif isa(S, 'function_handle')
                Se = S(t0 + obj.c(i) * dt);
            else
                Se = S;
            end

            %< Update LHS matrix if step size changed
            if obj.timeline.hasStepSizeChanged
                %< Evaluate LHS mass term
                if isempty(M)
                    Me = speye(length(U0));
                elseif isa(M, 'function_handle')
                    if isempty(obj.MLhs)
                        if nargin(M) == 1
                            obj.MLhs = M(t0 + obj.c(i) * dt);
                        elseif nargin(M) == 2
                            [~, obj.MLhs] = M(t0 + obj.c(i) * dt, t0 + obj.c(i) * dt);
                        end
                    end
                    Me = obj.MLhs;
                else
                    Me = M;
                end

                %< Construct LHS matrix
                obj.setLhs(i, Me);
            end

            %< Evaluate RHS mass term at previous time level
            if isempty(M)
                Me = [];
            elseif isa(M, 'function_handle')
                if nargin(M) == 1
                    if isempty(obj.MRhs{i})
                        obj.MRhs{i} = M(t0);
                    end
                    Me = obj.MRhs{i};
                elseif nargin(M) == 2
                    Me = cell(1, 2);
                    if isempty(obj.MRhs{i}) || obj.timeline.hasStepSizeChanged
                        [Me{1:2}] = M(t0 + obj.c(i) * dt, t0);
                        obj.MRhs{i} = Me{2};
                    else
                        Me{2} = obj.MRhs{i};
                        Me{1} = M(t0 + obj.c(i) * dt, t0);
                    end
                end
            else
                Me = M;
            end

            %< Construct RHS vector
            obj.setRhs(i, Me);

            %< Solve for stage approximation
            U = obj.solver.solve(obj.lhs{i}, obj.rhs{i});

            %< Evaluate nonlinear term at stage approximation
            if isempty(F)
                Fe = [];
            elseif isa(F, 'function_handle')
                Fe = F(U);
            else
                Fe = F;
            end

            %< Compute derivative approximation for this stage
            obj.K{i} = 0;
            if ~isempty(Le)
                obj.K{i} = obj.K{i} + Le * U;
            end
            if ~isempty(Fe)
                obj.K{i} = obj.K{i} + Fe;
            end
            if ~isempty(Se)
                obj.K{i} = obj.K{i} + Se;
            end
        end

        function U = step(obj, L, F, S, M)
            % STEP Advance one time step using the explicit RK method.
            %
            %   U = step(obj, L, F, S, M) performs a complete time step by
            %   evaluating all stages explicitly and combining them according
            %   to the Butcher tableau weights. Empty terms contribute zero
            %   to all stages.
            %
            % Inputs:
            %   obj - The ExrkIntegrator object
            %   L - Linear operator (matrix, function, or empty)
            %   F - Nonlinear operator (function or empty)
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
                [~] = obj.stage(i, L, F, S, M);
            end

            %< Update LHS matrix for final solve if step size changed
            if obj.timeline.hasStepSizeChanged
                %< Evaluate LHS mass term
                if isempty(M)
                    Me = speye(length(U0));
                elseif isa(M, 'function_handle')
                    if isempty(obj.MLhs)
                        if nargin(M) == 1
                            obj.MLhs = M(t0 + dt);
                        elseif nargin(M) == 2
                            [~, obj.MLhs] = M(t0 + dt, t0 + dt);
                        end
                    end
                    Me = obj.MLhs;
                else
                    Me = M;
                end
                obj.lhs{s + 1} = Me;
            end

            %< Evaluate mass term for final solve
            if isempty(M)
                Me = [];
            elseif isa(M, 'function_handle')
                if nargin(M) == 1
                    if isempty(obj.MRhs{s + 1})
                        obj.MRhs{s + 1} = M(t0);
                    end
                    Me = obj.MRhs{s + 1};
                elseif nargin(M) == 2
                    Me = cell(1, 2);
                    if isempty(obj.MRhs{s + 1}) || obj.timeline.hasStepSizeChanged
                        [Me{1:2}] = M(t0 + dt, t0);
                        obj.MRhs{s + 1} = Me{2};
                    else
                        Me{2} = obj.MRhs{s + 1};
                        Me{1} = M(t0 + dt, t0);
                    end
                end
            else
                Me = M;
            end

            %< Add mass term contribution from previous solution
            if isempty(Me)
                obj.rhs{s + 1} = U0;
            elseif iscell(Me)
                obj.rhs{s + 1} = Me{2} * U0;
                if ~isempty(Me{1})
                    obj.rhs{s + 1} = obj.rhs{s + 1} + Me{1};
                end
            else
                obj.rhs{s + 1} = Me * U0;
            end

            %< Add weighted stage derivatives
            for i = 1:s
                obj.rhs{s + 1} = obj.rhs{s + 1} + dt * obj.b(i) * obj.K{i};
            end

            %< Solve for final step approximation
            U = obj.solver.solve(obj.lhs{s + 1}, obj.rhs{s + 1});
        end
    end
end