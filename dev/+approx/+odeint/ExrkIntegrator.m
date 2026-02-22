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
        end
        
        function obj = setLhs(obj, i, options)
            % SETLhs Set the mass matrix for explicit solve at stage i.
            %
            %   obj = setLhs(obj, i, M=M) sets the left-hand side matrix to
            %   the mass operator for stage i. For explicit methods, this
            %   is used to solve M*K = Rhs for the derivative
            %   approximation.

            arguments
                obj approx.odeint.ExrkIntegrator
                i {mustBePositive, mustBeInteger} % Stage index
                options.M = [] % Mass operator
                options.dt = [] % Time step override
                options.t = [] % Time values for time-dependent operators
            end

            %< Evaluate mass operator for LHS
            Me = obj.evalMLhs(options.M, options.t);

            obj.Solver{i}.setLhs(Me);
            obj.Solver{i}.precompute();
        end
        
        function obj = setRhs(obj, i, options)
            % SETRhs Set the right-hand side vector for stage i.
            %
            %   obj = setRhs(obj, i, M=M) constructs the right-hand side
            %   vector using mass operator and accumulated previous stages.

            arguments
                obj approx.odeint.ExrkIntegrator
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

            %< Add accumulated term from previous stages
            for j = 1:(i - 1)
                Rhs = Rhs + dt * obj.A(i, j) * obj.K{j};
            end

            obj.Solver{i}.setRhs(Rhs);
        end
        
        function U = stage(obj, i, options)
            % STAGE Compute stage i of the explicit RK method.
            %
            %   U = stage(obj, i, L=L, F=F, S=S, M=M) computes the solution at stage i
            %   by evaluating all terms explicitly and computing the derivative
            %   approximation for use in subsequent stages.

            arguments
                obj approx.odeint.ExrkIntegrator
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

            %< Evaluate linear term at initial time
            Le = obj.evalOp(options.L, U0, [t0, t0]);

            %< Evaluate source term at stage time
            Se = obj.evalOp(options.S, U0, [t0+obj.c(i)*dt, t0]);

            %< Update Lhs matrix if step size changed or custom dt provided
            if obj.Timeline.HasStepSizeChanged || ~isempty(options.dt)
                %< Construct Lhs matrix
                obj.setLhs(i, M=options.M, dt=dt, t=t0+obj.c(i)*dt);
            end

            %< Construct Rhs vector
            obj.setRhs(i, M=options.M, dt=dt, t=[t0+obj.c(i)*dt, t0]);

            %< Solve for stage approximation using stateful API
            U = obj.Solver{i}.solve();

            %< Evaluate nonlinear term at stage approximation
            Fe = obj.evalOp(options.F, U, [t0+obj.c(i)*dt, t0]);

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
        
        function U = step(obj, options)
            % STEP Advance one time step using the explicit RK method.
            %
            %   U = step(obj, L=L, F=F, S=S, M=M) performs a complete time step by
            %   evaluating all stages explicitly and combining them according
            %   to the Butcher tableau weights. Empty terms contribute zero
            %   to all stages.

            arguments
                obj approx.odeint.ExrkIntegrator
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
            s = obj.NStages;

            %< Compute all stage approximations
            for i = 1:s
                [~] = obj.stage(i, L=options.L, F=options.F, S=options.S, M=options.M, dt=dt);
            end

            %< Update Lhs matrix for final solve if step size changed or custom dt provided
            if obj.Timeline.HasStepSizeChanged || ~isempty(options.dt)
                %< Evaluate Lhs mass term
                Me = obj.evalMLhs(options.M, t0+dt);

                %< Construct Lhs matrix
                obj.setLhs(s+1, M=Me, dt=dt, t=t0+dt);
            end

            %< Evaluate mass term for final solve
            Me = obj.evalMRhs(s+1, options.M, [t0+dt, t0]);

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

            %< Add weighted stage derivatives
            for i = 1:s
                Rhs = Rhs + dt * obj.b(i) * obj.K{i};
            end

            obj.Solver{s+1}.setRhs(Rhs);

            %< Solve for final step approximation using stateful API
            U = obj.Solver{s + 1}.solve();
        end
    end
end