classdef BdfIntegrator < approx.odeint.OdeIntegrator
    % BDFINTEGRATOR Backward Differentiation Formula integrator base class.
    %
    %   BdfIntegrator implements multi-step Backward Differentiation
    %   Formula (BDF) methods for solving ordinary differential equations
    %   of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + S,
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator, and
    %   \f$S\f$ is a source term. BDF methods use multiple previous time
    %   steps to achieve higher- order accuracy and are A-stable, making
    %   them suitable for stiff problems.
    %
    %   BDF methods are implicit and require solving linear systems at each
    %   time step. They are particularly effective for stiff ODEs where
    %   explicit methods would require prohibitively small time steps.
    %
    %   The integrator also supports augmented systems with constraints by
    %   providing coupling operators (P), constraint matrices (A, B), and
    %   constraint source terms (C). When constraint parameters are provided,
    %   the integrator automatically switches to augmented mode.
    %
    % See also:
    %   approx.odeint.OdeIntegrator, approx.odeint.Bdf2Integrator,
    %   approx.odeint.BeIntegrator

    properties (Access = public)
        alpha % Coefficients for previous solution values (row vector)
        beta  % Coefficient for current derivative evaluation (scalar)
    end

    methods
        function obj = BdfIntegrator(nSteps, final)
            % BDFINTEGRATOR Constructor for BdfIntegrator.
            %
            %   obj = BdfIntegrator(nSteps, final) creates a BDF integrator
            %   with the specified number of time steps and final time.

            arguments
                nSteps {mustBeInteger, mustBePositive}
                final {mustBeNumeric, mustBeNonnegative}
            end

            obj@approx.odeint.OdeIntegrator(nSteps, 1, final);
        end

        function obj = reset(obj, options)
            % RESET Reset the BDF integrator to initial state.
            %
            %   obj = reset(obj) calls the parent reset method and resets
            %   the underlying linear solver.

            arguments
                obj approx.odeint.BdfIntegrator
                options.linOpts struct = struct() % Linear solver options
            end

            reset@approx.odeint.OdeIntegrator(obj);

            args = namedargs2cell(options.linOpts);
            obj.Solver = core.linalg.LinearSolver(args{:});
        end

        function obj = setTimeStep(obj, options)
            % SETTIMESTEP Set time step and coefficients for the BDF
            % integrator.
            %
            %   obj = setTimeStep(obj, h=h, C=C) sets time step as dt = C*h
            %   for dynamic timelines.
            %
            %   obj = setTimeStep(obj, h=h, C=C, p=p) sets time step as dt
            %   = C*h^p for dynamic timelines.
            %
            %   obj = setTimeStep(obj, dt=dt) directly sets the time step
            %   for manual control (use with caution).

            arguments
                obj approx.odeint.BdfIntegrator
                options.dt {mustBeNumeric, mustBePositive} = []
                options.h {mustBeNumeric, mustBePositive} = []
                options.C {mustBeNumeric, mustBePositive} = 1
                options.p {mustBeNumeric, mustBePositive} = 1
            end

            varargin = namedargs2cell(options);
            setTimeStep@approx.odeint.OdeIntegrator(obj, varargin{:});
            obj.setCoefficients();
        end

        function U = step(obj, options)
            % STEP Advance one time step using the BDF method.
            %
            %   U = step(obj, L=L, S=S, M=M) performs a complete time step by
            %   solving the implicit BDF system. For BDF methods, this
            %   is equivalent to computing a single stage.
            %
            %   U = step(obj, L=L, S=S, M=M, P=P, A=A, B=B, C=C) performs a
            %   complete time step for augmented BDF systems with constraints.

            arguments
                obj approx.odeint.BdfIntegrator
                options.L = [] % Linear operator
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

            U = obj.stage(options.L, options.S, options.M, ...
                options.P, options.A, options.B, options.C, dt);
        end
    end

    methods (Access = protected)
        function obj = setLhs(obj, L, M, P, A, B, dt, t)
            % SETLhs Set the left-hand side matrix for the BDF system.

            %< Detect augmented mode
            isAugmented = ~isempty(P) || ~isempty(A) || ~isempty(B);

            %< Evaluate mass operator for LHS
            Me = obj.evalMLhs(M, t);

            %< Evaluate linear operator if provided
            if ~isempty(L)
                Le = obj.evalOp(L, [], t);
            else
                Le = [];
            end

            if isAugmented
                %< Construct augmented block matrix
                Lhs = cell(2, 2);

                %< Block (1,1): M - dt*beta*L
                Lhs{1, 1} = Me;
                if ~isempty(Le)
                    Lhs{1, 1} = Lhs{1, 1} - dt * obj.beta * Le;
                end

                %< Block (1,2): -dt*beta*P
                if ~isempty(P)
                    Pe = obj.evalOp(P, [], t);
                    Lhs{1, 2} = -dt * obj.beta * Pe;
                else
                    Lhs{1, 2} = [];
                end

                %< Block (2,1): A
                Lhs{2, 1} = A;

                %< Block (2,2): B
                Lhs{2, 2} = B;

                Lhs = core.linalg.block(Lhs);
            else
                %< Regular: M - dt*beta*L
                Lhs = Me;
                if ~isempty(Le)
                    Lhs = Lhs - dt * obj.beta * Le;
                end
            end

            obj.Solver.setLhs(Lhs);
            obj.Solver.precompute();
        end

        function obj = setRhs(obj, S, M, C, dt, t)
            % SETRhs Set the right-hand side vector for the BDF system.

            %< Detect augmented mode
            isAugmented = ~isempty(C);

            m = min(obj.Timeline.Count, obj.NSteps);

            %< Initialize with mass term contribution from previous solutions
            Rhs1 = 0;
            for j = 1:m
                if isempty(M)
                    R = obj.U0{j};
                elseif iscell(M)
                    n = size(M{j}{2}, 2);
                    R = M{j}{2} * obj.U0{j}(1:n);
                    if ~isempty(M{j}{1})
                        R = R + M{j}{1}(1:n);
                    end
                else
                    n = size(M, 2);
                    R = M * obj.U0{j}(1:n);
                end
                Rhs1 = Rhs1 - obj.alpha(j) * R;
            end

            %< Add source term contribution
            if ~isempty(S)
                Rhs1 = Rhs1 + dt * obj.beta * S;
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

        function U = stage(obj, L, S, M, P, A, B, C, dt)
            % STAGE Compute a single stage of the BDF method.

            %< Detect augmented mode
            isAugmented = ~isempty(P) || ~isempty(A) || ...
                ~isempty(B) || ~isempty(C);

            t0 = obj.Timeline.Now;
            U0 = obj.U0{1};

            if isAugmented
                %< Evaluate operators for augmented system
                Le = obj.evalOp(L, U0, [t0 + dt, t0]);
                Se = obj.evalOp(S, U0, [t0 + dt, t0]);
                Pe = obj.evalOp(P, U0, [t0 + dt, t0]);
                Ae = obj.evalOp(A, U0, [t0 + dt, t0]);
                Be = obj.evalOp(B, U0, [t0 + dt, t0]);
                Ce = obj.evalOp(C, U0, [t0 + dt, t0]);

                %< Update Lhs matrix if step size changed
                if obj.Timeline.HasStepSizeChanged
                    obj.setLhs(Le, M, Pe, Ae, Be, dt, t0+dt);
                end

                %< Evaluate Rhs mass term for all previous time levels
                if isempty(M)
                    Me = [];
                elseif isa(M, 'function_handle')
                    m = min(obj.Timeline.Count, obj.NSteps);
                    h = obj.Timeline.StepSizeQueue;
                    Me = cell(1, m);
                    for j = 1:m
                        tj = t0 + dt - sum(h(1:j));
                        Me{j} = obj.evalMRhs(j, M, [t0 + dt, tj]);
                    end
                else
                    Me = M;
                end

                %< Construct augmented Rhs vector
                obj.setRhs(Se, Me, Ce, dt, [t0+dt, t0]);
            else
                %< Regular BDF system
                Le = obj.evalOp(L, U0, [t0 + dt, t0]);
                Se = obj.evalOp(S, U0, [t0 + dt, t0]);

                %< Update Lhs matrix if step size changed
                if obj.Timeline.HasStepSizeChanged
                    obj.setLhs(Le, M, [], [], [], dt, t0+dt);
                end

                %< Evaluate Rhs mass term for all previous time levels
                if isempty(M)
                    Me = [];
                elseif isa(M, 'function_handle')
                    m = min(obj.Timeline.Count, obj.NSteps);
                    h = obj.Timeline.StepSizeQueue;
                    Me = cell(1, m);
                    for j = 1:m
                        tj = t0 + dt - sum(h(1:j));
                        Me{j} = obj.evalMRhs(j, M, [t0 + dt, tj]);
                    end
                else
                    Me = M;
                end

                %< Construct Rhs vector
                obj.setRhs(Se, Me, [], dt, [t0+dt, t0]);
            end

            %< Solve the linear system using stateful API
            U = obj.Solver.solve();
        end
    end

    methods (Abstract)
        % SETCOEFFICIENTS Set BDF coefficients for variable step sizes.
        obj = setCoefficients(obj)
    end
end
