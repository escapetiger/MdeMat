classdef BdfIntegrator < mdemat.common.ode.OdeIntegrator
    % BDFINTEGRATOR Backward Differentiation Formula integrator base class.
    %
    %   BdfIntegrator implements multi-step Backward Differentiation
    %   Formula (BDF) methods for solving ordinary differential equations
    %   of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + F(u,t) + S,
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator, and
    %   \f$S\f$ is a source term. BDF methods use multiple previous time
    %   steps to achieve higher-order accuracy. BDF2 is A-stable; higher
    %   orders have more restricted stability regions.
    %
    %   BDF methods solve the current nonlinear stage implicitly, with an
    %   independently checked JFNK residual. They are useful for stiff ODEs where
    %   explicit methods would require prohibitively small time steps.
    %
    %   The integrator also supports augmented systems with constraints by
    %   providing coupling operators (P), constraint matrices (A, B), and
    %   constraint source terms (C). When constraint parameters are provided,
    %   the integrator automatically switches to augmented mode.
    %
    % See also:
    %   mdemat.common.ode.OdeIntegrator, mdemat.common.ode.Bdf2Integrator,
    %   mdemat.common.ode.BeIntegrator

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

            obj@mdemat.common.ode.OdeIntegrator(nSteps, 1, final);
        end

        function obj = reset(obj, options)
            % RESET Reset the BDF integrator to initial state.
            %
            %   obj = reset(obj) calls the parent reset method and resets
            %   the underlying linear solver.

            arguments
                obj mdemat.common.ode.BdfIntegrator
                options.linOpts struct = struct() % Linear solver options
            end

            reset@mdemat.common.ode.OdeIntegrator(obj);

            args = namedargs2cell(options.linOpts);
            obj.Solver = mdemat.common.linalg.LinearSolver(args{:});
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
                obj mdemat.common.ode.BdfIntegrator
                options.dt {mustBeNumeric, mustBePositive} = []
                options.h {mustBeNumeric, mustBePositive} = []
                options.C {mustBeNumeric, mustBePositive} = 1
                options.p {mustBeNumeric, mustBePositive} = 1
            end

            varargin = namedargs2cell(options);
            setTimeStep@mdemat.common.ode.OdeIntegrator(obj, varargin{:});
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
                obj mdemat.common.ode.BdfIntegrator
                options.L = [] % Linear operator
                options.F = [] % Nonlinear current-state term (implicit)
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

            options.M = obj.completeConstraintMass(options.M, options.P, options.A, options.B, options.C);

            if ~isempty(options.dt)
                dt = options.dt;
                obj.Timeline.OldStepSizeQueue(:) = 0;
            else
                dt = obj.Timeline.StepSize;
            end

            U = obj.stage(options.L, options.F, options.S, options.M, ...
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

                Lhs = mdemat.common.linalg.block(Lhs);
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

        function obj = setRhs(obj, S, M, C, dt, ~)
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
                    entry = M{j};
                    if iscell(entry)
                        n = size(entry{2}, 2);
                        R = entry{2} * obj.U0{j}(1:n);
                        if ~isempty(entry{1}), R = R+entry{1}(1:n); end
                    elseif isempty(entry)
                        R = obj.U0{j};
                    else
                        n = size(entry, 2);
                        R = entry*obj.U0{j}(1:n);
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
                Rhs = mdemat.common.linalg.block(Rhs);
            else
                %< Regular RHS
                Rhs = Rhs1;
            end

            obj.Solver.setRhs(Rhs);
        end

        function U = stage(obj, L, F, S, M, P, A, B, C, dt)
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
                Fe = obj.evalOp(F, U0, [t0 + dt, t0]);
                if ~isempty(Fe)
                    if isempty(Se), Se = Fe; else, Se = Se+Fe; end
                end
                Pe = obj.evalOp(P, U0, [t0 + dt, t0]);
                Ae = obj.evalOp(A, U0, [t0 + dt, t0]);
                Be = obj.evalOp(B, U0, [t0 + dt, t0]);
                Ce = obj.evalOp(C, U0, [t0 + dt, t0]);

                %< Update Lhs matrix if step size changed
                if obj.Timeline.HasStepSizeChanged || isa(L, 'function_handle') || ...
                        isa(M, 'function_handle') || isa(P, 'function_handle')
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
                Fe = obj.evalOp(F, U0, [t0 + dt, t0]);
                if ~isempty(Fe)
                    if isempty(Se), Se = Fe; else, Se = Se+Fe; end
                end

                %< Update Lhs matrix if step size changed
                if obj.Timeline.HasStepSizeChanged || isa(L, 'function_handle') || isa(M, 'function_handle')
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
            if ~isempty(F) || isa(L, 'function_handle') || isa(S, 'function_handle') || isa(P, 'function_handle')
                mass = obj.evalMLhs(M, t0+dt);
                if isAugmented, prescribed = Ce; else, prescribed = []; end
                obj.setRhs([], Me, prescribed, dt, [t0+dt, t0]);
                base = obj.Solver.Rhs(1:size(mass, 1));
                residual = @(candidate) obj.implicitResidual(candidate, mass, base, ...
                    dt*obj.beta, L, F, S, P, A, B, C, [t0+dt, t0]);
                if ~isempty(F), U = U0; end
                U = obj.solveImplicitStage(residual, U);
            end
        end
    end

    methods (Abstract)
        % SETCOEFFICIENTS Set BDF coefficients for variable step sizes.
        obj = setCoefficients(obj)
    end
end
