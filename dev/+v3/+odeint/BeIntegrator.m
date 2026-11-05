classdef BeIntegrator < approx.odeint.OdeIntegrator
    % BEINTEGRATOR Backward Euler integrator.
    %
    %   BeIntegrator implements the Backward Euler method, a first-order
    %   implicit time integration scheme for solving ordinary differential
    %   equations of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + S,
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator, and
    %   \f$S\f$ is a source term. Any of these terms can be empty. For mass
    %   operator, it means it reduces to an identity. For linear and source
    %   term, it means they vanish from the equation. The method is
    %   A-stable and suitable for stiff problems.
    %
    %   The Backward Euler method is the simplest implicit method,
    %   requiring only one function evaluation per time step but involving
    %   the solution of a linear system at each step.
    %
    % Examples:
    %   % Basic usage
    %   integrator = approx.odeint.BeIntegrator(1.0);
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
    %   approx.odeint.OdeIntegrator, approx.odeint.BdfIntegrator,
    %   approx.odeint.FeIntegrator

    properties (Constant)
        ORDER = 1 % Accuracy order of the Backward Euler method
    end

    methods
        function obj = BeIntegrator(final)
            % BEINTEGRATOR Constructor for BeIntegrator.
            %
            %   obj = BeIntegrator(final) creates a Backward Euler
            %   integrator with the specified final time.
            %
            % Inputs:
            %   final - Final time for integration (positive scalar)
            %
            % Outputs:
            %   obj - Constructed BeIntegrator object

            obj@approx.odeint.OdeIntegrator(1, 1, final);
        end
       
        function obj = setLhs(obj, L, M)
            % SETLHS Set the left-hand side matrix for the Backward Euler
            % system.
            %
            %   obj = setLhs(obj, L, M) constructs the left-hand side
            %   matrix as M - dt*L for the implicit Backward Euler system.
            %
            % Inputs:
            %   obj - The BeIntegrator object
            %   L - Linear operator (matrix or empty)
            %   M - Mass operator (matrix or empty)
            %
            % Outputs:
            %   obj - The BeIntegrator object

            dt = obj.timeline.dt;

            %< Add mass term
            obj.lhs = M;

            %< Add linear term contribution
            if ~isempty(L)
                obj.lhs = obj.lhs - dt * L;
            end
        end

        function obj = setRhs(obj, S, M)
            % SETRHS Set the right-hand side vector for the Backward Euler
            % system.
            %
            %   obj = setRhs(obj, S, M) constructs the right-hand side
            %   vector using the mass operator and source term.
            %
            % Inputs:
            %   obj - The BeIntegrator object
            %   S - Source term (vector, function, or empty)
            %   M - Mass operator (matrix, function, or empty)
            %
            % Outputs:
            %   obj - The BeIntegrator object

            dt = obj.timeline.dt;
            U0 = obj.U0{1};

            %< Add mass term contribution from previous solution
            if isempty(M)
                obj.rhs = U0;
            elseif iscell(M)
                obj.rhs = M{2} * U0;
                if ~isempty(M{1})
                    obj.rhs = obj.rhs + M{1};
                end
            else
                obj.rhs = M * U0;
            end

            %< Add source term contribution
            if ~isempty(S)
                obj.rhs = obj.rhs + dt * S;
            end
        end

        function U = stage(obj, L, S, M)
            % STAGE Compute the single stage of the Backward Euler method.
            %
            %   U = stage(obj, L, S, M) computes the solution at the
            %   current time step by solving the implicit Backward Euler
            %   system. For Backward Euler, there is only one stage per
            %   time step.
            %
            % Inputs:
            %   obj - The BeIntegrator object
            %   L - Linear operator (matrix, function, or empty)
            %   S - Source term (vector, function, or empty)
            %   M - Mass operator (matrix, function, or empty)
            %
            % Outputs:
            %   U - Solution approximation at current stage

            t0 = obj.timeline.now;
            dt = obj.timeline.dt;
            U0 = obj.U0{1};

            %< Evaluate linear term at current time
            if isempty(L)
                Le = [];
            elseif isa(L, 'function_handle')
                Le = L(t0 + dt);
            else
                Le = L;
            end

            %< Evaluate source term at current time
            if isempty(S)
                Se = [];
            elseif isa(S, 'function_handle')
                Se = S(t0 + dt);
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
                            obj.MLhs = M(t0 + dt);
                        elseif nargin(M) == 2
                            [~, obj.MLhs] = M(t0 + dt, t0 + dt);
                        end
                    end
                    Me = obj.MLhs;
                else
                    Me = M;
                end

                %< Construct LHS matrix
                obj.setLhs(Le, Me);
            end

            %< Evaluate RHS mass term at previous time level
            if isempty(M)
                Me = [];
            elseif isa(M, 'function_handle')
                if nargin(M) == 1
                    if isempty(obj.MRhs)
                        obj.MRhs = M(t0);
                    end
                    Me = obj.MRhs;
                elseif nargin(M) == 2
                    Me = cell(1, 2);
                    if isempty(obj.MRhs) || obj.timeline.hasStepSizeChanged                        
                        [Me{1:2}] = M(t0 + dt, t0);
                        obj.MRhs = Me{2};
                    else
                        Me{2} = obj.MRhs;
                        Me{1} = M(t0 + dt, t0);
                    end
                end
            else
                Me = M;
            end

            %< Construct RHS vector
            obj.setRhs(Se, Me);
            
            %< Solve the linear system
            U = obj.solver.solve(obj.lhs, obj.rhs);
        end

        function U = step(obj, L, S, M)
            % STEP Advance one time step using the Backward Euler method.
            %
            %   U = step(obj, L, S, M) performs a complete time step by
            %   solving the implicit Backward Euler system. For Backward
            %   Euler, this is equivalent to computing a single stage.
            %
            % Inputs:
            %   obj - The BeIntegrator object
            %   L - Linear operator (matrix, function, or empty)
            %   S - Source term (vector, function, or empty)
            %   M - Mass operator (matrix, function, or empty)
            %
            % Outputs:
            %   U - Solution at the next time step

            U = obj.stage(L, S, M);
        end
    end
end