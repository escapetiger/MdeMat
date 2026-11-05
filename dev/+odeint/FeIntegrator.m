classdef FeIntegrator < approx.odeint.OdeIntegrator
    % FEINTEGRATOR Forward Euler integrator.
    %
    %   FeIntegrator implements the Forward Euler method, a first-order
    %   explicit time integration scheme for solving ordinary differential
    %   equations of the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + F(u) + S(t),
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator,
    %   \f$F(u)\f$ is a nonlinear operator, and \f$S(t)\f$ is a source
    %   term. Any of these terms can be empty. For mass term, it means it
    %   reduces to an identity. For linear, nonlinear and source term, it
    %   means they vanish from the equation. This is the simplest explicit
    %   method, using only one evaluation per time step.
    %
    %   Forward Euler has a stability region of \f$|1 + z| \leq 1\f$, which
    %   restricts the time step for stability. For eigenvalue
    %   \f$\lambda\f$, the stability condition is \f$|1 + \lambda \Delta t|
    %   \leq 1\f$.
    %
    % Examples:
    %   % Basic usage with all terms
    %   integrator = approx.odeint.FeIntegrator(1.0);
    %   U_new = integrator.step(L, F, S, M);
    %   integrator.update(U_new);
    %   integrator.advance();
    %
    %   % Pure diffusion (only linear term)
    %   U_new = integrator.step(L, [], [], []);
    %
    %   % Pure advection (only nonlinear term)
    %   U_new = integrator.step([], F, [], []);
    %
    %   % Forced system (linear + source)
    %   U_new = integrator.step(L, [], S, []);
    %
    %   % With mass matrix
    %   U_new = integrator.step(L, F, S, M);
    %
    % Notes:
    %   Forward Euler is conditionally stable and may require very small
    %   time steps for stiff problems. Use implicit methods for better
    %   stability properties with stiff systems.
    %
    % See also:
    %   approx.odeint.OdeIntegrator, approx.odeint.BeIntegrator,
    %   approx.odeint.ExrkIntegrator

    properties (Constant)
        ORDER = 1 % Accuracy order of the Forward Euler method
    end

    methods
        function obj = FeIntegrator(final)
            % FEINTEGRATOR Constructor for FeIntegrator.
            %
            %   obj = FeIntegrator(final) creates a Forward Euler integrator
            %   with the specified final time.
            %
            % Inputs:
            %   final - Final time for integration (positive scalar)
            %
            % Outputs:
            %   obj - Constructed FeIntegrator object
            
            obj@approx.odeint.OdeIntegrator(1, 1, final);
        end

        function obj = setLhs(obj, M)
            % SETLHS Set the left-hand side matrix for the Forward Euler
            % system.
            %
            %   obj = setLhs(obj, M) sets the left-hand side matrix to the
            %   mass operator.
            %
            % Inputs:
            %   obj - The FeIntegrator object
            %   M - Mass operator (matrix or empty)
            %
            % Outputs:
            %   obj - The FeIntegrator object

            obj.lhs = M;
        end
        
        function obj = setRhs(obj, L, F, S, M)
            % SETRHS Set the right-hand side vector for the Forward Euler
            % system.
            %
            %   obj = setRhs(obj, L, F, S, M) constructs the right-hand
            %   side vector using the explicit Forward Euler formula: M*U0
            %   + dt*(L*U0 + F + S).
            %
            % Inputs:
            %   obj - The FeIntegrator object
            %   L - Linear operator (matrix, function, or empty)
            %   F - Nonlinear term (function, vector, or empty)
            %   S - Source term (vector, function, or empty)
            %   M - Mass operator (matrix, function, or empty)
            %
            % Outputs:
            %   obj - The FeIntegrator object

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

            %< Add linear term contribution: dt*L*U0
            if ~isempty(L)
                obj.rhs = obj.rhs + dt * L * U0;
            end
            
            %< Add nonlinear term contribution: dt*F(U0)
            if ~isempty(F)
                obj.rhs = obj.rhs + dt * F;
            end
            
            %< Add source term contribution: dt*S(t0)
            if ~isempty(S)
                obj.rhs = obj.rhs + dt * S;
            end
        end

        function U = stage(obj, L, F, S, M)
            % STAGE Compute the single stage of Forward Euler method.
            %
            %   U = stage(obj, L, F, S, M) computes the solution at the next
            %   time step using the explicit Forward Euler formula. For
            %   Forward Euler, there is only one stage per time step.
            %
            % Inputs:
            %   obj - The FeIntegrator object
            %   L - Linear operator (matrix, function, or empty)
            %   F - Nonlinear operator (function or empty)
            %   S - Source term (vector, function, or empty)
            %   M - Mass operator (matrix, function, or empty)
            %
            % Outputs:
            %   U - Solution approximation at current stage
            
            t0 = obj.timeline.now;
            dt = obj.timeline.dt;
            U0 = obj.U0{1};

            %< Evaluate linear term at initial time
            Le = obj.evaluateOperator(L, U0, t0);

            %< Evaluate nonlinear term at current solution
            Fe = obj.evaluateOperator(F, U0, t0);

            %< Evaluate source term at current time
            Se = obj.evaluateOperator(S, U0, t0);

            %< Update LHS matrix if step size changed
            if obj.timeline.hasStepSizeChanged
                %< Evaluate LHS mass term
                Me = obj.evaluateLhsMass(M, t0 + dt);

                %< Construct LHS matrix
                obj.setLhs(Me);
            end

            %< Evaluate RHS mass term at current time level
            Me = obj.evaluateRhsMass([], M, [t0+dt, t0]);

            %< Construct RHS vector
            obj.setRhs(Le, Fe, Se, Me);
            
            %< Solve the linear system
            U = obj.solver.solve(obj.lhs, obj.rhs);
        end

        function U = step(obj, L, F, S, M)
            % STEP Advance one time step using Forward Euler method.
            %
            %   U = step(obj, L, F, S, M) performs a complete time step
            %   using the explicit Forward Euler method. For Forward Euler,
            %   this is equivalent to computing a single stage.
            %
            % Inputs:
            %   obj - The FeIntegrator object
            %   L - Linear operator (matrix, function, or empty)
            %   F - Nonlinear operator (function or empty)
            %   S - Source term (vector, function, or empty)
            %   M - Mass operator (matrix, function, or empty)
            %
            % Outputs:
            %   U - Solution at the next time step

            U = obj.stage(L, F, S, M);
        end
    end
end