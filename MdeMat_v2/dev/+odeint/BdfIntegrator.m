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
    % Examples:
    %   % Cannot instantiate abstract class directly
    %   % Use concrete subclasses instead:
    %   integrator = approx.odeint.Bdf2Integrator(1.0);
    %   
    %   % Full implicit problem with all terms
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
            %
            % Inputs:
            %   nSteps - Number of time steps for the BDF method (positive integer)
            %   final - Final time for integration (positive scalar)
            %
            % Outputs:
            %   obj - Constructed BdfIntegrator object

            obj@approx.odeint.OdeIntegrator(nSteps, 1, final);
        end
        
        function obj = setLhs(obj, L, M)
            % SETLHS Set the left-hand side matrix for the BDF system.
            %
            %   obj = setLhs(obj, L, M) constructs the left-hand side
            %   matrix as M - dt*beta*L for the implicit BDF system.
            %
            % Inputs:
            %   obj - The BdfIntegrator object
            %   L - Linear operator (matrix or empty)
            %   M - Mass operator (matrix or empty)
            %
            % Outputs:
            %   obj - The BdfIntegrator object

            dt = obj.timeline.dt;

            %< Start with mass term
            obj.lhs = M;

            %< Add linear term contribution
            if ~isempty(L)
                obj.lhs = obj.lhs - dt * obj.beta * L;
            end
        end

        function obj = setRhs(obj, S, M)
            % SETRHS Set the right-hand side vector for the BDF system.
            %
            %   obj = setRhs(obj, S, M) constructs the right-hand side
            %   vector using the BDF approximation and source term.
            %
            % Inputs:
            %   obj - The BdfIntegrator object
            %   S - Source term (vector, function, or empty)
            %   M - Mass operator (matrix, function, or empty)
            %
            % Outputs:
            %   obj - The BdfIntegrator object

            dt = obj.timeline.dt;
            m = min(obj.timeline.count, obj.nSteps);

            %< Initialize with mass term contribution from previous
            %< solutions
            obj.rhs = 0;
            for j = 1:m
                if isempty(M)
                    R = obj.U0{j};
                elseif iscell(M)
                    R = M{j}{2} * obj.U0{j};
                    if ~isempty(M{j}{1})
                        R = R + M{j}{1};
                    end
                else
                    R = M * obj.U0{j};
                end
                obj.rhs = obj.rhs - obj.alpha(j) * R;
            end

            %< Add source term contribution
            if ~isempty(S)
                obj.rhs = obj.rhs + dt * obj.beta * S;
            end
        end

        function U = stage(obj, L, S, M)
            % STAGE Compute a single stage of the BDF method.
            %
            %   U = stage(obj, L, S, M) computes the solution at the current
            %   time step by solving the implicit BDF system.
            %
            % Inputs:
            %   obj - The BdfIntegrator object
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
            Le = obj.evaluateOperator(L, U0, [t0 + dt, t0]);

            %< Evaluate source term at current time
            Se = obj.evaluateOperator(S, U0, [t0 + dt, t0]);

            %< Update LHS matrix if step size changed
            if obj.timeline.hasStepSizeChanged
                %< Evaluate LHS mass term
                Me = obj.evaluateLhsMass(M, t0+dt);

                %< Construct LHS matrix
                obj.setLhs(Le, Me);
            end

            %< Evaluate RHS mass term for all previous time levels
            if isempty(M)
                Me = [];
            elseif isa(M, 'function_handle')
                m = min(obj.timeline.count, obj.nSteps);
                Me = cell(1, m);
                for j = 1:m
                    tj = t0 + dt - sum(obj.timeline.h(1:j));
                    Me{j} = obj.evaluateRhsMass(j, M, [t0 + dt, tj]);
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
            % STEP Advance one time step using the BDF method.
            %
            %   U = step(obj, L, S, M) performs a complete time step by
            %   solving the implicit BDF system. For BDF methods, this
            %   is equivalent to computing a single stage.
            %
            % Inputs:
            %   obj - The BdfIntegrator object
            %   L - Linear operator (matrix, function, or empty)
            %   S - Source term (vector, function, or empty)
            %   M - Mass operator (matrix, function, or empty)
            %
            % Outputs:
            %   U - Solution at the next time step

            U = obj.stage(L, S, M);
        end
    end

    methods (Abstract)
        % SETCOEFFICIENTS Set BDF coefficients for variable step sizes.
        obj = setCoefficients(obj)
    end
end