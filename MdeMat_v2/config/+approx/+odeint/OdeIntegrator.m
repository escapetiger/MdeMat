classdef OdeIntegrator < handle
    % ODEINTEGRATOR Abstract base class for ODE integrators.
    %
    %   OdeIntegrator provides a common interface and shared functionality
    %   for all ODE integration methods. This base class manages solution
    %   history, timeline operations, and linear system solving
    %   infrastructure. Concrete integrators must implement the abstract
    %   methods for specific integration schemes.
    %
    %   The class focuses on solving ordinary differential equations of
    %   the form:
    %
    %   \f[
    %     \frac{\mathrm{d}(Mu)}{\mathrm{d}t} = Lu + F(u) + S(t),
    %   \f]
    %
    %   where \f$M\f$ is a mass operator, \f$L\f$ is a linear operator,
    %   \f$F(u)\f$ is a nonlinear operator, and \f$S(t)\f$ is a source term.
    %   Any of these terms can be empty. In particular, M = [] means the
    %   mass operator reduces to an identity operator; {L, F, S} = [] means
    %   they vanish from the equation.
    %
    %   The contribution of mass operator M is determined by the physical
    %   simulation framework:
    %   * In the Eulerian framework, M is space-only so that it can be
    %     represented by a unique static matrix.
    %   * In the Semi-Lagrangian framework, M is time-dependent. The LHS
    %     mass term is forced to be a unique static matrix as in Eulerian
    %     framework. However, the RHS mass term should be updated as time
    %     evolves. In particular, M is a matrix under periodic boundary
    %     conditions, while M is a pair formed by a matrix and a vector
    %     under Dirichlet boundary conditions.
    %
    % See also:
    %   approx.odeint.BeIntegrator, approx.odeint.BdfIntegrator,
    %   approx.odeint.DirkIntegrator, approx.odeint.ExrkIntegrator

    properties
        nSteps % Number of time steps required for the method
        nStages % Number of stages required for the method
        timeline % Timeline object for time step management
        U0 % Cell array of solution history
        lhs % Left-hand side matrix/matrices for linear systems
        rhs % Right-hand side vector/vectors for linear systems
        MLhs % Left-hand side mass term
        MRhs % Right-hand side mass term(s)
        solver % Linear solver for all integrators
    end

    methods
        function obj = OdeIntegrator(nSteps, nStages, final)
            % ODEINTEGRATOR Constructor for OdeIntegrator.
            %
            %   obj = OdeIntegrator(nSteps, nStages, final) creates an ODE
            %   integrator with the specified parameters. This constructor
            %   is called by concrete subclasses to initialize the base
            %   class functionality.
            %
            % Inputs:
            %   nSteps - Number of time steps required for the method (positive integer)
            %   nStages - Number of stages required for the method (positive integer)
            %   final - Final time for integration (positive scalar)
            %
            % Outputs:
            %   obj - Constructed OdeIntegrator object

            obj.nSteps = nSteps;
            obj.nStages = nStages;
            obj.timeline = approx.mesh.DynamicTimeline(final, nSteps);
            obj.reset();
        end

        function obj = reset(obj)
            % RESET Reset the integrator to initial state.
            %
            %   obj = reset(obj) clears solution history, resets the
            %   timeline, and initializes storage arrays for left-hand side
            %   matrices and right-hand side vectors. Also creates a fresh
            %   linear solver instance. This method is called automatically
            %   during construction and can be called to restart
            %   integration.
            %
            % Inputs:
            %   obj - The OdeIntegrator object
            %
            % Outputs:
            %   obj - The OdeIntegrator object

            obj.U0 = cell(obj.nStages, obj.nSteps);
            obj.timeline.reset();

            %< Initialize storage based on method characteristics
            if obj.nStages > 1 && obj.nSteps == 1
                %< Multi-stage, single-step methods (RK-type)
                obj.lhs = cell(obj.nStages+1, obj.nSteps);
                obj.rhs = cell(obj.nStages+1, obj.nSteps);
            elseif obj.nStages == 1 && obj.nSteps > 1
                %< Single-stage, multi-step methods (BDF-type)
                obj.lhs = cell(1, obj.nSteps);
                obj.rhs = cell(1, obj.nSteps);
            else
                %< Single-stage, single-step methods
                obj.lhs = [];
                obj.rhs = [];
            end

            %< Initialize mass term storage
            if obj.nStages > 1 && obj.nSteps == 1
                obj.MLhs = [];
                obj.MRhs = cell(obj.nStages+1, 1);
            elseif obj.nStages == 1 && obj.nSteps > 1
                obj.MLhs = [];
                obj.MRhs = cell(1, obj.nSteps);
            else
                obj.MLhs = [];
                obj.MRhs = [];
            end

            %< Create linear solver instance
            obj.solver = core.linalg.LinearSolver();
        end

        function obj = update(obj, U)
            % UPDATE Update the solution history with a new solution.
            %
            %   obj = update(obj, U) shifts the existing solution history
            %   and adds the new solution at the most recent position. For
            %   multi-step methods, this maintains the required number of
            %   previous solutions. For single-step methods, simply stores
            %   the current solution.
            %
            % Inputs:
            %   obj - The OdeIntegrator object
            %   U - The new solution to add to history
            %
            % Outputs:
            %   obj - The OdeIntegrator object

            %< Shift history for multi-step methods
            if obj.nSteps > 1
                obj.U0 = circshift(obj.U0, [0, 1]);
            end

            %< Store new solution at all stage positions
            for i = 1:obj.nStages
                obj.U0{i, 1} = U;
            end
        end

        function obj = advance(obj)
            % ADVANCE Advance the timeline to the next time step.
            %
            %   obj = advance(obj) updates the timeline's internal time
            %   counter and step size information. Should be called after
            %   each successful integration step to maintain proper time
            %   synchronization.
            %
            % Inputs:
            %   obj - The OdeIntegrator object
            %
            % Outputs:
            %   obj - The OdeIntegrator object

            obj.timeline.advance();
        end

        function obj = setTimeStep(obj, varargin)
            % SETTIMESTEP Set time step for the integrator timeline.
            %
            %   obj = setTimeStep(obj, h, C) sets time step as dt = C*h
            %   for dynamic timelines.
            %
            %   obj = setTimeStep(obj, h, C, p) sets time step as dt =
            %   C*h^p for dynamic timelines.
            %
            %   obj = setTimeStep(obj, dt) directly sets the time step for
            %   manual control (use with caution).
            %
            % Inputs:
            %   obj - The OdeIntegrator object
            %   varargin - Input arguments
            %
            % Outputs:
            %   obj - The OdeIntegrator object

            if isa(obj.timeline, 'approx.mesh.DynamicTimeline')
                obj.timeline.setTimeStep(varargin{:});
            elseif nargin == 2 && isscalar(varargin{1})
                dt = varargin{1};
                core.except.assert(dt >= 0, 'InvalidInput', ...
                    'Time step must be nonnegative.');

                if isa(obj.timeline, 'approx.mesh.DynamicTimeline')
                    obj.timeline.h(1) = dt;
                    obj.timeline.next = obj.timeline.now + dt;
                else
                    core.except.warning('StaticTimeline', ...
                        'Cannot modify time step for static timeline.');
                end
            else
                core.except.error('InvalidInput', ...
                    'Invalid arguments for setTimeStep method.');
            end
        end
    end

    methods (Access = protected)
        function Me = evaluateLhsMass(obj, M, t)
            % EVALUATELHSMASS Evaluate mass operator for left-hand side.
            %
            %   Me = evaluateLhsMass(obj, M, t) evaluates the mass operator
            %   for use in the left-hand side matrix construction. Handles
            %   empty, constant, and time-dependent mass operators.
            %
            % Inputs:
            %   obj - The OdeIntegrator object
            %   M - Mass operator (empty, matrix, or function handle)
            %   t - Time(s) for evaluation (scalar or vector)
            %       Scalar: single evaluation time
            %       Vector: [t1, t2, ...] uses first element for LHS
            %
            % Outputs:
            %   Me - Evaluated mass operator (matrix)

            core.except.assert(isscalar(t), 'InvalidInput', ...
                'LHS mass requires a scalar time.');

            if isempty(M)
                Me = speye(length(obj.U0{1}));
            elseif isa(M, 'function_handle')
                switch nargin(M)
                    case 1
                        obj.MLhs = M(t);
                    case 2
                        [~, obj.MLhs] = M(t, t);
                    otherwise
                        core.except.assert(0, 'InvalidInput', ...
                            'Mass function requires 1 or 2 arguments.');
                end
                Me = obj.MLhs;
            else
                Me = M;
            end
        end

        function Me = evaluateRhsMass(obj, i, M, t)
            % EVALUATERHSMASS Evaluate mass operator for right-hand side.
            %
            %   Me = evaluateRhsMass(obj, i, M, t) evaluates the mass
            %   operator for use in right-hand side vector construction.
            %   Handles caching for efficiency.
            %
            % Inputs:
            %   obj - The OdeIntegrator object
            %   i - Index for caching (stage or step)
            %   M - Mass operator (empty, matrix, or function handle)
            %   t - Time(s) for evaluation (scalar or vector)
            %       Scalar: single time for evaluation
            %       Vector: [tCurrent, tPrevious, ...] for multi-time evaluation
            %
            % Outputs:
            %   Me - Evaluated mass operator (matrix, cell, or empty)

            if isempty(M)
                Me = [];
            elseif isa(M, 'function_handle')
                switch nargin(M)
                    case 1
                        if isscalar(t)
                            te = t;
                        else
                            te = t(2);
                        end
                        if obj.shouldUpdateRhsMass(i)
                            obj.setRhsMass(i, M(te));
                        end
                        Me = obj.getRhsMass(i);
                    case 2
                        core.except.assert(length(t) >= 2, 'InvalidInput', ...
                            'Mass function requires at least 2 time values.');
                        Me = cell(1, 2);
                        if obj.shouldUpdateRhsMass(i)
                            [Me{1:2}] = M(t(1), t(2));
                            obj.setRhsMass(i, Me{2});
                        else
                            Me{2} = obj.getRhsMass(i);
                            Me{1} = M(t(1), t(2));
                        end
                    otherwise
                        core.except.assert(0, 'InvalidInput', ...
                            'Mass function requires 1 or 2 arguments.');
                end
            else
                Me = M;
            end
        end

        function Oe = evaluateOperator(~, O, U, t)
            % EVALUATEOPERATOR Evaluate a general operator (L, F, or S).
            %
            %   Oe = evaluateOperator(obj, O, U, t) evaluates linear,
            %   nonlinear, or source operators with flexible argument
            %   handling.
            %
            % Inputs:
            %   obj - The OdeIntegrator object
            %   O - Operator to evaluate (empty, matrix, or function handle)
            %   U - Current solution vector
            %   t - Time(s) for evaluation (scalar or vector)
            %       Scalar: single evaluation time
            %       Vector: [t1, t2, ...] for multi-time evaluation
            %
            % Outputs:
            %   Oe - Evaluated operator result

            if isempty(O)
                Oe = [];
            elseif isa(O, 'function_handle')
                switch nargin(O)
                    case 1
                        Oe = O(U);
                    case 2
                        if isscalar(t)
                            te = t;
                        else
                            te = t(1);
                        end
                        Oe = O(U, te);
                    case 3
                        core.except.assert(length(t) >= 2, 'InvalidInput', ...
                            'Operator function requires at least 2 time values.');
                        Oe = O(U, t(1), t(2));
                    otherwise
                        core.except.assert(0, 'InvalidInput', ...
                            'Operator function requires 1, 2 or 3 arguments.');
                end
            else
                Oe = O;
            end
        end

        function tf = shouldUpdateRhsMass(obj, i)
            % SHOULDUPDATERHSMASS Check if RHS mass cache should be
            % updated.
            %
            %   tf = shouldUpdateRhsMass(obj, i) determines whether the RHS
            %   mass operator cache needs updating based on step size
            %   changes and cache state.
            %
            % Inputs:
            %   obj - The OdeIntegrator object
            %   i - Index for cache lookup (stage or step index)
            %
            % Outputs:
            %   tf - True if cache should be updated (logical)

            if isempty(i)
                tf = isempty(obj.MRhs) || obj.timeline.hasStepSizeChanged;
            else
                tf = isempty(obj.MRhs{i}) || obj.timeline.hasStepSizeChanged;
            end
        end

        function obj = setRhsMass(obj, i, M)
            % SETRHSMASS Update the RHS mass operator cache.
            %
            %   obj = setRhsMass(obj, i, M) stores the evaluated RHS mass
            %   operator in the appropriate cache location.
            %
            % Inputs:
            %   obj - The OdeIntegrator object
            %   i - Index for cache storage (stage or step index)
            %   M - Evaluated RHS mass operator to cache
            %
            % Outputs:
            %   obj - The OdeIntegrator object

            if isempty(i)
                obj.MRhs = M;
            else
                obj.MRhs{i} = M;
            end
        end

        function M = getRhsMass(obj, i)
            % GETRHSMASS Retrieve mass operator from cache.
            %
            %   M = getRhsMass(obj, i) retrieves the cached RHS mass
            %   operator value for the specified index.
            %
            % Inputs:
            %   obj - The OdeIntegrator object
            %   i - Index for cache retrieval (stage or step index)
            %
            % Outputs:
            %   M - Cached RHS mass operator value

            if isempty(i)
                M = obj.MRhs;
            else
                M = obj.MRhs{i};
            end
        end
    end

    methods (Abstract)
        % SETLHS Set left-hand side matrix/matrices for linear systems.
        obj = setLhs(obj, varargin)

        % SETRHS Set right-hand side vector/vectors for linear systems.
        obj = setRhs(obj, varargin)

        % STAGE Advance one stage of the integration method.
        U = stage(obj, varargin)

        % STEP Advance one complete time step.
        U = step(obj, varargin)
    end
end