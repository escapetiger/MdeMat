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
    %   \f$F(u)\f$ is a nonlinear operator, and \f$S(t)\f$ is a source
    %   term. Any of these terms can be empty. In particular, M = [] means
    %   the mass operator reduces to an identity operator; {L, F, S} = []
    %   means they vanish from the equation.
    %
    %   The contribution of mass operator M is determined by the physical
    %   simulation framework:
    %   * In the Eulerian framework, M is space-only so that it can be
    %     represented by a unique static matrix.
    %   * In the Semi-Lagrangian framework, M is time-dependent. The Lhs
    %     mass term is forced to be a unique static matrix as in Eulerian
    %     framework. However, the Rhs mass term should be updated as time
    %     evolves. In particular, M is a matrix under periodic boundary
    %     conditions, while M is a pair formed by a matrix and a vector
    %     under Dirichlet boundary conditions.
    %
    % See also:
    %   approx.odeint.BeIntegrator, approx.odeint.BdfIntegrator,
    %   approx.odeint.DirkIntegrator, approx.odeint.ExrkIntegrator
    
    properties
        NSteps % Number of time steps required for the method
        NStages % Number of stages required for the method
        Timeline % Timeline object for time step management
        U0 % Cell array of solution history
        MLhs % Left-hand side mass term
        MRhs % Right-hand side mass term(s)
        Solver % Stage-specific linear solvers
    end
    
    methods
        function obj = OdeIntegrator(nSteps, nStages, final)
            % ODEINTEGRATOR Constructor for OdeIntegrator.
            %
            %   obj = OdeIntegrator(nSteps, nStages, final) creates an ODE
            %   integrator with the specified parameters. This constructor
            %   is called by concrete subclasses to initialize the base
            %   class functionality.
            
            arguments
                nSteps {mustBeInteger, mustBePositive}
                nStages {mustBeInteger, mustBePositive}
                final {mustBeNumeric, mustBeNonnegative}
            end
            
            obj.NSteps = nSteps;
            obj.NStages = nStages;
            obj.Timeline = approx.mesh.DynamicTimeline(final, nSteps=nSteps);
            obj.reset();
        end
        
        function obj = reset(obj)
            % RESET Reset the integrator to initial state.
            %
            %   obj = reset(obj) clears solution history, resets the
            %   timeline, and initializes storage arrays for left-hand side
            %   matrices and right-hand side vectors. This method is called 
            %   automatically during construction and can be called to restart
            %   integration.
            
            arguments
                obj approx.odeint.OdeIntegrator
            end
            
            obj.U0 = cell(obj.NStages, obj.NSteps);
            obj.Timeline.reset();
            
            %< Initialize mass term storage
            if obj.NStages > 1 && obj.NSteps == 1
                obj.MLhs = [];
                obj.MRhs = cell(obj.NStages+1, 1);
            elseif obj.NStages == 1 && obj.NSteps > 1
                obj.MLhs = [];
                obj.MRhs = cell(1, obj.NSteps);
            else
                obj.MLhs = [];
                obj.MRhs = [];
            end
        end
        
        function obj = update(obj, U)
            % UPDATE Update the solution history with a new solution.
            %
            %   obj = update(obj, U) shifts the existing solution history
            %   and adds the new solution at the most recent position. For
            %   multi-step methods, this maintains the required number of
            %   previous solutions. For single-step methods, simply stores
            %   the current solution.
            
            arguments
                obj approx.odeint.OdeIntegrator
                U {mustBeNumeric}
            end
            
            %< Shift history for multi-step methods
            if obj.NSteps > 1
                obj.U0 = circshift(obj.U0, [0, 1]);
            end
            
            %< Store new solution at all stage positions
            obj.add(U);
        end
        
        function obj = add(obj, U)
            % ADD Add a new solution at the newest position of history.
            %
            %   obj = add(obj, U) adds the new solution at the most recent
            %   position.
            
            arguments
                obj approx.odeint.OdeIntegrator
                U {mustBeNumeric}
            end
            
            for i = 1:obj.NStages
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
            
            arguments
                obj approx.odeint.OdeIntegrator
            end
            
            obj.Timeline.advance();
        end
        
        function obj = setTimeStep(obj, options)
            % SETTIMESTEP Set time step for the integrator timeline.
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
                obj approx.odeint.OdeIntegrator
                options.dt {mustBeNonnegative} = []
                options.h {mustBePositive} = []
                options.C {mustBePositive} = 1
                options.p {mustBePositive} = 1
            end
            
            if isa(obj.Timeline, 'approx.mesh.DynamicTimeline')
                if ~isempty(options.dt)
                    obj.Timeline.setTimeStep(dt=options.dt);
                else
                    obj.Timeline.setTimeStep(h=options.h, C=options.C, p=options.p);
                end
            else
                core.except.assert(0, 'InvalidInput', ...
                    'Must specify either dt or h for time step.');
            end
        end
    end
    
    methods (Access = protected)
        function Me = evalMLhs(obj, M, t)
            % EVALMLhS Evaluate mass operator for left-hand side.
            %
            %   Me = evalMLhs(obj, M, t) evaluates the mass operator
            %   for use in the left-hand side matrix construction. Handles
            %   empty, constant, and time-dependent mass operators.
            
            core.except.assert(isscalar(t), 'InvalidInput', ...
                'Lhs mass requires a scalar time.');
            
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
        
        function Me = evalMRhs(obj, i, M, t)
            % EVALMRhS Evaluate mass operator for right-hand side.
            %
            %   Me = evalMRhs(obj, i, M, t) evaluates the mass
            %   operator for use in right-hand side vector construction.
            %   Handles caching for efficiency.
            
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
                        if obj.shouldUpdateMRhs(i)
                            obj.setMRhs(i, M(te));
                        end
                        Me = obj.getMRhs(i);
                    case 2
                        core.except.assert(length(t) >= 2, 'InvalidInput', ...
                            'Mass function requires at least 2 time values.');
                        Me = cell(1, 2);
                        if obj.shouldUpdateMRhs(i)
                            [Me{1:2}] = M(t(1), t(2));
                            obj.setMRhs(i, Me{2});
                        else
                            Me{2} = obj.getMRhs(i);
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
        
        function Oe = evalOp(~, O, U, t)
            % EVALOP Evaluate a general operator (L, F, or S).
            %
            %   Oe = evalOp(obj, O, U, t) evaluates linear,
            %   nonlinear, or source operators with flexible argument
            %   handling.
            
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
        
        function tf = shouldUpdateMRhs(obj, i)
            % SHOULDUPDATEMRhS Check if Rhs mass cache should be
            % updated.
            %
            %   tf = shouldUpdateMRhs(obj, i) determines whether the Rhs
            %   mass operator cache needs updating based on step size
            %   changes and cache state.
            
            if isempty(i)
                tf = isempty(obj.MRhs) || obj.Timeline.HasStepSizeChanged;
            else
                tf = isempty(obj.MRhs{i}) || obj.Timeline.HasStepSizeChanged;
            end
        end
        
        function obj = setMRhs(obj, i, M)
            % SETMRhS Update the Rhs mass operator cache.
            %
            %   obj = setMRhs(obj, i, M) stores the evaluated Rhs mass
            %   operator in the appropriate cache location.
            
            if isempty(i)
                obj.MRhs = M;
            else
                obj.MRhs{i} = M;
            end
        end
        
        function M = getMRhs(obj, i)
            % GETMRhS Retrieve mass operator from cache.
            %
            %   M = getMRhs(obj, i) retrieves the cached Rhs mass
            %   operator value for the specified index.
            
            if isempty(i)
                M = obj.MRhs;
            else
                M = obj.MRhs{i};
            end
        end
    end
    
    methods (Abstract)
        % STEP Advance one complete time step.
        U = step(obj, options)
    end
end