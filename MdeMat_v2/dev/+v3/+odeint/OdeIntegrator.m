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
                obj.lhs = cell(obj.nStages + 1, obj.nSteps);
                obj.rhs = cell(obj.nStages + 1, obj.nSteps);
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
                obj.MRhs = cell(obj.nStages + 1, 1);
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