classdef SldgScheme < physics.advection.AdvectionScheme
    % SLDGSCHEME Semi-Lagrangian discontinuous Galerkin scheme.
    %
    %   SldgScheme implements a Semi-Lagrangian discontinuous Galerkin (SLDG)
    %   method for solving advection equations. The method combines the stability
    %   advantages of Lagrangian methods with the flexibility of discontinuous
    %   Galerkin spatial discretizations.
    %
    %   The SLDG method traces characteristics backward in time and uses
    %   high-order interpolation to evaluate solutions at departure points,
    %   making it particularly effective for transport-dominated problems
    %   with large time steps and minimal numerical diffusion.
    %
    %   The semi-Lagrangian approach allows for larger time steps compared
    %   to traditional Eulerian methods, making it computationally
    %   efficient for advection-dominated flows.
    %
    % See Also:
    %   physics.OdeScheme, physics.advection.UwdgScheme

    properties (Constant)
        Name = 'SLDG' % Scheme name
    end

    properties
        Assembly % Semi-Lagrangian assembly
    end

    methods
        function obj = SldgScheme(options)
            % SEMILAGRANGIANSCHEME Construct an instance of SemiLagrangianScheme.
            %
            %   obj = SemiLagrangianScheme() creates a semi-Lagrangian scheme
            %   with default configuration.
            %
            %   obj = SemiLagrangianScheme(file=file) creates a semi-Lagrangian
            %   scheme and loads configuration from the specified @a filename.
            %
            %   obj = SemiLagrangianScheme(config=config) creates a
            %   semi-Lagrangian scheme with the specified @a config
            %   configuration.
            
            arguments
                options.file {mustBeTextScalar} = ''
                options.config struct = struct()
            end
            
            obj@physics.advection.AdvectionScheme(file=options.file, config=options.config);
            obj.Assembly = [];
        end

        function state = initialize(obj, state)
            % INITIALIZE Initialize the SL scheme.
            %
            %   state = initialize(obj, state) sets up the semi-Lagrangian
            %   discontinuous Galerkin scheme by configuring finite element
            %   data, assembling the advection operator, fitting initial
            %   conditions, and preparing visualization components.
            
            arguments
                obj physics.advection.SldgScheme
                state physics.state.SpatialState
            end

            %< Fit initial condition
            state.setDof('U', state.XDisc.fit(obj.Config.ic));

            %< Create Semi-Lagrangian assembly
            obj.Assembly = approx.assembly.SemiLagrangianAssembly( ...
                state.XDisc, obj.Config.bcType);

            %< Set mass terms
            obj.setMassTerm(state);
        end

        function state = preStep(obj, state)
            % PRESTEP Preprocess operations before time step.
            %
            %   state = preStep(obj, state) performs preprocessing
            %   including time step size computation based on CFL condition
            %   and coefficient updates for multi-step methods.
            
            arguments
                obj physics.advection.SldgScheme
                state physics.state.SpatialState
            end

            if obj.hasConfig('cfl') && ~isempty(obj.Config.cfl)
                h = state.XDisc.Mesh.computeMeasure();
                obj.TDisc.setTimeStep(h = h, C = obj.Config.cfl, p = 1);
            end

            if obj.hasConfig('dt') && ~isempty(obj.Config.dt)
                obj.TDisc.setTimeStep(dt = obj.Config.dt);
            end
        end


        function state = step(obj, state)
            % STEP Perform one semi-Lagrangian time step.
            %
            %   state = step(obj, state) advances the solution by one time
            %   step using the semi-Lagrangian transport operator. The
            %   method updates the solution history and applies the
            %   characteristic-based transport with appropriate boundary
            %   condition.
            
            arguments
                obj physics.advection.SldgScheme
                state physics.state.SpatialState
            end

            %< Update ODE integrator with current solution
            obj.TDisc.update(state.Dofs.U);

            %< Perform ODE integration step
            state.Dofs.U = obj.TDisc.step(M=obj.M);
        end

        function state = postStep(obj, state)
            % POSTSTEP Postprocess operations after time step.
            %
            %   state = postStep(obj, state) performs postprocessing
            %   including visualization updates and timeline advancement
            %   for the next integration step.
            
            arguments
                obj physics.advection.SldgScheme
                state physics.state.SpatialState
            end

            %< Advance timeline for next step
            obj.TDisc.advance();
        end

        function state = finalize(obj, state)
            % FINALIZE Finalize simulation with scheme and state reporting.
            %
            %   state = finalize(obj, state) performs final reporting
            %   including scheme configuration details and state
            %   information summary.
            
            arguments
                obj physics.advection.SldgScheme
                state physics.state.SpatialState
            end

            %< Report state information
            state.report();

            fprintf('[R] ODE Integrator: %s\n', class(obj.TDisc));
            fprintf('[R] Final Time: %.6g\n', obj.Config.tFinal);
            if obj.hasConfig('cfl') && ~isempty(obj.Config.cfl)
                fprintf('[R] CFL number: %.3f\n', obj.Config.cfl);
            end
            if obj.hasConfig('dt') && ~isempty(obj.Config.dt)
                fprintf('[R] Time step size: %.3f\n', obj.Config.dt);
            end

            %< Report boundary condition type
            fprintf('[R] Boundary Conditions: %s\n', obj.Config.bcType);
        end
    end

    methods (Access = protected) 
        function obj = setup(obj)
            % SETUP Setup default configuration options for SLDG schemes.

            %< Call parent setup
            obj = setup@physics.advection.AdvectionScheme(obj);

            %< Add SLDG-specific required configuration options
            obj.addConfig('xBasisOrder', default=1, ...
                validator=@(x) isnumeric(x) && isscalar(x) && (x > 0));
            obj.addConfig('xBasisType', default='nodal', ...
                validator=@(x) (ischar(x) || isstring(x)) && ismember(lower(x), {'modal', 'nodal'}));
            obj.addConfig('xBasisPattern', default='Q', ...
                validator=@(x) (ischar(x) || isstring(x)) && ismember(upper(x), {'P', 'Q'}));
        end

        function obj = configure(obj)
            % CONFIGURE Generate SLDG-specific configuration dependencies.

            %< Call parent configure
            obj = configure@physics.advection.AdvectionScheme(obj);

            %< Determine spatial discretization identifier
            xId = sprintf('%s%d', upper(obj.Name), obj.Config.xBasisOrder);
            obj.setConfig('xId', xId);

            %< Set plot title
            eId = obj.Config.eId;
            obj.setConfig('titlePrefix', sprintf('%s with %s', eId, xId));

            %< Set visualization options
            visualizer = obj.Config.visualizer;
            visualizer.setTitlePrefix(obj.Config.titlePrefix);
        end

        function obj = setMassTerm(obj, ~)
            % SETMASSTERM Set up mass term.

            obj.M = @(tBegin, tEnd) obj.computeM(tBegin, tEnd);
        end

        function [V, M] = computeM(obj, tBegin, tEnd)
            % COMPUTEM Compute mass term.

            field = obj.Config.advection;
            
            obj.Assembly.reset();

            if nargout >= 2
                M = obj.Assembly.assembleMatrix(field, tBegin, tEnd);
            end

            if nargout >= 1
                V = [];
                if strcmpi(obj.Config.bcType, 'dirichlet')
                    if abs(tBegin-tEnd) < 1e-8, return; end
                    V = obj.Assembly.assembleVector(field, tBegin, tEnd, obj.Config.bc);
                end
            end
        end
    end
end