classdef UwdgScheme < physics.advection.AdvectionScheme
    % UWDGSCHEME Upwind discontinuous Galerkin scheme.
    %
    %   UwdgScheme implements an upwind discontinuous Galerkin scheme for
    %   solving advection equations. The method uses upwind numerical fluxes
    %   to handle the hyperbolic nature of advection problems, providing
    %   stability and sharp resolution of solution features.
    %
    %   The upwind flux selection provides natural stabilization for
    %   hyperbolic problems but may require smaller time steps compared
    %   to semi-Lagrangian methods for stability.
    %
    % See also:
    %   physics.scheme.MolScheme, physics.advection.SldgScheme

    properties (Constant)
        Name = 'UWDG' % Scheme name
    end

    properties
        Assembly % Upwind assembly
    end

    methods
        function obj = UwdgScheme(options)
            % UWDGSCHEME Construct an instance of UwdgScheme.
            %
            %   obj = UwdgScheme() creates an upwind DG scheme with default
            %   configuration.
            %
            %   obj = UwdgScheme(file=file) creates an upwind DG scheme and
            %   loads configuration from the specified @a filename.
            %
            %   obj = UwdgScheme(config=config) creates an upwind DG scheme
            %   with the specified @a config configuration.
            %
            %   obj = UwdgScheme(file=file, config=config) loads
            %   configuration from file and then overrides with struct
            %   values.
            
            arguments
                options.file {mustBeTextScalar} = ''
                options.config struct = struct()
            end
            
            obj@physics.advection.AdvectionScheme(file=options.file, config=options.config);
            obj.Assembly = [];
        end

        function state = initialize(obj, state)
            % INITIALIZE Initialize the upwind scheme.
            %
            %   state = initialize(obj, state) sets up the upwind scheme by
            %   assembling the advection operator, fitting initial
            %   conditions, and preparing visualization components.
            
            arguments
                obj physics.advection.UwdgScheme
                state physics.state.SpatialState
            end

            %< Fit initial condition
            state.setDof('U', state.XDisc.fit(obj.Config.ic));

            %< Create upwind operator
            obj.Assembly = approx.assembly.UpwindAssembly(state.XDisc, obj.Config.bcType);

            %< Set linear and source terms
            obj.setLinearTerm(state);
            obj.setSourceTerm(state);
        end

        function state = preStep(obj, state)
            % PRESTEP Preprocess operations before time step.
            %
            %   state = preStep(obj, state) performs preprocessing
            %   including time step size computation based on CFL condition
            %   and coefficient updates for multi-step methods.
            
            arguments
                obj physics.advection.UwdgScheme
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
            % STEP Perform one time step.
            %
            %   state = step(obj, state) advances the solution by one time
            %   step using the upwind scheme. The method updates the
            %   solution history and applies the advection operator with
            %   appropriate boundary conditions.
            
            arguments
                obj physics.advection.UwdgScheme
                state physics.state.SpatialState
            end

            %< Update ODE integrator with current solution
            obj.TDisc.update(state.Dofs.U);

            %< Perform ODE integration step
            state.Dofs.U = obj.TDisc.step(L=obj.L, S=obj.S);
        end

        function state = postStep(obj, state)
            % POSTSTEP Postprocess operations after time step.
            %
            %   state = postStep(obj, state) performs postprocessing
            %   including visualization updates and timeline advancement
            %   for the next integration step.
            
            arguments
                obj physics.advection.UwdgScheme
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
                obj physics.advection.UwdgScheme
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
            % SETUP Setup default configuration options for UWDG schemes.
            
            %< Call parent setup
            obj = setup@physics.advection.AdvectionScheme(obj);
            
            %< Add UWDG-specific required configuration options
            obj.addConfig('xBasisOrder', default=1, ...
                validator=@(x) isnumeric(x) && isscalar(x) && (x > 0));
            obj.addConfig('xBasisType', default='nodal', ...
                validator=@(x) (ischar(x) || isstring(x)) && ismember(lower(x), {'modal', 'nodal'}));
            obj.addConfig('xBasisPattern', default='Q', ...
                validator=@(x) (ischar(x) || isstring(x)) && ismember(upper(x), {'P', 'Q'}));
        end

        function obj = configure(obj)
            % CONFIGURE Generate UWDG-specific configuration
            % dependencies.
            
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

        function obj = setLinearTerm(obj, ~)
            % SETLINEARTERM Assemble linear operator sparse matrix.

            field = obj.Config.advection;
            obj.L = -obj.Assembly.assembleMatrix(field);
        end

        function obj = setSourceTerm(obj, ~)
            % SETSOURCETERM Set up source term for boundary conditions.

            switch obj.Config.bcType
                case 'periodic'
                    obj.S = [];
                case 'dirichlet'
                    obj.S = @(U, t) obj.computeSourceTerm(U, t);
            end
        end

        function S = computeSourceTerm(obj, ~, t)
            % COMPUTESOURCETERM Compute source term for boundary
            % conditions.

            func = obj.Config.bc;
            field = obj.Config.advection;
            direction = obj.Config.advection;
            S = -obj.Assembly.assembleVector(field, direction, func, args={'t', t});
        end
    end
end