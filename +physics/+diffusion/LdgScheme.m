classdef LdgScheme < physics.diffusion.DiffusionScheme
    % LDGSCHEME Local discontinuous Galerkin scheme for diffusion equations.
    %
    %   LdgScheme implements a local discontinuous Galerkin method for solving
    %   diffusion equations using discontinuous Galerkin discretization.
    %   The scheme supports both periodic and Dirichlet boundary conditions
    %   with elliptic assembly operations.
    %
    % See also:
    %   physics.diffusion.DiffusionScheme 

    properties (Constant)
        Name = 'LDG' % Scheme name
    end

    properties
        Assembly % Diffusion assembly
    end

    methods
        function obj = LdgScheme(options)
            % LDGSCHEME Constructor for LdgScheme.
            %
            %   obj = LdgScheme() creates a local discontinuous Galerkin
            %   scheme with default configuration.
            %
            %   obj = LdgScheme(file=file) creates a local discontinuous
            %   Galerkin scheme and loads configuration from the specified
            %   @a filename.
            %
            %   obj = LdgScheme(config=config) creates a local
            %   discontinuous Galerkin scheme with the specified @a config
            %   configuration.
            %
            %   obj = LdgScheme(file=file, config=config) loads
            %   configuration from file and then overrides with struct
            %   values.
            %
            arguments
                options.file {mustBeTextScalar} = ''
                options.config struct = struct()
            end
            
            obj@physics.diffusion.DiffusionScheme(file=options.file, config=options.config);
            obj.Assembly = [];
        end

        function state = initialize(obj, state)
            % INITIALIZE Initialize the local scheme.
            %
            %   state = initialize(obj, state) sets up the LDG scheme by
            %   assembling the local diffusion operator, fitting initial
            %   conditions, and preparing visualization components.
            %
            arguments
                obj physics.diffusion.LdgScheme
                state physics.state.SpatialState
            end

            %< Fit initial condition
            state.setDof('U', state.XDisc.fit(obj.Config.ic));

            %< Create diffusion assembly
            obj.Assembly = approx.assembly.EllipticAssembly( ...
                state.XDisc, obj.Config.bcType, obj.Config.xPenaltyType);

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
            %
            arguments
                obj physics.diffusion.LdgScheme
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
            %   solution history and applies the diffusion operator with
            %   appropriate boundary conditions.
            %
            arguments
                obj physics.diffusion.LdgScheme
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
            %
            arguments
                obj physics.diffusion.LdgScheme
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
            %
            arguments
                obj physics.diffusion.LdgScheme
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
            % SETUP Setup default configuration options for LDG schemes.
            
            %< Call parent setup
            obj = setup@physics.diffusion.DiffusionScheme(obj);
            
            %< Add LDG-specific required configuration options
            obj.addConfig('xBasisOrder', default=1, ...
                validator=@(x) isnumeric(x) && isscalar(x) && (x > 0));
            obj.addConfig('xBasisType', default='nodal', ...
                validator=@(x) (ischar(x) || isstring(x)) && ismember(lower(x), {'modal', 'nodal'}));
            obj.addConfig('xBasisPattern', default='Q', ...
                validator=@(x) (ischar(x) || isstring(x)) && ismember(upper(x), {'P', 'Q'}));
            obj.addConfig('xPenaltyType', default={'', 'right'; 'left', ''}, ...
                validator=@(x) iscell(x) && isequal(size(x), [2, 2]));
        end

        function obj = configure(obj)
            % CONFIGURE Generate LDG-specific configuration dependencies.
            
            %< Call parent configure
            obj = configure@physics.diffusion.DiffusionScheme(obj);

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
    end

    methods (Access = private) %< Rule 246: helper methods
        function obj = setLinearTerm(obj, state) %#ok<INUSD>
            % SETLINEARTERM Assemble linear operator sparse matrix.

            field = obj.Config.diffusion;
            obj.L = obj.Assembly.assembleMatrix(field);
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
            % COMPUTESOURCETERM Compute source term for boundary conditions.

            func = obj.Config.bc;
            field = obj.Config.diffusion;
            S = obj.Assembly.assembleVector(field, func, args={'t',t});
        end
    end
end
