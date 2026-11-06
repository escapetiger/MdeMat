classdef PdgScheme < mdemat.solver.elliptic.PoissonScheme
    % PDGSCHEME Base class for PDG schemes solving Poisson equations.
    %
    %   PdgScheme implements a PDG scheme for solving Poisson
    %   equations. The PDG scheme is a variant of the DG scheme
    %   that solves only the primal formulation of the Poisson equation.
    %   It can be derived from the LDG scheme by eliminating the
    %   dual variables and constraints.
    %
    % See also:
    %   mdemat.solver.elliptic.PoissonScheme

    properties (Constant)
        Name = 'PDG' % Scheme name
    end
    
    properties
        Laplace % Laplace assembly
        Source % Source assembly
    end
    
    methods
        function state = initialize(obj, state)
            % INITIALIZE Initialize the PDG scheme.
            %
            %   state = initialize(obj, state) sets up the PDG scheme
            %   by assembling the elliptic operator, setting up source
            %   terms, and preparing the linear system for solving.
            
            arguments
                obj mdemat.solver.elliptic.PdgScheme
                state mdemat.common.fem.SpatialState
            end
            
            %< Laplace operator
            obj.Laplace = mdemat.common.fem.EllipticAssembly( ...
                state.XDisc, obj.Config.bcType, obj.Config.xPenaltyType);
            
            %< Source operator
            obj.Source = mdemat.common.fem.SourceAssembly(state.XDisc);
            
            %< Add degrees of freedom and coefficients
            state.setDof('U', []);
            
            %< Set LHS and RHS terms
            timer = obj.getObserver('timer');
            if ~isempty(timer)
                [~, interval] = timer.start(record='Assembly');
                assemblyGuard = onCleanup(@() timer.cancel(interval));
            end
            obj.setLhsTerm(state);
            obj.setRhsTerm(state);
            if ~isempty(timer), timer.stop(interval=interval); end
        end
        
        function state = step(obj, state)
            % STEP Solve the Poisson equation.
            %
            %   state = step(obj, state) solves the linear system using
            %   the assembled matrices. The method applies appropriate
            %   boundary conditions and returns the updated state with
            %   the solution coefficients.
            
            arguments
                obj mdemat.solver.elliptic.PdgScheme
                state mdemat.common.fem.SpatialState
            end
            
            switch obj.Config.bcType
                case 'periodic'
                    U = obj.Solver.solve(A=obj.A, b=obj.b);
                    state.Dofs.U = U(1:end-1);
                case 'dirichlet'
                    state.Dofs.U = obj.Solver.solve(A=obj.A, b=obj.b);
            end
        end
        
        function state = finalize(obj, state)
            % FINALIZE Finalize simulation with scheme and state reporting.
            %
            %   state = finalize(obj, state) performs final reporting
            %   including scheme configuration details and state
            %   information summary.
            
            arguments
                obj mdemat.solver.elliptic.PdgScheme
                state mdemat.common.fem.SpatialState
            end
            
            %< Report state information
            state.report();
            
            %< Report boundary condition type
            fprintf('[R] Boundary Conditions: %s\n', obj.Config.bcType);
        end
    end
    
    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options for DG schemes.
            
            %< Call parent setup
            obj = setup@mdemat.solver.elliptic.PoissonScheme(obj);
            
            %< Add LDG-specific required configuration options
            obj.addConfig('xBasisOrder', default=1, ...
                validator=@(x) isnumeric(x) && isscalar(x) && (x > 0));
            obj.addConfig('xBasisType', default='nodal', ...
                validator=@(x) (ischar(x) || isstring(x)) && ismember(lower(x), {'modal', 'nodal'}));
            obj.addConfig('xBasisPattern', default='Q', ...
                validator=@(x) (ischar(x) || isstring(x)) && ismember(upper(x), {'P', 'Q'}));
            obj.addConfig('xPenaltyType', default={'boundary', 'right'; 'left', ''}, ...
                validator=@(x) iscell(x) && isequal(size(x), [2, 2]));
        end

        function obj = configure(obj)
            % CONFIGURE Generate DG-specific configuration
            % dependencies.

            %< Call parent configure
            obj = configure@mdemat.solver.elliptic.PoissonScheme(obj);
            
            %< Determine spatial discretization identifier
            xId = sprintf('%s%d', upper(obj.Name), obj.Config.xBasisOrder);
            obj.setDerived('xId', xId);

            %< Set plot title
            eId = obj.Config.eId;
            obj.setDerived('titlePrefix', sprintf('%s with %s', eId, xId));

            %< Set components
            components = obj.Config.components;
            components.U = 1;

            obj.setConfig('components', components);
        end

        function obj = setLhsTerm(obj, state)
            % SETLHSTERM Assemble left-hand side matrix.

            obj.A = -obj.Laplace.assembleMatrix(eye(state.nXDims));
            if strcmpi(obj.Config.bcType, 'periodic')
                c = obj.Source.assembleVector(@(x) ones(1, size(x, 2)));
                ell = state.XDisc.integrationWeights();
                obj.A = [obj.A, c; ell.', 0];
            end
        end
        
        function obj = setRhsTerm(obj, state)
            % SETRHSTERM Set up right-hand side vector.

            f = obj.Source.assembleVector(obj.Config.force);
            obj.b = f;
            if strcmpi(obj.Config.bcType, 'periodic')
                obj.b = [obj.b; 0];
                return; 
            end
            b0 = obj.Laplace.assembleVector(eye(state.nXDims), obj.Config.bc);
            obj.b = obj.b + b0;
        end
    end
end
