classdef PdgScheme < physics.poisson.PoissonScheme
    % PDGSCHEME Base class for PDG schemes solving Poisson equations.
    %
    %   PdgScheme implements a PDG scheme for solving Poisson
    %   equations. The PDG scheme is a variant of the DG scheme
    %   that solves only the primal formulation of the Poisson equation.
    %   It can be derived from the LDG scheme by eliminating the
    %   dual variables and constraints.
    %
    % See also:
    %   physics.scheme.PoissonScheme

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
                obj physics.poisson.PdgScheme
                state physics.state.SpatialState
            end
            
            %< Laplace operator
            obj.Laplace = approx.assembly.EllipticAssembly( ...
                state.XDisc, obj.Config.bcType, obj.Config.xPenaltyType);
            
            %< Source operator
            obj.Source = approx.assembly.SourceAssembly(state.XDisc);
            
            %< Add degrees of freedom and coefficients
            state.setDof('U', []);
            
            %< Set LHS and RHS terms
            obj.Config.timer.start(record='Assembly');
            obj.setLhsTerm(state);
            obj.setRhsTerm(state);
            obj.Config.timer.stop(record='Assembly');
        end
        
        function state = step(obj, state)
            % STEP Solve the Poisson equation.
            %
            %   state = step(obj, state) solves the linear system using
            %   the assembled matrices. The method applies appropriate
            %   boundary conditions and returns the updated state with
            %   the solution coefficients.
            
            arguments
                obj physics.poisson.PdgScheme
                state physics.state.SpatialState
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
                obj physics.poisson.PdgScheme
                state physics.state.SpatialState
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
            obj = setup@physics.poisson.PoissonScheme(obj);
            
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
            % CONFIGURE Generate DG-specific configuration
            % dependencies.

            %< Call parent configure
            obj = configure@physics.poisson.PoissonScheme(obj);
            
            %< Determine spatial discretization identifier
            xId = sprintf('%s%d', upper(obj.Name), obj.Config.xBasisOrder);
            obj.setConfig('xId', xId);

            %< Set plot title
            eId = obj.Config.eId;
            obj.setConfig('titlePrefix', sprintf('%s with %s', eId, xId));

            %< Set components
            obj.Config.components.U = 1;

            %< Set visualization options
            visualizer = obj.Config.visualizer;
            visualizer.setComponents(obj.Config.components);
            visualizer.addExact('U', obj.Config.exact.U);
            visualizer.setTitlePrefix(obj.Config.titlePrefix);

            %< Set analysis options
            analyzer = obj.Config.analyzer;
            analyzer.setComponents(obj.Config.components);
            if isempty(fieldnames(analyzer.Reductions))
                compFields = fieldnames(obj.Config.components);
                red = struct();
                for i = 1:length(compFields)
                    red.(compFields{i}) = struct('x', {obj.Config.errorTypes}, 'v', '');
                end
                analyzer.setReductions(red);
            end
            analyzer.addExact('U', obj.Config.exact.U);
        end
        
        function obj = setLhsTerm(obj, state)
            % SETLHSTERM Assemble left-hand side matrix.

            obj.A = -obj.Laplace.assembleMatrix(eye(state.nXDims));
            if strcmpi(obj.Config.bcType, 'periodic')
                c = obj.Source.assembleVector(@(x) ones(1, size(x, 2)));
                obj.A = [obj.A, c; c.', 0];
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