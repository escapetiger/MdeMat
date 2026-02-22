classdef LdgScheme < physics.poisson.PoissonScheme
    % LDGSCHEME Base class for all LDG schemes solving Poisson equations.
    %
    %   LdgScheme implements an LDG scheme for solving Poisson 
    %   equations. 
    %
    % See also:
    %   physics.scheme.PoissonScheme

    properties (Constant)
        Name = 'LDG' % Scheme name
    end
    
    properties
        Laplace % Laplace assembly
        Source % Source assembly
    end
    
    methods
        function state = initialize(obj, state)
            % INITIALIZE Initialize the Poisson scheme.
            %
            %   state = initialize(obj, state) sets up the Poisson scheme
            %   by assembling the elliptic operator, setting up source
            %   terms, and preparing the linear system for solving.
            
            arguments
                obj physics.poisson.LdgScheme
                state physics.state.SpatialState
            end

            %< Laplace operator
            obj.Laplace = approx.assembly.AdjointAssembly( ...
                state.XDisc, obj.Config.bcType, obj.Config.xPenaltyType);
            
            %< Source operator
            obj.Source = approx.assembly.SourceAssembly(state.XDisc);
            
            %< Add degrees of freedom and coefficients
            state.setDof('U', []);
            state.setDof('Q', []);
            
            %< Set LHS and RHS terms
            obj.Config.timer.start(record='Assembly');
            obj.setLhsTerm(state);
            obj.Config.timer.stop(record='Assembly');
            obj.setRhsTerm(state);
        end
        
        function state = step(obj, state)
            % STEP Solve the Poisson equation.
            %
            %   state = step(obj, state) solves the linear system using
            %   the assembled matrices. The method applies appropriate
            %   boundary conditions and returns the updated state with
            %   the solution coefficients.
            
            arguments
                obj physics.poisson.LdgScheme
                state physics.state.SpatialState
            end
            
            switch obj.Config.bcType
                case 'periodic'
                    W = obj.Solver.solve(A=obj.A, b=obj.b);
                    W = reshape(W(1:end-1), [], state.XDisc.NDims+1);
                    state.Dofs.U = W(:, 1);
                    state.Dofs.Q = W(:, 2:end);
                case 'dirichlet'
                    W = obj.Solver.solve(A=obj.A, b=obj.b);
                    W = reshape(W, [], state.XDisc.NDims+1);
                    state.Dofs.U = W(:, 1);
                    state.Dofs.Q = W(:, 2);
            end
        end
        
        function state = finalize(obj, state)
            % FINALIZE Finalize simulation with scheme and state reporting.
            %
            %   state = finalize(obj, state) performs final reporting
            %   including scheme configuration details and state
            %   information summary.
            
            arguments
                obj physics.poisson.PoissonScheme
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
            % SETUP Setup default configuration options for Poisson schemes.
            
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
            % CONFIGURE Generate Poisson-specific configuration
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
            if obj.Config.hasAuxiliaryVariable
                obj.Config.components.Q = 1;
            end

            %< Set visualization options
            visualizer = obj.Config.visualizer;
            visualizer.setComponents(obj.Config.components);
            visualizer.addExact('U', obj.Config.exact.U);
            visualizer.addExact('Q', obj.Config.exact.Q);
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
            analyzer.addExact('Q', obj.Config.exact.Q);
        end
        
        function obj = setLhsTerm(obj, state)
            % SETLHSTERM Assemble left-hand side matrix.

            D = obj.Laplace.assembleMatrix(eye(state.nXDims));
            nd = state.XDisc.NDims;
            ng = state.XDisc.NGlobalDofs;
            obj.A = repmat({sparse(ng, ng)}, nd+1, nd+1);
            for d = 1:nd
                obj.A{1, 1} = obj.A{1, 1} + D{1, 1}{d};
                obj.A{1, d+1} = -D{1, 2}{d};
                obj.A{d+1, 1} = -D{2, 1}{d};
            end
            for d = 1:nd
                obj.A{d+1, d+1} = speye(ng, ng);
            end
            obj.A = core.linalg.block(obj.A);

            if strcmpi(obj.Config.bcType, 'periodic')
                c = obj.Source.assembleVector(@(x) ones(1, size(x, 2)));
                c = [c; repmat(sparse(ng, 1), nd, 1)];
                obj.A = [obj.A, c; c.', 0];
            end
        end
        
        function obj = setRhsTerm(obj, state)
            % SETRHSTERM Set up right-hand side vector.

            ng = state.XDisc.NGlobalDofs;
            obj.b = repmat({sparse(ng, 1)}, state.XDisc.NDims + 1, 1);
            f = obj.Source.assembleVector(obj.Config.force);
            obj.b{1} = f;
            if strcmpi(obj.Config.bcType, 'periodic') 
                obj.b = vertcat(obj.b{:});
                obj.b = [obj.b; 0];
                return; 
            end
            b0 = obj.Laplace.assembleVector(eye(state.nXDims), obj.Config.bc);
            for d = 1:state.XDisc.NDims
                obj.b{1} = obj.b{1} - b0{1, 1}{d};
                obj.b{d+1} = b0{2, 1}{d};
            end
            obj.b = vertcat(obj.b{:});
        end
    end
end