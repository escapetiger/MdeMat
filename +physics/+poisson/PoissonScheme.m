classdef PoissonScheme < physics.scheme.SteadyScheme
    % POISSONSCHEME Base class for all schemes solving poisson equations.
    %
    %   PoissonScheme implements a numerical scheme for solving poisson
    %   equations. This base class defines the basic validation and parse
    %   logic of configuration struct.
    %
    % See also:
    %   physics.scheme.SteadyScheme

    properties
        Solver % Linear solver
        A % LHS matrix
        b % RHS vector
    end
    
    methods
        function obj = PoissonScheme(options)
            % POISSONSCHEME Construct an instance of PoissonScheme.
            %
            %   obj = PoissonScheme() creates a PoissonScheme with default
            %   configuration.
            %
            %   obj = PoissonScheme(file=file) creates a PoissonScheme and
            %   loads configuration from the specified @a filename.
            %
            %   obj = PoissonScheme(config=config) creates a PoissonScheme
            %   with the specified @a config configuration.
            %
            %   obj = PoissonScheme(file=file, config=config) loads
            %   configuration from file and then overrides with struct
            %   values.
            
            arguments
                options.file {mustBeTextScalar} = ''
                options.config struct = struct()
            end

            obj@physics.scheme.SteadyScheme(file=options.file, config=options.config);

            obj.Solver = core.linalg.LinearSolver(nnzTh=1e6, condTh=Inf);
            obj.A = [];
            obj.b = [];
        end
    end

    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options for Poisson schemes.
            
            %< Call parent setup
            obj = setup@physics.scheme.SteadyScheme(obj);
            
            %< Add Poisson-specific required configuration options
            obj.addConfig('bc', default=[], ...
                validator=@(x) true);
            obj.addConfig('force', default=@(x) 0, ...
                validator=@(x) isa(x, 'function_handle'));
            obj.addConfig('eId', default='Poisson', ...
                validator=@(x) (ischar(x) || isstring(x)));
            obj.addConfig('errorTypes', default={'L1', 'L2', 'max'}, ...
                validator=@(x) iscell(x) && all(cellfun(@(y) ischar(y) || isstring(y), x)));
            obj.addConfig('exact', default=[], ...
                validator=@(x) isstruct(x) || isempty(x));
            obj.addConfig('components', default=[], ...
                validator=@(x) isstruct(x) && all(structfun(@(y) isempty(y) || isnumeric(y) && isscalar(y) && (y > 0), x)));
        end

        function obj = configure(obj)
            % CONFIGURE Generate Poisson-specific configuration
            % dependencies.

            %< Call parent configure
            obj = configure@physics.scheme.SteadyScheme(obj);
            
            %< Determine boundary condition type based on bc value
            if isempty(obj.Config.bc)
                obj.setConfig('bcType', 'periodic');
            else
                obj.setConfig('bcType', 'dirichlet');
            end
        end
    end
end