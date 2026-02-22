classdef Scheme < physics.scheme.ConfigParser
    % SCHEME Base class for all numerical schemes.
    %
    %   Scheme provides the fundamental framework for numerical schemes
    %   with built-in configuration management capabilities. The class
    %   manages the complete simulation lifecycle through abstract methods
    %   that must be implemented by subclasses: initialize, preStep, step,
    %   postStep, and finalize.
    %
    %   Configuration can be loaded from files or provided as structures.
    %   All schemes inherit powerful configuration management including
    %   validation, defaults, and file parsing capabilities.
    %
    % See also:
    %   physics.scheme.ConfigParser, physics.scheme.OdeScheme,
    %   physics.scheme.MolScheme
    
    methods
        function obj = Scheme(options)
            % SCHEME Construct an instance of Scheme.
            %
            %   obj = Scheme() creates a scheme with default configuration.
            %
            %   obj = Scheme(file=filename) creates a scheme and loads
            %   configuration from the specified @a filename.
            %
            %   obj = Scheme(config=config) creates a scheme with the
            %   specified @a configStruct configuration.
            %
            %   obj = Scheme(file=file, config=config) loads configuration
            %   from file and then overrides with struct values.
            
            arguments
                options.file {mustBeTextScalar} = ''
                options.config struct = struct()
            end
            
            %< Call ConfigParser constructor
            obj@physics.scheme.ConfigParser();
            
            %< Setup default configuration options
            obj.setup();
            
            %< Load from file if specified
            obj.load(file=options.file, config=options.config);
        end

        function obj = load(obj, options)
            % LOAD Load configuration from file and/or struct.
            %
            %   obj = load(obj, file=file, config=config) loads configuration from
            %   the specified @a file and/or struct @a config, applying
            %   validation and defaults as necessary.

            arguments
                obj physics.scheme.Scheme
                options.file {mustBeTextScalar} = ''
                options.config struct = struct()
            end

            %< Load from file if specified
            if ~isempty(options.file)
                obj.loadFromFile(options.file);
            end
            
            %< Load from struct (overrides file values)
            if ~isempty(fieldnames(options.config))
                obj.loadFromStruct(options.config);
            end
            
            %< Configure scheme-specific settings
            obj.configure();
        end
        
        function obj = loadFromStruct(obj, config)
            % LOADFROMSTRUCT Load configuration from struct.
            %
            %   obj = loadFromStruct(obj, config) loads configuration
            %   values from the struct @a config, applying validation
            %   automatically through setConfig method.
            
            arguments
                obj physics.scheme.Scheme
                config struct
            end
            
            fieldNames = fieldnames(config);
            for i = 1:length(fieldNames)
                name = fieldNames{i};
                value = config.(name);
                obj.setConfig(name, value);
            end
        end

        function state = initialize(obj, state)
            % INITIALIZE Initialize the scheme.
            
            arguments
                obj physics.scheme.Scheme
                state physics.state.State
            end
            
            core.except.assert(0, 'NotImplemented', ...
                'Method initialize is not implemented by class %s.', class(obj));
        end
        
        function state = preStep(obj, state)
            % PRESTEP Preprocess operations before each step.
            
            arguments
                obj physics.scheme.Scheme
                state physics.state.State
            end
            
            core.except.assert(0, 'NotImplemented', ...
                'Method preStep is not implemented by class %s.', class(obj));
        end
        
        function state = step(obj, state)
            % STEP Perform one step.
            
            arguments
                obj physics.scheme.Scheme
                state physics.state.State
            end
            
            core.except.assert(0, 'NotImplemented', ...
                'Method step is not implemented by class %s.', class(obj));
        end
        
        function state = postStep(obj, state)
            % POSTSTEP Postprocess operations after each step.
            
            arguments
                obj physics.scheme.Scheme
                state physics.state.State
            end
            
            core.except.assert(0, 'NotImplemented', ...
                'Method postStep is not implemented by class %s.', class(obj));
        end
        
        function state = finalize(obj, state)
            % FINALIZE Finalize scheme and solution state.
            
            arguments
                obj physics.scheme.Scheme
                state physics.state.State
            end
            
            core.except.assert(0, 'NotImplemented', ...
                'Method finalize is not implemented by class %s.', class(obj));
        end
    end
    
    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options.
            
            obj.addConfig('verbose', default=0, validator=@(x) x >= 0);
        end

        function obj = configure(obj)
            % CONFIGURE Generate configuration dependencies and validate
            % setup.
            %
            %   obj = configure(obj) processes the loaded configuration
            %   to set up derived properties, validate requirements, and
            %   prepare the scheme for use. Override in subclasses to add
            %   specific configuration logic.
            
            % Base class has no additional configuration requirements
            % Subclasses should override this method for specific setup
        end
    end
end