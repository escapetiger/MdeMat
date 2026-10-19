classdef Scheme < handle
    % SCHEME Base class for all numerical schemes.
    %
    %   Scheme provides the fundamental framework for numerical schemes.
    %   The class manages the complete simulation lifecycle through
    %   abstract methods that must be implemented by subclasses:
    %   initialize, preStep, step, postStep, and finalize.

    properties
        config % Configuration structure
    end

    methods
        function obj = Scheme(config)
            % SCHEME Constructor for Scheme.
            %
            %   obj = Scheme(config) creates a base scheme object with the
            %   specified configuration structure.
            %
            % Inputs:
            %   config - Configuration structure
            %
            % Outputs:
            %   obj - Constructed Scheme object

            obj.config = config;
            obj.parse();
        end
    end

    methods (Access = protected)
        function obj = parse(obj)
            % PARSE Check and parse the contained configuration struct.
            %
            %   obj = parse(obj) checks and parses the contained
            %   configuration struct.
            %
            % Inputs:
            %   obj - The Scheme object
            %
            % Outputs:
            %   obj - The Scheme object
        
            core.except.assert(isstruct(obj.config), 'InvalidInput', ...
                'config must be a struct.');
            
            %< Set default verbose level if not specified
            if ~isfield(obj.config, 'verbose')
                obj.config.verbose = 0;
            end
        end   
    end

    methods (Abstract)
        % INITIALIZE Initialize scheme and solution state.
        state = initialize(obj, state)

        % PRESTEP Preprocess operations before each step.
        state = preStep(obj, state)

        % STEP Perform one step.
        state = step(obj, state)

        % POSTSTEP Postprocess operations after each step.
        state = postStep(obj, state)

        % FINALIZE Finalize scheme and solution state.
        state = finalize(obj, state)
    end
end