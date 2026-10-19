classdef Scheme < handle
    % SCHEME Base class for all numerical schemes.
    %
    %   Scheme provides the fundamental framework for numerical schemes.
    %   The class manages the complete simulation lifecycle through
    %   abstract methods that must be implemented by subclasses:
    %   initialize, preStep, step, postStep, and finalize.

    properties (SetAccess = private, GetAccess = public)
        config % Configuration structure
    end

    properties (Abstract, Constant)
        REQUIRED_CONFIGS % Required configurations
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
            [tf, mf] = core.validate.hasFields(config, obj.REQUIRED_CONFIGS);
            core.except.assert(tf, 'InvalidInput', ...
                'Required configuration %s is missing.', mf);
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