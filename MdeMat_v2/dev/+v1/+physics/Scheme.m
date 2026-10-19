classdef Scheme < handle
    % SCHEME Base class for physics simulation schemes.
    %
    %   The Scheme class provides a common framework for physics simulation
    %   schemes, managing time integration loops, performance monitoring,
    %   visualization, and analysis workflows. It defines the essential
    %   structure for numerical solution procedures while allowing specific
    %   implementations through abstract method definitions.
    %
    %   Key features include:
    %   - Standardized time integration workflow
    %   - Integrated performance timing and monitoring
    %   - Visualization and analysis pipeline management
    %   - Configurable verbose output and reporting
    %   - Extensible architecture for different physics models
    %
    % Notes:
    %   This is an abstract base class that cannot be instantiated
    %   directly. Concrete subclasses must implement all abstract methods
    %   to define specific physics and numerical schemes.

    properties (Constant)
        RECORDS = {'Integration', 'Visualization', 'Analysis'}; % Timer record categories
    end

    properties
        config % Configuration structure containing simulation parameters
        tDisc % Temporal discretization object for time integration
        timer % Performance timer for monitoring execution phases
        visualizer % Visualization object for solution plotting
        analyzer % Analysis object for error computation and assessment
    end

    methods
        function obj = Scheme(config, tDisc, visualizer, analyzer)
            % SCHEME Constructor for Scheme.
            %
            %   obj = Scheme(config, tDisc, visualizer, analyzer) creates
            %   a base scheme object with the specified configuration,
            %   temporal discretization, visualization, and analysis
            %   components.
            %
            % Inputs:
            %   config - Configuration structure with simulation parameters
            %   tDisc - Temporal discretization object
            %   visualizer - Visualization object for solution plotting
            %   analyzer - Analysis object for error computation
            %
            % Outputs:
            %   obj - Constructed Scheme object

            obj.config = config;
            obj.tDisc = tDisc;
            obj.timer = core.chrono.Timer(obj.RECORDS, config.verbose);
            obj.visualizer = visualizer;
            obj.analyzer = analyzer;
        end

        function state = run(obj, state)
            % RUN Execute the complete time integration loop.
            %
            %   state = run(obj, state) performs the full time integration
            %   from initial to final time, managing preprocessing,
            %   integration steps, and postprocessing with performance
            %   monitoring throughout the simulation.
            %
            % Inputs:
            %   obj - The Scheme object
            %   state - Initial solution state
            %
            % Outputs:
            %   state - Final solution state after time integration

            %< Time integration loop
            while ~obj.tDisc.timeline.isFinished
                %< Preprocess step
                state = obj.preStep(state);

                %< Integration step
                msg = sprintf( ...
                    "[%d] now = %f; step size = %f; runtime = %.2f s;", ...
                    obj.tDisc.timeline.count, ...
                    obj.tDisc.timeline.now, ...
                    obj.tDisc.timeline.h(1), ...
                    obj.timer.records(1).duration);
                obj.timer.start(1, msg);
                state = obj.step(state);
                obj.timer.stop(1);

                %< Postprocess step
                state = obj.postStep(state);
            end
        end
    end

    methods (Abstract)
        % INITIALIZE Initialize the scheme and solution state.
        state = initialize(obj, state)

        % PRESTEP Preprocess operations before each time step.
        state = preStep(obj, state)

        % STEP Perform one time integration step.
        state = step(obj, state)

        % POSTSTEP Postprocess operations after each time step.
        state = postStep(obj, state)

        % FINALIZE Finalize simulation and perform final analysis.
        state = finalize(obj, varargin)
    end
end