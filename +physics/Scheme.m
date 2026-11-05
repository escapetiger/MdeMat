classdef Scheme < handle
    % SCHEME Base class for all numerical schemes.
    %
    %   Scheme provides the fundamental framework for numerical simulation
    %   schemes including time integration, visualization, analysis, and
    %   performance monitoring. It defines the common interface and
    %   workflow for all concrete scheme implementations.
    %
    %   The class manages the complete simulation lifecycle through
    %   abstract methods that must be implemented by subclasses:
    %   initialize, preStep, step, postStep, and finalize. It coordinates
    %   temporal discretization, performance timing, visualization, and
    %   error analysis components.
    %
    % Notes:
    %   Subclasses must implement all abstract methods. The run method
    %   provides the main time integration loop with automatic performance
    %   monitoring and step management.
    %
    % See also:
    %   approx.odeint.OdeIntegrator, core.chrono.Timer,
    %   profilers.analysis.Analyzer, profilers.visual.Visualizer

    properties (Constant)
        RECORDS = {'Integration', 'Visualization', 'Analysis'}; % Timer record categories
    end

    properties
        config     % Configuration structure containing simulation parameters
        tDisc      % Temporal discretization object for time integration
        timer      % Performance timer for monitoring execution phases
        visualizer % Visualization object for solution plotting and output
        analyzer   % Analysis object for error computation and assessment
    end

    properties
        isFullyImplicit % Check fully implicit temporal discretization
        isBdf % Check BDF temporal discretization
        isDirk % Check DIRK temporal discretization
    end

    methods
        function obj = Scheme(config)
            % SCHEME Constructor for Scheme.
            %
            %   obj = Scheme(config) creates a base scheme object with the
            %   specified configuration structure. Initializes all
            %   component properties to empty, requiring explicit setup via
            %   setter methods.
            %
            % Inputs:
            %   config - Configuration structure
            %
            % Outputs:
            %   obj - Constructed Scheme object

            obj.config = config;
            obj.tDisc = [];
            obj.timer = [];
            obj.visualizer = [];
            obj.analyzer = [];
        end

        function tf = get.isFullyImplicit(obj)
            % GET.ISFULLYIMPLICIT True if fully implicit temporal
            % discretization is applied.

            tf1 = isa(obj.tDisc, 'approx.odeint.BeIntegrator');
            tf2 = isa(obj.tDisc, 'approx.odeint.BdfIntegrator');
            tf3 = isa(obj.tDisc, 'approx.odeint.DirkIntegrator');
            tf = tf1 | tf2 | tf3;
        end

        function tf = get.isBdf(obj)
            % GET.ISBDF True if BDF temporal discretization is applied.

            tf1 = isa(obj.tDisc, 'approx.odeint.BeIntegrator');
            tf2 = isa(obj.tDisc, 'approx.odeint.BdfIntegrator');
            tf = tf1 | tf2;
        end

        function tf = get.isDirk(obj)
            % GET.ISDIRK True if DIRK temporal discretization is applied.

            tf1 = isa(obj.tDisc, 'approx.odeint.BeIntegrator');
            tf2 = isa(obj.tDisc, 'approx.odeint.DirkIntegrator');
            tf = tf1 | tf2;
        end

        function obj = setTimer(obj, verbose)
            % SETTIMER Set performance timer with specified verbosity.
            %
            %   obj = setTimer(obj, verbose) creates and configures a
            %   performance timer for monitoring simulation phases
            %   including integration, visualization, and analysis
            %   operations.
            %
            % Inputs:
            %   obj - The Scheme object
            %   verbose - Verbosity flag for timer output (logical)
            %
            % Outputs:
            %   obj - The Scheme object with configured timer

            obj.timer = core.chrono.Timer(obj.RECORDS, verbose);
        end

        function obj = setTimeDiscretization(obj, tDiscId, tFinal)
            % SETTIMEDISCRETIZATION Set temporal discretization method.
            %
            %   obj = setTimeDiscretization(obj, tDiscId, tFinal) creates
            %   and configures the appropriate time integration method
            %   based on the specified identifier. Supports various
            %   implicit methods including backward Euler, BDF, and SDIRK
            %   schemes.
            %
            % Inputs:
            %   obj - The Scheme object
            %   tDiscId - Time discretization identifiers (string or cell)
            %   tFinal - Final integration time (positive scalar)
            %
            % Outputs:
            %   obj - The Scheme object with configured time discretization

            if ~iscell(tDiscId)
                splitting = [];
                tDiscId = {tDiscId};
            else
                splitting = tDiscId{1};
                tDiscId = tDiscId(2:end);
            end

            integrators = cell(size(tDiscId));
            for i = 1:length(tDiscId)
                 switch lower(tDiscId{i})
                    case 'be'
                        cls = 'BeIntegrator';
                    case 'bdf2'
                        cls = 'Bdf2Integrator';
                    case 'bdf3'
                        cls = 'Bdf3Integrator';
                    case 'sdirk2'
                        cls = 'Sdirk2Integrator';
                    case 'sdirk3'
                        cls = 'Sdirk3Integrator';
                    case 'fe'
                        cls = 'FeIntegrator';
                    case 'heun'
                        cls = 'HeunIntegrator';
                    case 'ssprk2'
                        cls = 'Ssprk2Integrator';
                    case 'ssprk3'
                        cls = 'Ssprk3Integrator';
                 end
                 integrators{i} = approx.odeint.(cls)(tFinal);
            end

            if ~isempty(splitting)
                switch lower(splitting)
                    case 'lie_trotter'
                        cls = 'LieTrotterIntegrator';
                    case 'strang'
                        cls = 'StrangIntegrator';
                end
                obj.tDisc = approx.odeint.(cls)(integrators);
            else
                obj.tDisc = integrators{1};
            end
        end

        function obj = setAnalyzer(obj, varargin)
            % SETANALYZER Set analysis component for error assessment.
            %
            %   obj = setAnalyzer(obj, varargin) creates and configures an
            %   analyzer object for computing solution errors, convergence
            %   rates, and other quantitative assessments during
            %   simulation.
            %
            % Inputs:
            %   obj - The Scheme object
            %   varargin - Name-value parameter pairs
            %
            % Outputs:
            %   obj - The Scheme object

            obj.analyzer = profilers.analysis.Analyzer(varargin{:});
        end

        function obj = setVisualizer(obj, varargin)
            % SETVISUALIZER Set visualization component for solution output.
            %
            %   obj = setVisualizer(obj, varargin) creates and configures
            %   a visualizer object for plotting solutions, generating
            %   animations, and producing graphical output during simulation.
            %
            % Inputs:
            %   obj - The Scheme object
            %   varargin - Name-value parameter pairs
            %
            % Outputs:
            %   obj - The Scheme object

            obj.visualizer = profilers.visual.Visualizer(varargin{:});
        end

        function obj = addDataset(obj, varargin)
            obj.visualizer.addDataset(varargin{:});
        end

        function obj = addStyle(obj, varargin)
            obj.visualizer.addStyle(varargin{:});
        end

        function state = run(obj, state)
            % RUN Execute the complete time integration loop.
            %
            %   state = run(obj, state) performs the full time integration
            %   from initial to final time, managing preprocessing,
            %   integration steps, and postprocessing with performance
            %   monitoring throughout the simulation. The method
            %   coordinates the abstract step methods defined by
            %   subclasses.
            %
            % Inputs:
            %   obj - The Scheme object
            %   state - Initial solution state structure
            %
            % Outputs:
            %   state - Final solution state after time integration

            %< Time integration loop
            while ~obj.tDisc.timeline.isFinished
                %< Preprocess step
                state = obj.preStep(state);

                %< Integration step with performance monitoring
                msg = sprintf( ...
                    "[%d] now = %f; step size = %f; runtime = %.2f s;", ...
                    obj.tDisc.timeline.count, ...
                    obj.tDisc.timeline.now, ...
                    obj.tDisc.timeline.dt, ...
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