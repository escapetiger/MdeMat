classdef SteadyScheme < physics.Scheme
    % STEADYSCHEME Steady simulation controller.
    %
    %   SteadyScheme provides a high-level interface for executing steady
    %   simulations. It inherits from Scheme and implements the run
    %   method that encapsulates the common simulation pattern using the
    %   template method design pattern.
    %
    %   Subclasses must implement the abstract methods: initialize,
    %   preStep, step, postStep, and finalize to define the specific
    %   physics and numerical methods.
    %
    %   All simulation parameters and profiling tools are configured
    %   through the config structure.
    %
    % See Also:
    %   physics.Scheme

    methods
        function state = run(obj, state)
            % RUN Execute Method of Lines simulation.
            %
            %   state = run(obj, state) executes a single simulation run
            %   using the template method pattern. The method orchestrates
            %   the simulation workflow by calling the abstract methods
            %   that must be implemented by subclasses.
            %
            % Inputs:
            %   obj - The SteadyScheme object
            %   state - Initial solution state
            %
            % Outputs:
            %   state - Final solution state

            %< Input validation (Rule 235)
            core.except.assert(isa(state, 'physics.State'), ...
                'InvalidInput', 'state must be a physics.State object');

            %< Reset timer if available
            if isfield(obj.config, 'timer')
                obj.config.timer.reset();
            end

            %< Initialize scheme and state
            state = obj.initialize(state);

            %< Initialize visualizer
            if isfield(obj.config, 'visualizer')
                obj.config.visualizer.initialize();
            end

            %< Preprocess
            state = obj.preStep(state);

            %< Integration step
            if isfield(obj.config, 'timer')
                obj.config.timer.start(1, '[M] Integration.');
            end
            state = obj.step(state);
            if isfield(obj.config, 'timer')
                obj.config.timer.stop(1);
            end

            %< Visualization
            if isfield(obj.config, 'visualizer') && obj.config.visualizer.isEnabled
                if isfield(obj.config, 'timer')
                    obj.config.timer.start(2, '[M] Visualization.');
                end
                obj.config.visualizer.plot(state.xDisc, state.dofs);
                if isfield(obj.config, 'timer')
                    obj.config.timer.stop(2);
                end
            end

            %< Postprocess
            state = obj.postStep(state);

            %< Finalize simulation
            state = obj.finalize(state);
        end

        function state = converge(obj, state)
            % CONVERGE Execute convergence study with mesh refinement.
            %
            %   state = converge(obj, state) executes multiple simulation
            %   runs with progressively refined meshes to study convergence
            %   behavior. Uses the run method for each level and handles
            %   error analysis with proper timing.
            %
            % Inputs:
            %   obj - The SteadyScheme object
            %   state - Initial state with coarsest mesh
            %
            % Outputs:
            %   state - Final state from finest mesh level

            %< Input validation
            core.except.assert(isa(state, 'physics.SpatialState'), ...
                'InvalidInput', 'Convergence study requires SpatialState');
            core.except.assert(isfield(obj.config, 'analyzer'), ...
                'InvalidOperation', 'Convergence study requires analyzer in config');

            %< Initialize analyzer for convergence study
            analyzer = obj.config.analyzer;
            analyzer.initialize();
            nLevels = analyzer.nLevels;
            coarse = [];

            if ~analyzer.hasExact                
                nLevels = nLevels + 1;
            end

            %< Refinement loop
            for iLevel = 1:nLevels
                fprintf('[M] Refinement Level %d/%d \n', iLevel, nLevels);

                %< Run simulation
                state = obj.run(state);

                %< Error analysis with timing
                if isfield(obj.config, 'timer')
                    obj.config.timer.start(3, '[M] Error Evaluation.');
                end

                %< Add new mesh level
                analyzer.addLevel(state.xDisc.mesh);

                %< Error evaluation
                if analyzer.hasExact
                    analyzer.addAbsolute(state.xDisc, state.dofs);
                else
                    if iLevel > 1
                        analyzer.addRelative(state.xDisc, state.dofs, coarse.xDisc, coarse.dofs);
                    end
                    coarse = state;
                end

                if isfield(obj.config, 'timer')
                    obj.config.timer.stop(3);
                end

                %< Report timing information including error evaluation
                if isfield(obj.config, 'timer') && obj.config.verbose > 0
                    obj.config.timer.report();
                end

                %< Refine state for next iteration
                if iLevel < nLevels
                    state = state.refine(1);
                end

                fprintf('=====================================================\n');
            end

            %< Analyze convergence results
            analyzer.analyze(obj.config.verbose);
        end
    end

    methods (Access = protected)
        function obj = parse(obj)
            % PARSE Check and parse the contained configuration struct.
            %
            %   obj = parse(obj) checks and parses the contained
            %   configuration struct, calling parent parse method and
            %   adding SteadyScheme-specific validation.
            %
            % Inputs:
            %   obj - The SteadyScheme object
            %
            % Outputs:
            %   obj - The SteadyScheme object

            %< Call parent parse method
            parse@physics.Scheme(obj);

            %< Validate optional profiler objects
            if isfield(obj.config, 'timer')
                tf = isa(obj.config.timer, 'core.chrono.Timer');
                core.except.assert(tf, 'InvalidInput', ...
                    'timer must be a Timer object.');
            end

            if isfield(obj.config, 'visualizer')
                tf = isa(obj.config.visualizer, 'physics.Visualizer');
                core.except.assert(tf, 'InvalidInput', ...
                    'visualizer must be a Visualizer object.');
            end

            if isfield(obj.config, 'analyzer')
                tf = isa(obj.config.analyzer, 'profilers.Analyzer');
                core.except.assert(tf, 'InvalidInput', ...
                    'analyzer must be an Analyzer object.');
            end
        end
    end
end