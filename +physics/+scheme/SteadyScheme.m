classdef SteadyScheme < physics.scheme.Scheme
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
    % See also:
    %   physics.scheme.Scheme
    
    methods
        function state = run(obj, state)
            % RUN Execute Method of Lines simulation.
            %
            %   state = run(obj, state) executes a single simulation run
            %   using the template method pattern. The method orchestrates
            %   the simulation workflow by calling the abstract methods
            %   that must be implemented by subclasses.
            
            arguments
                obj physics.scheme.SteadyScheme
                state physics.state.SpatialState
            end

            timer = obj.Config.timer;
            visualizer = obj.Config.visualizer;
            
            %< Initialize timer
            timer.reset();
            
            %< Initialize scheme and state
            state = obj.initialize(state);
            
            %< Initialize visualizer
            visualizer.initialize();

            %< Integration step
            fprintf('[M] Integration.\n');
            timer.start(record='Integration');
            state = obj.step(state);
            timer.stop(record='Integration');

            %< Visualization
            fprintf('[M] Visualization.\n');
            timer.start(record='Visualization');
            obj.plot(state);
            timer.stop(record='Visualization');
            
            %< Finalize simulation
            state = obj.finalize(state);
        end
        
        function [state, table] = converge(obj, state)
            % CONVERGE Execute convergence study with mesh refinement.
            %
            %   [state, table] = converge(obj, state) executes multiple
            %   simulation runs with progressively refined meshes to study
            %   convergence behavior. Uses the run method for each level
            %   and handles error analysis with proper timing.
            
            arguments
                obj physics.scheme.SteadyScheme
                state physics.state.SpatialState
            end

            timer = obj.Config.timer;
            analyzer = obj.Config.analyzer;

            analyzer.initialize();

            if ~analyzer.IsEnabled
                return;
            end
            
            nLevels = analyzer.NLevels;
            coarse = [];

            if ~analyzer.HasExact
                nLevels = nLevels + 1;
            end
            
            %< Refinement loop
            for iLevel = 1:nLevels
                fprintf('[M] Refinement Level %d/%d \n', iLevel, nLevels);
                
                %< Run simulation
                state = obj.run(state);
                
                %< Start timing
                fprintf('[M] Error Evaluation.\n');
                timer.start(record='Analysis');

                %< Add new mesh level
                analyzer.addLevel(state.XDisc.Mesh);
                
                %< Error evaluation
                if analyzer.HasExact
                    analyzer.addAbsolute(state);
                else
                    if iLevel > 1
                        analyzer.addRelative(state, coarse);
                    end
                    coarse = state;
                end
                
                %< Stop timing
                timer.stop(record='Analysis');

                %< Report timing information
                if obj.Config.verbose > 0
                    timer.report();
                end
                
                %< Refine state for next iteration
                if iLevel < nLevels
                    state = state.refine(1);
                end
                
                fprintf('=====================================================\n');
            end
            
            %< Analyze convergence results
            table = obj.Config.analyzer.analyze(obj.Config.verbose);
        end
    end
    
    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options for steady schemes.
            
            %< Call parent setup
            setup@physics.scheme.Scheme(obj);
            
            %< Add steady-specific configuration options
            obj.addConfig('timer', default=core.chrono.Timer(), ...
                validator=@(x) isa(x, 'core.chrono.Timer'));
            
            obj.addConfig('visualizer', default=physics.visual.Visualizer(), ...
                validator=@(x) isa(x, 'physics.visual.Visualizer'));
            
            obj.addConfig('analyzer', default=physics.analysis.Analyzer(), ...
                validator=@(x) isa(x, 'physics.analysis.Analyzer'));
        end

        function obj = plot(obj, state)
            % PLOT Generate visualization for current solution state.
            %
            %   obj = plot(obj, state) creates plots for the current
            %   solution state using default styling based on dataset
            %   types.

            arguments
                obj physics.scheme.SteadyScheme
                state physics.state.SpatialState
            end

            visualizer = obj.Config.visualizer;

            if ~visualizer.IsEnabled
                return;
            end

            space = state.XDisc;
            dofs = state.Dofs;

            visualizer.setRefXNodes(space.Element.Geometry);
            visualizer.setPhyXNodes(space.Mesh);

            visualizer.update(space, dofs);
            visualizer.render();
        end
    end
end