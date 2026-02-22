classdef MolScheme < physics.scheme.OdeScheme
    % MOLSCHEME Method of Lines simulation controller.
    %
    %   MolScheme provides a high-level interface for executing Method of
    %   Lines simulations. It inherits from OdeScheme and implements the
    %   run method that encapsulates the common simulation pattern using
    %   the template method design pattern.
    %
    %   Subclasses must implement the abstract methods: initialize,
    %   preStep, step, postStep, and finalize to define the specific
    %   physics and numerical methods.
    %
    %   All simulation parameters and profiling tools are configured
    %   through the config structure.
    %
    % See also:
    %   physics.scheme.OdeScheme

    properties (Dependent)
        ShouldCheckpointing % Check if current step requires checkpointing
    end
    
    methods
        function state = run(obj, state)
            % RUN Execute Method of Lines simulation.
            %
            %   state = run(obj, state) executes a single simulation run
            %   using the template method pattern. The method orchestrates
            %   the simulation workflow by calling the abstract methods
            %   that must be implemented by subclasses.
            
            arguments
                obj physics.scheme.MolScheme
                state physics.state.SpatialState
            end
            
            timer = obj.Config.timer;
            visualizer = obj.Config.visualizer;
            
            %< Initialize timer
            timer.reset();
            
            %< Initialize visualizer
            visualizer.initialize();

            %< Reset ODE integrator
            obj.TDisc.reset();
            
            %< Initialize scheme and state
            fprintf('[M] Initialization.\n');
            timer.start(record='Initialization');
            state = obj.initialize(state);
            timer.stop(record='Initialization');

            %< Initial visualization
            fprintf('[M] Visualization.\n');
            timer.start(record='Visualization');
            obj.plot(state);
            timer.stop(record='Visualization');
            
            %< Main time stepping loop
            while ~obj.TDisc.Timeline.IsExhausted
                %< Preprocess
                state = obj.preStep(state);
                
                %< Integration step
                fprintf('[M] Integration.\n');
                timer.start(record='Integration');
                state = obj.step(state);
                timer.stop(record='Integration');
                record = obj.Config.timer.find('Integration');
                fprintf('%s.\n', obj.echo(record.Duration));

                %< Visualization
                fprintf('[M] Visualization.\n');
                timer.start(record='Visualization');
                obj.plot(state);
                timer.stop(record='Visualization');

                %< Checkpointing
                if obj.ShouldCheckpointing
                    obj.saveCkptState(state);
                end
                
                %< Postprocess
                state = obj.postStep(state);
            end
            
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
                obj physics.scheme.MolScheme
                state physics.state.SpatialState
            end
            
            timer = obj.Config.timer;
            analyzer = obj.Config.analyzer;
            
            analyzer.initialize();
            
            if ~analyzer.IsEnabled
                table = [];
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
                    analyzer.addAbsolute(state, args={'t', obj.Config.tFinal});
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
            table = analyzer.analyze(obj.Config.verbose);
        end

        function tf = get.ShouldCheckpointing(obj)
            % GET.SHOULDCHECKPOINTING Check if current step requires
            % checkpointing.

            if isempty(obj.Config.ckptTimeStamps)
                tf = false;
                return;
            end

            tol = 1e-8;
            tf = any(abs(obj.TDisc.Timeline.Next - obj.Config.ckptTimeStamps) < tol);
        end
    end
    
    methods (Access = protected)
        function obj = setup(obj)
            % SETUP Setup default configuration options for MoL schemes.
            
            %< Call parent setup
            setup@physics.scheme.OdeScheme(obj);
            
            %< Add MoL-specific configuration options
            obj.addConfig('timer', default=core.chrono.Timer(), ...
                validator=@(x) isa(x, 'core.chrono.Timer'));
            
            obj.addConfig('visualizer', default=physics.visual.Visualizer(), ...
                validator=@(x) isa(x, 'physics.visual.Visualizer'));
            
            obj.addConfig('analyzer', default=physics.analysis.Analyzer(), ...
                validator=@(x) isa(x, 'physics.analysis.Analyzer'));

            obj.addConfig('dt', default=[], ...
                validator=@(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));

            obj.addConfig('cfl', default=[], ...
                validator=@(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));

            obj.addConfig('ckptTimeStamps', default=[], ...
                validator=@(x) isempty(x) || (isvector(x) && all(x >= 0)));
            
            obj.addConfig('ckptDir', default='ckpts', ...
                validator=@(x) ischar(x) || isstring(x));
            
            obj.addConfig('ckptPrefix', default='ckpt', ...
                validator=@(x) ischar(x) || isstring(x));

            obj.addConfig('ckptPostfix', default='', ...
                validator=@(x) ischar(x) || isstring(x));

            obj.addConfig('ckptStateName', default='state', ...
                validator=@(x) ischar(x) || isstring(x));

            obj.addConfig('ckptFigureName', default='figure', ...
                validator=@(x) ischar(x) || isstring(x));

            obj.addConfig('ckptFigureFormat', default='pdf', ...
                validator=@(x) ischar(x) || isstring(x));
        end

        function obj = configure(obj)
            % CONFIGURE Configure MoL scheme with current configuration
            % struct.
            
            %< Call parent configure
            obj = configure@physics.scheme.OdeScheme(obj);

            core.except.assert(~isempty(obj.Config.cfl) || ~isempty(obj.Config.dt), ...
                'InvalidInput', ...
                'Either dt or cfl must be specified in configuration.');
        end

        function obj = plot(obj, state)
            % PLOT Generate visualization for current solution state.
            %
            %   obj = plot(obj, state) creates plots for the current
            %   solution state using default styling based on dataset
            %   types.
            
            arguments
                obj physics.scheme.MolScheme
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
            
            [V0, n0] = obj.precompute(space, dofs);
            
            while visualizer.Timeline.isBefore(obj.TDisc.Timeline.Next)
                visualizer.update(space, dofs, ...
                    timeline=obj.TDisc.Timeline, V0 = V0, n0 = n0);
                visualizer.render();
                visualizer.Timeline.advance();
            end

            %< Checkpointing
            if obj.ShouldCheckpointing
                obj.saveCkptFigure(gcf);
            end
        end

        function obj = saveCkptState(obj, state)
            % SAVECKPTSTATE Save current state as checkpoint.
            
            fprintf('[M] Saving checkpoint state.\n');
            ckptStateName = obj.Config.ckptStateName;
            ckptDir = obj.Config.ckptDir;
            ckptPrefix = obj.Config.ckptPrefix;
            next = obj.TDisc.Timeline.Next;
            next = sprintf('t%.1f', next);
            next = strrep(next, '.', 'p');
            if isempty(obj.Config.ckptPostfix)
                ckptPostfix = next;
            else
                ckptPostfix = strjoin({obj.Config.ckptPostfix, next}, '_');
            end
            ckptName = sprintf('%s_%s_%s.mat', ckptPrefix, ckptStateName, ckptPostfix);
            state.save(fullfile(ckptDir, ckptName));
        end

        function obj = saveCkptFigure(obj, fig)
            % SAVECKPTFIGURE Save current figure as checkpoint.
            
            fprintf('[M] Saving checkpoint figure.\n');
            ckptFigureName = obj.Config.ckptFigureName;
            ckptFigureFormat = obj.Config.ckptFigureFormat;
            ckptDir = obj.Config.ckptDir;
            ckptPrefix = obj.Config.ckptPrefix;
            next = obj.TDisc.Timeline.Next;
            next = sprintf('t%.1f', next);
            next = strrep(next, '.', 'p');
            if isempty(obj.Config.ckptPostfix)
                ckptPostfix = next;
            else
                ckptPostfix = strjoin({obj.Config.ckptPostfix, next}, '_');
            end
            ckptName = sprintf('%s_%s_%s.%s', ...
                ckptPrefix, ckptFigureName, ckptPostfix, ckptFigureFormat);
            exportgraphics(fig, fullfile(ckptDir, ckptName), ... 
                'ContentType', 'vector');
        end
        
        function msg = echo(obj, runtime)
            % ECHO Generate status message for current time step.
            
            msg = sprintf('[%d] now = %.6g; step size = %.6g; runtime = %.2f s', ...
                obj.TDisc.Timeline.Count, ...
                obj.TDisc.Timeline.Next, ...
                obj.TDisc.Timeline.StepSize, ...
                runtime);
        end
        
        function [V0, n0] = precompute(obj, space, dofs)
            % TIMEINTERPOLATE Compute temporal interpolation data.
            
            arguments
                obj physics.scheme.MolScheme
                space approx.space.FiniteElementSpace
                dofs struct
            end
            
            t0 = obj.TDisc.Timeline.Next;
            t1 = obj.TDisc.Timeline.Next;
            
            if t1 <= t0
                V0 = [];
                n0 = 0;
                return;
            end
            
            n0 = min(obj.TDisc.Timeline.NSteps, obj.TDisc.Timeline.Count);
            V0 = cell(1, n0);
            for i = 1:n0
                W0 = reshape(obj.TDisc.U0{i}, space.NGlobalDofs, []);
                dofNames = fieldnames(dofs);
                offset = 0;
                for j = 1:length(dofNames)
                    dofName = dofNames{j};
                    if isempty(obj.Config.visualizer.Components.(dofName))
                        continue;
                    end
                    nComponentsPerDof = size(dofs.(dofName), 2);
                    V0{i}.(dofName) = W0(:, offset+(1:nComponentsPerDof));
                    offset = offset + nComponentsPerDof;
                end
            end
        end
    end
end