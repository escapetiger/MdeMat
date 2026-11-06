classdef MolScheme < mdemat.internal.execution.OdeScheme
    %MOLSCHEME Execute the numerical lifecycle with optional value observations.
    %   Numerical state/history belongs to the scheme and ODE integrator.
    %   Tools receive only sampled values; plotting is optional and never
    %   supplies a required numerical input. Checkpoints use accepted steps.
    properties (Dependent)
        ShouldCheckpointing
    end
    properties (Access = private)
        PreviousPlot = struct()
        PreviousPlotTime = []
        PlotCursor = 1
        FigureFrame = 0
    end
    methods
        function state = run(obj, state, options)
            arguments
                obj mdemat.internal.execution.MolScheme
                state mdemat.common.fem.SpatialState
                options.resume (1,1) logical = false
            end
            obj.prepareObservers(false);
            timer = obj.getObserver('timer');
            timer.reset();
            obj.PreviousPlot = struct(); obj.PreviousPlotTime = [];
            obj.PlotCursor = 1; obj.FigureFrame = 0;
            if options.resume
                callbacks = coefficientCallbacks(state.Coefs);
                accepted = state.fromRecord(state.snapshot(omitCallbacks=true), callbacks=callbacks);
                workspace = mdemat.common.ode.restoreIntegrator( ...
                    mdemat.common.ode.snapshotIntegrator(obj.TDisc));
                target = state;
                resumeGuard = onCleanup(@() restoreAccepted(obj,target,accepted,workspace));
            else
                obj.TDisc.reset();
            end
            state = timed(timer,'Initialization',@() obj.initialize(state));
            if options.resume
                clear resumeGuard
                state = target;
            end
            timed(timer,'Visualization',@() obj.plot(state));
            if obj.checkpointAt(obj.TDisc.Timeline.Now), obj.saveCkptState(state); end
            while ~obj.TDisc.Timeline.IsExhausted
                state = obj.preStep(state);
                state = timed(timer,'Integration',@() obj.step(state));
                timed(timer,'Visualization',@() obj.plot(state));
                checkpoint = obj.ShouldCheckpointing;
                state = obj.postStep(state);
                if checkpoint, obj.saveCkptState(state); end
            end
            state = obj.finalize(state);
        end

        function saveCheckpoint(obj, filename, state, options)
            %SAVECHECKPOINT Persist an accepted state and complete ODE history.
            arguments
                obj mdemat.internal.execution.MolScheme
                filename {mustBeTextScalar}
                state mdemat.common.fem.SpatialState
                options.overwrite (1,1) logical = false
                options.createDirectory (1,1) logical = false
                options.omitCallbacks (1,1) logical = false
            end
            record = struct('format','mdemat.molCheckpoint','schemaVersion',2, ...
                'phase','completed-step','scheme',class(obj), ...
                'state',state.snapshot(omitCallbacks=options.omitCallbacks), ...
                'integrator',mdemat.common.ode.snapshotIntegrator(obj.TDisc), ...
                'configuration',mdemat.internal.execution.checkpointConfiguration(obj.Configuration.toStruct()));
            mdemat.tools.io.writeRecord(filename,record,overwrite=options.overwrite, ...
                createDirectory=options.createDirectory);
        end

        function state = loadCkptState(obj, filename, prototype, options)
            %LOADCKPTSTATE Validate everything before committing restored history.
            arguments
                obj mdemat.internal.execution.MolScheme
                filename {mustBeTextScalar}
                prototype mdemat.common.fem.SpatialState
                options.callbacks (1,1) struct = struct()
            end
            record = mdemat.internal.execution.readCheckpoint(filename);
            saved = mdemat.internal.execution.normalizeCheckpointConfiguration(record.configuration);
            expected = mdemat.internal.execution.checkpointConfiguration(obj.Configuration.toStruct());
            if ~strcmp(record.scheme,class(obj)) || ~isequaln(saved,expected)
                error('mdemat:persistence:InvalidCheckpoint', ...
                    'The checkpoint numerical model or explicit callback slots differ.');
            end
            state = prototype.fromRecord(record.state,callbacks=options.callbacks);
            integrator = mdemat.common.ode.restoreIntegrator(record.integrator);
            if ~strcmp(class(integrator),class(obj.TDisc)) || ...
                    integrator.Timeline.Final ~= obj.TDisc.Timeline.Final
                error('mdemat:persistence:InvalidCheckpoint','The ODE method/horizon differs.');
            end
            obj.TDisc = integrator;
        end

        function [state,tables] = converge(obj,state)
            [state,tables] = mdemat.internal.execution.runConvergence(obj,state);
        end

        function fields = samplePlotFields(obj,state,plan)
            fields = mdemat.common.fem.plotFieldSamples( ...
                state.XDisc,state.Dofs,obj.Observation.Plot.Components,plan);
        end

        function fields = samplePlotReferences(obj,state,plan,time)
            fields = mdemat.internal.execution.plotReferenceFields(state.XDisc, ...
                obj.Observation.Plot.Exacts,obj.Observation.Plot.Components,plan,time);
        end

        function tf = get.ShouldCheckpointing(obj)
            tf = obj.checkpointAt(obj.TDisc.Timeline.Next);
        end
    end
    methods (Access = protected)
        function obj = setup(obj)
            setup@mdemat.internal.execution.OdeScheme(obj);
            mdemat.internal.execution.defineStudyInputs(obj,'spatial');
            obj.addConfig('dt',default=[],validator=@(x) isempty(x) || ...
                (isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0));
            obj.addConfig('cfl',default=[],validator=@(x) isempty(x) || ...
                (isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0));
            obj.addConfig('ckptTimeStamps',default=[],validator=@(x) isempty(x) || ...
                (isnumeric(x)&&isreal(x)&&isvector(x)&&all(isfinite(x)&x>=0)));
            obj.addConfig('ckptDir',default='ckpts',validator=@isText);
            obj.addConfig('ckptOverwrite',default=false,validator=@(x) islogical(x)&&isscalar(x));
            obj.addConfig('ckptPrefix',default='ckpt',validator=@isText);
            obj.addConfig('ckptPostfix',default='',validator=@isText);
            obj.addConfig('ckptStateName',default='state',validator=@isText);
            obj.addConfig('ckptFigureName',default='figure',validator=@isText);
            obj.addConfig('ckptFigureFormat',default='pdf',validator=@isText);
        end

        function obj = configure(obj)
            configure@mdemat.internal.execution.OdeScheme(obj);
            if isempty(obj.Config.cfl) && isempty(obj.Config.dt)
                error('mdemat:internal:execution:MolScheme:InvalidInput', ...
                    'Specify a positive finite dt or cfl.');
            end
        end

        function obj = plot(obj,state)
            options = obj.Observation.Plot;
            if isempty(options.Modes) || isempty(options.PointsPerDimension) || ...
                    any(options.PointsPerDimension<=0), return; end
            currentTime = obj.TDisc.Timeline.Next;
            current = mdemat.internal.execution.sampleObservation(obj,state,currentTime);
            mdemat.internal.execution.configurePlotLayouts(obj,current);
            session = obj.getObserver('visualizer');
            if isempty(session) || ~session.IsEnabled, return; end
            if isempty(obj.PreviousPlotTime)
                previous = current; previousTime = currentTime;
            else
                previous = obj.PreviousPlot; previousTime = obj.PreviousPlotTime;
            end
            times = options.Times;
            if isempty(times)
                obj.renderFrame(session,current,currentTime,options.Title);
            else
                if ~isnumeric(times)||~isvector(times)||any(~isfinite(times))||any(diff(times)<=0)
                    error('mdemat:execution:MolScheme:ObservationTimes', ...
                        'Observation times must be a strictly increasing finite vector.');
                end
                tolerance = 32*eps(max(1,abs(currentTime)));
                while obj.PlotCursor<=numel(times) && times(obj.PlotCursor)<=currentTime+tolerance
                    time = times(obj.PlotCursor); obj.PlotCursor = obj.PlotCursor+1;
                    if time<previousTime-tolerance, continue; end
                    if currentTime>previousTime
                        fraction = min(1,max(0,(time-previousTime)/(currentTime-previousTime)));
                        datasets = mdemat.internal.execution.interpolateObservation(previous,current,fraction);
                    else
                        datasets = current;
                    end
                    % Accepted numeric frames are already sampled; only the
                    % time-dependent reference needs evaluation here.
                    exact = mdemat.internal.execution.sampleObservation(obj,state,time,false);
                    keys = fieldnames(exact);
                    for key=1:numel(keys)
                        if strcmp(exact.(keys{key}).Type,'exact')
                            datasets.(keys{key})=exact.(keys{key});
                        end
                    end
                    obj.renderFrame(session,datasets,time,options.Title);
                end
            end
            obj.PreviousPlot=current; obj.PreviousPlotTime=currentTime;
            if obj.ShouldCheckpointing
                obj.renderFrame(session,current,currentTime,options.Title);
                handles=session.FigureHandles; keys=fieldnames(handles);
                for key=1:numel(keys)
                    obj.saveCkptFigure(handles.(keys{key}),keys{key});
                end
            end
        end

        function obj = saveCkptState(obj,state)
            step = obj.TDisc.Timeline.Count-1;
            suffix = sprintf('step%09d',step);
            if ~isempty(obj.Config.ckptPostfix), suffix=[char(obj.Config.ckptPostfix),'_',suffix]; end
            name=sprintf('%s_%s_%s.mat',obj.Config.ckptPrefix,obj.Config.ckptStateName,suffix);
            obj.saveCheckpoint(fullfile(obj.Config.ckptDir,name),state, ...
                overwrite=obj.Config.ckptOverwrite,createDirectory=true,omitCallbacks=true);
        end

        function obj = saveCkptFigure(obj,fig,layoutName)
            if nargin<3,layoutName='profile';end
            suffix=sprintf('frame%09d_%s',obj.FigureFrame,layoutName);
            if ~isempty(obj.Config.ckptPostfix),suffix=[char(obj.Config.ckptPostfix),'_',suffix];end
            name=sprintf('%s_%s_%s.%s',obj.Config.ckptPrefix,obj.Config.ckptFigureName, ...
                suffix,obj.Config.ckptFigureFormat);
            destination=fullfile(obj.Config.ckptDir,name);
            if isfile(destination)&&~obj.Config.ckptOverwrite
                error('mdemat:io:FileExists','Refusing to overwrite %s.',destination);
            end
            if ~isfolder(obj.Config.ckptDir)
                [ok,message]=mkdir(obj.Config.ckptDir);
                if ~ok,error('mdemat:io:WriteFailed','%s',message);end
            end
            exportgraphics(fig,destination,'ContentType','vector');
        end

        function [V0, n0] = precompute(obj, space, dofs, options)
            % TIMEINTERPOLATE Compute temporal interpolation data.
            
            arguments
                obj mdemat.internal.execution.MolScheme
                space mdemat.common.fem.FiniteElementSpace
                dofs struct
                options.reconstruct = []
            end
            
            t0 = obj.TDisc.Timeline.Now;
            t1 = obj.TDisc.Timeline.Next;
            
            if t1 <= t0
                V0 = [];
                n0 = 0;
                return;
            end
            
            n0 = min([obj.TDisc.Timeline.NSteps, ...
                obj.TDisc.Timeline.Count, size(obj.TDisc.U0, 2)]);
            V0 = cell(1, n0);
            for i = 1:n0
                % Rows are RK stages, columns are accepted-state history.
                if ~isempty(options.reconstruct)
                    V0{i} = options.reconstruct(obj.TDisc.U0{1, i});
                    continue;
                end
                W0 = reshape(obj.TDisc.U0{1, i}, space.NGlobalDofs, []);
                V0{i} = struct();
                dofNames = fieldnames(dofs);
                offset = 0;
                for j = 1:length(dofNames)
                    dofName = dofNames{j};
                    nComponentsPerDof = size(dofs.(dofName), 2);
                    columns = offset+(1:nComponentsPerDof);
                    offset = offset + nComponentsPerDof;
                    if ~isfield(obj.Observation.Plot.Components, dofName) || ...
                            isempty(obj.Observation.Plot.Components.(dofName))
                        continue;
                    end
                    V0{i}.(dofName) = W0(:, columns);
                end
            end
        end
    end
    methods (Access = private)
        function tf = checkpointAt(obj,time)
            stamps=obj.Config.ckptTimeStamps;
            tf=~isempty(stamps)&&any(abs(time-stamps)<32*eps(max(1,abs(time))));
        end
        function renderFrame(obj,session,datasets,time,titleText)
            session.render(datasets,time=time,title=titleText);
            obj.FigureFrame=obj.FigureFrame+1;
        end
    end
end

function callbacks = coefficientCallbacks(coefficients)
callbacks = struct();
names = fieldnames(coefficients);
for index = 1:numel(names)
    if isa(coefficients.(names{index}), 'function_handle')
        callbacks.(names{index}) = coefficients.(names{index});
    end
end
end

function restoreAccepted(scheme, target, accepted, workspace)
target.Dofs = accepted.Dofs;
target.Coefs = accepted.Coefs;
target.History = accepted.History;
target.XDisc = accepted.XDisc;
if isprop(accepted, 'VDisc'), target.VDisc = accepted.VDisc; end
if isprop(accepted, 'Levels'), target.Levels = accepted.Levels; end
if isprop(accepted, 'Decomposition'), target.Decomposition = accepted.Decomposition; end
if isprop(accepted, 'Cache') && ismethod(accepted, 'reset'), target.reset(); end
scheme.TDisc = workspace;
end

function output = timed(timer,name,operation)
[~,interval]=timer.start(record=name);
guard=onCleanup(@()timer.cancel(interval));
output=operation();
timer.stop(interval=interval);
end
function tf = isText(value)
tf=(ischar(value)&&(isrow(value)||isempty(value))) || ...
    (isstring(value)&&isscalar(value)&&~ismissing(value));
end
