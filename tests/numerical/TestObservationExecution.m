classdef TestObservationExecution < matlab.unittest.TestCase
    % Numerical execution supplies explicit samples and owns accepted-step identity.
    methods (Test)
        function unequalStepFramesUseTheirPairedTimes(testCase)
            previousFigure = get(groot,'CurrentFigure');
            testCase.addTeardown(@() restoreFigure(previousFigure));
            captured = [];
            panel = struct('Mode','1d','Field','U','ComponentId',1, ...
                'DatasetIds',"NUMERIC",'RenderFn',@capture);
            layout = struct('Name','profile','Size',[1,1],'Panels',panel);
            session = mdemat.tools.visual.Visualizer(layout,visible="off");
            testCase.addTeardown(@() session.finalize(force=true));
            scheme = ObservationSchemeProbe();
            scheme.setObserver('visualizer',session);
            scheme.setObservation('plot',struct('PointsPerDimension',2, ...
                'Components',struct('U',1),'DatasetTypes',struct('NUMERIC','numeric'), ...
                'Times',0.2,'Modes',{{'1d'}},'Title','paired samples'));
            scheme.run(mdemat.common.fem.SpatialState(makeSpace()));
            testCase.verifyEqual(scheme.AcceptedTimes,[0.1,0.3],'AbsTol',4*eps);
            testCase.verifySize(captured,[1,4]);
            testCase.verifyEqual(captured,0.2*ones(1,4),'AbsTol',2e-14);
            plotAxes = session.Axes.profile{1};
            testCase.verifyTrue(contains(string(plotAxes.Parent.Title.String),'0.2'));
            function handles = capture(targetAxes,group)
                view = group.NUMERIC.Data.U;
                captured(end+1,:) = view.Values;
                handles = plot(targetAxes,view.Axes{1},view.Values);
            end
        end

        function physicalXAndVFieldsCannotOverwriteTheirOwnGeometry(testCase)
            space = makeSpace();
            plan = mdemat.common.fem.plotSamplingPlan(space,2);
            originalCoordinates = plan.Coordinates;
            dofs = struct('x',space.fit(@(x) 10+x),'v',space.fit(@(x) 100+2*x));
            fields = mdemat.common.fem.plotFieldSamples(space,dofs,struct('x',1,'v',1),plan);
            datasets = struct('NUMERIC',struct('Type','numeric','Data',fields));
            x = mdemat.tools.visual.selectComponent(datasets,'x',1);
            v = mdemat.tools.visual.selectComponent(datasets,'v',1);
            coordinates = plan.Axes{1};
            testCase.verifyEqual(x.NUMERIC.Data.x.Axes,plan.Axes);
            testCase.verifyEqual(v.NUMERIC.Data.v.Axes,plan.Axes);
            testCase.verifyEqual(x.NUMERIC.Data.x.Values,10+coordinates,'AbsTol',2e-13);
            testCase.verifyEqual(v.NUMERIC.Data.v.Values,100+2*coordinates,'AbsTol',2e-13);
            testCase.verifyEqual(plan.Coordinates,originalCoordinates);
        end

        function failedIntegrationCancelsOnlyItsActiveTimingInterval(testCase)
            scheme = ObservationSchemeProbe();
            timer = mdemat.tools.analyze.Timer(verbose=false);
            scheme.setObserver('timer',timer);
            scheme.FailStep = true;
            state = mdemat.common.fem.SpatialState(makeSpace());
            testCase.verifyError(@() scheme.run(state),'mdemat:test:ObservationStepFailure');
            snapshot = timer.snapshot();
            testCase.verifyFalse(snapshot.IsActive);
            integration = timer.find('Integration');
            testCase.verifyEqual(integration.Duration,0);
            testCase.verifyFalse(integration.IsActive);
            [~,interval] = timer.start(record='retry');
            timer.stop(interval=interval);
            afterRetry = timer.snapshot();
            testCase.verifyFalse(afterRetry.IsActive);
        end

        function automaticCheckpointsKeepDistinctNearbyAcceptedSteps(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            config = modelConfig();
            config.dt = 0.11;
            config.tFinal = 0.14;
            config.ckptTimeStamps = [0.11,0.14];
            config.ckptDir = fixture.Folder;
            config.ckptPrefix = 'model';
            scheme = mdemat.solver.hyperbolic.UwdgScheme(config=config);
            state = scheme.run(mdemat.common.fem.SpatialState(makeSpace()));
            files = dir(fullfile(fixture.Folder,'model_state_step*.mat'));
            testCase.verifyNumElements(files,2);
            testCase.verifyNumElements(unique({files.name}),2);
            times = zeros(1,numel(files));
            steps = zeros(size(times));
            for index = 1:numel(files)
                record = mdemat.internal.execution.readCheckpoint(fullfile(files(index).folder,files(index).name));
                times(index) = record.integrator.timeline.Now;
                match = regexp(files(index).name,'step(\d+)\.mat$','tokens','once');
                steps(index) = str2double(match{1});
                testCase.verifyEqual(steps(index),record.integrator.timeline.Count-1);
                restored = mdemat.common.fem.SpatialState.fromRecord(record.state);
                testCase.verifyEqual(restored.Dofs.U,state.Dofs.U,'AbsTol',2e-13);
            end
            testCase.verifyEqual(sort(times),[0.11,0.14],'AbsTol',4*eps);
            testCase.verifyEqual(sort(steps),[1,2]);
        end

        function oldV1CheckpointIgnoresChangedObservationButRejectsChangedModel(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            file = fullfile(fixture.Folder,'old-v1.json');
            config = modelConfig();
            source = mdemat.solver.hyperbolic.UwdgScheme(config=config);
            state = source.initialize(mdemat.common.fem.SpatialState(makeSpace()));
            state = source.preStep(state);
            state = source.step(state);
            state = source.postStep(state);
            source.saveCheckpoint(file,state);
            record = mdemat.tools.io.readRecord(file);
            record.schemaVersion = 1;
            description = record.configuration;
            description = rmfield(description,'schemaVersion');
            for index = 1:numel(description.callbackNames)
                name = description.callbackNames{index};
                if isvarname(name) && isfield(description.values,name)
                    description.values = rmfield(description.values,name);
                end
            end
            description.values.components = struct('U',1);
            description.values.exact = struct('U',1);
            description.values.errorTypes = {'L1','L2','max'};
            description.values.titlePrefix = 'old observation';
            record.configuration = description;
            mdemat.tools.io.writeRecord(file,record,overwrite=true);
            changed = config;
            changed.components = struct('U',[]);
            changed.exact = @(x,varargin) ones(1,size(x,2));
            changed.errorTypes = {'LInf'};
            resumed = mdemat.solver.hyperbolic.UwdgScheme(config=changed);
            resumed.setObservation('plot',struct('Title','new observation'));
            restored = resumed.loadCkptState(file,mdemat.common.fem.SpatialState(makeSpace()));
            testCase.verifyEqual(restored.Dofs,state.Dofs);
            testCase.verifyEqual(resumed.TDisc.Timeline.Now,source.TDisc.Timeline.Now);
            testCase.verifyEqual(resumed.TDisc.U0,source.TDisc.U0);
            changed.advection = 1;
            differentModel = mdemat.solver.hyperbolic.UwdgScheme(config=changed);
            workspace = differentModel.TDisc;
            testCase.verifyError(@() differentModel.loadCkptState(file,restored), ...
                'mdemat:persistence:InvalidCheckpoint');
            testCase.verifySameHandle(differentModel.TDisc,workspace);
        end
    end
end

function space = makeSpace()
space = mdemat.common.fem.FiniteElementSpace(mdemat.data.mesh.UniformGrid([0,1],2), ...
    mdemat.common.fem.BH1OrthotopeElement.nodal(1,2));
end

function config = modelConfig()
config = struct('advection',0,'ic',@(x) ones(1,size(x,2)), ...
    'dt',0.1,'cfl',[],'tFinal',0.2,'components',struct('U',1));
end

function restoreFigure(previous)
if isempty(previous) || isgraphics(previous,'figure'), set(groot,'CurrentFigure',previous); end
end
