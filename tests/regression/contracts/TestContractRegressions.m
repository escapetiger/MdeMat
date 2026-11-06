classdef TestContractRegressions < matlab.unittest.TestCase
    % Independent regression checks for configuration and observation contracts.
    methods (Test)
        function parseMathematicsAndEscapedText(testCase)
            file = writeConfig({ ...
                'OPTION=matrix=[1 2;3 4]', ...
                'OPTION=weight=2/3', ...
                'OPTION=trigonometric=sin(pi/4)', ...
                'OPTION=precedence=-2^2', ...
                'OPTION=cells={''all'', ''right''; ''left'', ''''}', ...
                'OPTION=label=''a=b#c'' # actual comment'});
            guard = onCleanup(@() delete(file)); %#ok<NASGU>
            values = mdemat.tools.environment.readConfig(file);
            testCase.verifyEqual(values.matrix, [1, 2; 3, 4]);
            testCase.verifyEqual(values.weight, 2/3);
            testCase.verifyEqual(values.trigonometric, 1/sqrt(2), 'AbsTol', 8*eps);
            testCase.verifyEqual(values.precedence, -4);
            testCase.verifyEqual(values.cells, {'all', 'right'; 'left', ''});
            testCase.verifyEqual(values.label, 'a=b#c');
        end

        function rejectExecutableExpressions(testCase)
            file = writeConfig({'OPTION=value=feval(''sin'', pi)'});
            guard = onCleanup(@() delete(file)); %#ok<NASGU>
            testCase.verifyError(@() mdemat.tools.environment.readConfig(file), ...
                'mdemat:tools:environment:readConfig:UnsupportedExpression');
        end

        function validatedLoadIsAtomic(testCase)
            config = mdemat.tools.environment.Config();
            config = config.define('first', 1);
            config = config.define('positive', 2, @(x) isscalar(x) && x > 0);
            file = writeConfig({'OPTION=first=99', 'OPTION=positive=-1'});
            guard = onCleanup(@() delete(file)); %#ok<NASGU>
            values = mdemat.tools.environment.readConfig(file);
            testCase.verifyError(@() config.update(values), ...
                'mdemat:tools:environment:Config:InvalidValue');
            testCase.verifyEqual(config.get('first'), 1);
            testCase.verifyEqual(config.get('positive'), 2);
        end

        function everyApplicationConfigurationRetainsItsValues(testCase)
            root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
            names = {'hyperbolic', 'parabolic', 'elliptic', 'radiation', 'plasma'};
            methods = {'uwdg', 'ldg', 'ldg', 'mmdg', 'dghs'};
            orders = [2, 3, 3, 2, 2];
            original = path;
            for index = 1:numel(names)
                file = fullfile(root, 'toolbox', 'examples', names{index}, 'config.txt');
                values = mdemat.tools.environment.readConfig(file);
                testCase.verifyEqual(values.schemeName, methods{index});
                testCase.verifyEqual(values.xBasisOrder, orders(index));
                if isfield(values, 'xPenaltyType')
                    testCase.verifySize(values.xPenaltyType, [2, 2]);
                    testCase.verifyEqual(values.xPenaltyType{1, 2}, 'right');
                end
                if strcmp(names{index}, 'plasma')
                    testCase.verifyEqual(values.dealiasingCutoff, 2/3);
                    testCase.verifyEqual(values.tOdeIntName, {'sdirk2', 'sdirk2'});
                end
            end
            testCase.verifyEqual(path, original);
        end

        function gridFactorOutputUsesPointThenElementOrder(testCase)
            grid = mdemat.data.mesh.Grid({[0, 0.25, 0.75, 1]});
            row = grid.collocate({[-0.5, 0, 0.5]}, [2; 3]);
            column = grid.collocate({[-0.5; 0; 0.5]}, [2; 3]);
            expected = [0.25; 0.5; 0.75; 0.75; 0.875; 1];
            testCase.verifyEqual(row{1}, expected);
            testCase.verifyEqual(column{1}, expected);
        end

        function rejectNonintegerConnectivity(testCase)
            testCase.verifyError(@() mdemat.data.mesh.Graph([0; 1], [1, 1.5]), ...
                'mdemat:data:mesh:Graph:InvalidEdgeIndices');
            testCase.verifyError(@() mdemat.data.mesh.SeparableGraph({[0; 1]}, {[1, 1.5]}), ...
                'mdemat:data:mesh:SeparableGraph:InvalidEdgeIndices');
        end

        function structuredCoordinatesRetainComponentLayout(testCase)
            coordinates = {[-1,1],[0,2,4]};
            view = struct('Axes',{coordinates},'Values',[1:6;11:16],'ComponentIds',[1,2]);
            datasets = struct('sample',struct('Type','numeric','Data',struct('u',view)));
            component = mdemat.tools.visual.selectComponent(datasets,'u',2);
            values = component.sample.Data.u;
            testCase.verifyEqual(values.Axes,coordinates);
            testCase.verifyEqual(reshape(values.Values,[2,3]),[11,13,15;12,14,16]);
            testCase.verifyEqual(values.ComponentIds,2);
        end

        function visualSamplingRespectsEachAxis(testCase)
            space = mdemat.common.fem.FiniteElementSpace( ...
                mdemat.data.mesh.UniformGrid([0,3,-2,2],[1,1]), ...
                mdemat.common.fem.BH1OrthotopeElement.modal(2,2));
            plan = mdemat.common.fem.plotSamplingPlan(space,[2;3]);
            testCase.verifyEqual(cellfun(@numel,plan.Axes),[2,3]);
            testCase.verifyEqual(plan.Axes{1},[1,2],'AbsTol',eps);
            testCase.verifyEqual(plan.Axes{2},[-1,0,1],'AbsTol',eps);
            testCase.verifySize(plan.Nodes,[2,6]);
        end

        function componentCountsAndDatasetContentAreExplicit(testCase)
            components = struct('u',[1,2],'v',1);
            data = struct('u',struct('Values',[4;9],'ComponentIds',[1,2],'Coordinates',0), ...
                'v',struct('Values',3,'ComponentIds',1,'Coordinates',0));
            values = struct('Type','numeric','Data',data);
            testCase.verifyEqual(cellfun(@numel,struct2cell(components)).',[2,1]);
            testCase.verifyEqual(values.Type,'numeric');
            testCase.verifyEqual(values.Data.u.Values,[4;9]);
            selected = mdemat.tools.visual.selectComponent(struct('solution',values),'u',2);
            testCase.verifyEqual(selected.solution.Data.u.Values,9);
        end

        function datasetsCanSampleWithoutPlottersOrFigures(testCase)
            figures = findall(groot, 'Type', 'figure');
            currentFigure = get(groot, 'CurrentFigure');
            space = mdemat.common.fem.FiniteElementSpace( ...
                mdemat.data.mesh.UniformGrid([-1,1],1), ...
                mdemat.common.fem.BH1OrthotopeElement.nodal(1,2));
            plan = mdemat.common.fem.plotSamplingPlan(space,2);
            coefficients = space.fit(@(x) 2+3*x);
            fields = mdemat.common.fem.plotFieldSamples( ...
                space,struct('u',coefficients),struct('u',1),plan);
            dataset = struct('Type','numeric','Data',fields);
            testCase.verifyEqual(dataset.Data.u.Values.',[1;3],'AbsTol',4*eps);
            testCase.verifyTrue(isfield(dataset.Data,'u'));
            % Sampling is ready without enabling any graphics resource owner.
            visualizer = mdemat.tools.visual.Visualizer();
            testCase.verifyFalse(visualizer.IsEnabled);
            testCase.verifyEmpty(fieldnames(visualizer.FigureHandles));
            testCase.verifyEqual(findall(groot, 'Type', 'figure'), figures);
            testCase.verifyEqual(get(groot, 'CurrentFigure'), currentFigure);
        end

        function plotterPreservesLayoutDataAndUnrelatedFigures(testCase)
            previousFigure = get(groot, 'CurrentFigure');
            testCase.addTeardown(@() restoreFigure(previousFigure));
            unrelated = figure('Visible', 'off');
            testCase.addTeardown(@() delete(unrelated));
            target = figure('Visible', 'off');
            testCase.addTeardown(@() delete(target));
            set(groot, 'CurrentFigure', unrelated);
            panels = repmat(struct('Mode','custom','RenderFn',@drawSamples, ...
                'Tiles',1,'Spans',[1,1]),1,3);
            panels(1).Spans = [1,2];
            panels(2).Tiles = 3;
            panels(3).Tiles = 4;
            specification = struct('Name','profile','Size',[2,2], ...
                'Panels',panels,'Figure',target);
            session = mdemat.tools.visual.Visualizer(specification,visible="off");
            testCase.addTeardown(@() session.finalize(force=true));
            session.initialize();
            axes = [session.Axes.profile{:}];
            testCase.verifyEqual([session.Layouts.Panels.Tiles],[1,3,4]);
            testCase.verifyEqual(axes(1).Layout.TileSpan, [1, 2]);
            session.render(struct(),title="Owned panels");
            line = findobj(axes(1), 'Type', 'line');
            testCase.verifyEqual(line.XData, [0, 1]);
            testCase.verifyEqual(line.YData, [2, 3]);
            testCase.verifyEqual(string(axes(1).Parent.Title.String),"Owned panels");
            testCase.verifyEqual(get(groot, 'CurrentFigure'), unrelated);
            session.finalize();
            testCase.verifyFalse(any(isgraphics(axes)));
            testCase.verifyTrue(isgraphics(target, 'figure'));
            testCase.verifyTrue(isgraphics(unrelated, 'figure'));
        end

        function figureManagerReplacesOwnedFiguresAndRestoresRootState(testCase)
            previousFigure = get(groot, 'CurrentFigure');
            testCase.addTeardown(@() restoreFigure(previousFigure));
            previousUnits = get(groot, 'Units');
            testCase.addTeardown(@() set(groot, 'Units', previousUnits));
            set(groot, 'Units', 'pixels');
            expected = get(groot, 'ScreenSize');
            set(groot, 'Units', 'inches');
            unrelated = figure('Visible', 'off');
            testCase.addTeardown(@() delete(unrelated));
            panel = struct('Mode','custom','RenderFn',@drawSamples);
            specification = struct('Name','first','Size',[1,1],'Panels',panel);
            session = mdemat.tools.visual.Visualizer(specification,visible="off");
            testCase.addTeardown(@() session.finalize(force=true));
            set(groot,'CurrentFigure',unrelated);
            session.initialize();
            first = session.getFigure('first');
            testCase.verifyEqual(first.Units,'pixels');
            testCase.verifyLessThanOrEqual(first.Position(3:4),expected(3:4));
            testCase.verifyEqual(get(groot,'Units'),'inches');
            specification(2) = specification;
            specification(2).Name = 'second';
            session.setLayouts(specification);
            session.initialize();
            second = [session.getFigure('first'),session.getFigure('second')];
            testCase.verifyFalse(isgraphics(first));
            testCase.verifyTrue(all(isgraphics(second)));
            testCase.verifyTrue(isgraphics(unrelated));
            testCase.verifyEqual(get(groot, 'CurrentFigure'), unrelated);
        end
    end
end

function file = writeConfig(lines)
file = [tempname, '.txt'];
stream = fopen(file, 'wt');
assert(stream >= 0);
guard = onCleanup(@() fclose(stream)); %#ok<NASGU>
fprintf(stream, '%s\n', lines{:});
end

function handles = drawSamples(axes, datasets) %#ok<INUSD>
handles = plot(axes(1), [0, 1], [2, 3]);
end

function restoreFigure(previousFigure)
if isempty(previousFigure) || isgraphics(previousFigure, 'figure')
    set(groot, 'CurrentFigure', previousFigure);
end
end
