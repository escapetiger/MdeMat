classdef TestVisualizer < matlab.unittest.TestCase
    methods (Test)
        function emptySessionDoesNotCreateFigures(test)
            before = findall(groot, 'Type', 'figure'); current = get(groot, 'CurrentFigure');
            vis = mdemat.tools.visual.Visualizer();
            vis.render(struct());
            test.verifyFalse(vis.IsEnabled); test.verifyEmpty(fieldnames(vis.FigureHandles));
            test.verifyEqual(findall(groot, 'Type', 'figure'), before);
            test.verifyEqual(get(groot, 'CurrentFigure'), current);
        end
        function fixedBindingsDoNotShiftWhenAReferenceIsMissing(test)
            panel = struct('Mode', '2d', 'Field', 'u', 'ComponentId', 4, ...
                'DatasetIds', ["REFERENCE", "NUMERIC"]);
            layout = struct('Name', 'field', 'Size', [1, 2], 'Panels', panel);
            vis = mdemat.tools.visual.Visualizer(layout, visible='off');
            guard = onCleanup(@() vis.finalize(force=true));
            db = samples(); db.REFERENCE = db.NUMERIC;
            vis.render(db); axes = vis.Axes.field{1};
            test.verifyTrue(contains(axes(1).Title.String, 'REFERENCE'));
            test.verifyTrue(contains(axes(2).Title.String, 'NUMERIC'));
            vis.render(rmfield(db, 'REFERENCE'));
            test.verifyTrue(contains(axes(1).Title.String, 'no data'));
            test.verifyTrue(contains(axes(2).Title.String, 'NUMERIC'));
            image = findall(axes(2), 'Type', 'surface');
            test.verifyEqual(image.CData, [40, 42; 41, 43]);
            vis.render(rmfield(db, 'REFERENCE'));
            test.verifyNumElements(findall(axes(2), 'Type', 'surface'), 1);
        end
        function failedCustomCallbackReleasesItsAxisAndCanRetry(test)
            external = figure('Visible', 'off'); guardExternal = onCleanup(@() delete(external));
            outsideAxes = axes(external); outside = line(outsideAxes, [0, 1], [3, 4]);
            panel = struct('Mode', 'custom', 'RenderFn', @failAfterDrawing);
            layout = struct('Name', 'custom', 'Size', [1, 1], 'Panels', panel, 'Figure', external);
            vis = mdemat.tools.visual.Visualizer(layout, visible='off');
            guard = onCleanup(@() vis.finalize(force=true));
            test.verifyError(@() vis.render(struct()), 'mdemat:test:InjectedRenderFailure');
            test.verifyFalse(vis.InRender); test.verifyTrue(isgraphics(outside));
            test.verifyFalse(any(isgraphics(vis.Axes.custom{1})));
            layout.Panels.RenderFn = @drawSamples;
            vis.setLayouts(layout); vis.render(struct());
            test.verifyNumElements(findall(vis.Axes.custom{1}, 'Type', 'line'), 1);
            test.verifyTrue(isgraphics(outside));
        end
        function customCallbackAcceptsDataPlotHandleLayers(test)
            graph = mdemat.data.mesh.Graph([0, 0; 1, 1], [1, 2]);
            p = struct('Mode', 'custom', 'RenderFn', @(aa, ~) plot(graph, axes=aa(1)));
            layout = struct('Name', 'graph', 'Size', [1, 1], 'Panels', p);
            vis = mdemat.tools.visual.Visualizer(layout, visible='off');
            guard = onCleanup(@() vis.finalize(force=true));
            vis.render(struct()); vis.render(struct());
            test.verifyNumElements(findall(vis.Axes.graph{1}, 'Type', 'line'), 1);
            test.verifyNumElements(findall(vis.Axes.graph{1}, 'Type', 'scatter'), 1);
        end
        function callbackCannotClaimAnotherAxesArtist(test)
            f = figure('Visible', 'off'); guardF = onCleanup(@() delete(f));
            ax = axes(f); outside = line(ax, [0, 1], [4, 5]);
            panel = struct('Mode', 'custom', 'RenderFn', @(varargin) outside);
            layout = struct('Name', 'custom', 'Size', [1, 1], 'Panels', panel);
            vis = mdemat.tools.visual.Visualizer(layout, visible='off');
            guard = onCleanup(@() vis.finalize(force=true));
            test.verifyError(@() vis.render(struct()), 'mdemat:tools:visual:InvalidCallbackArtists');
            test.verifyTrue(isgraphics(outside)); test.verifyFalse(vis.InRender);
        end
        function finalizeAndVetoHaveAClosedLifecycle(test)
            panel = struct('Mode', '1d', 'Field', 'u', 'ComponentId', 2, 'DatasetIds', "NUMERIC");
            layout = struct('Name', 'profile', 'Size', [1, 1], 'Panels', panel);
            vis = mdemat.tools.visual.Visualizer(layout, visible='off');
            guard = onCleanup(@() vis.finalize(force=true));
            db = samples(); db.NUMERIC.Data.u.Axes = {[0, 1, 2, 3]};
            vis.render(db); figure = vis.getFigure('profile');
            figure.CloseRequestFcn = @(varargin) [];
            report = vis.finalize();
            test.verifyEqual(report.Vetoed, "profile"); test.verifyTrue(isgraphics(figure));
            test.verifyEqual(vis.getFigure('profile'), figure);
            test.verifyError(@() vis.render(db), 'mdemat:tools:visual:FinalizedSession');
            test.verifyError(@() vis.reset(), 'mdemat:tools:visual:CloseVetoed');
            vis.reset(force=true); test.verifyFalse(isgraphics(figure));
            vis.render(db); test.verifyTrue(isgraphics(vis.getFigure('profile')));
        end
        function closingAnOwnedFigureRebuildsTheSession(test)
            panel = struct('Mode', '1d', 'Field', 'u', 'ComponentId', 2, 'DatasetIds', "NUMERIC");
            layout = struct('Name', 'profile', 'Size', [1, 1], 'Panels', panel);
            vis = mdemat.tools.visual.Visualizer(layout, visible='off');
            guard = onCleanup(@() vis.finalize(force=true));
            db = samples(); db.NUMERIC.Data.u.Axes = {[0, 1, 2, 3]};
            vis.render(db); previous = vis.getFigure('profile'); delete(previous);
            vis.render(db); test.verifyTrue(isgraphics(vis.getFigure('profile')));
            test.verifyFalse(isgraphics(previous));
        end
    end
end
function h = failAfterDrawing(axes, ~)
h = line(axes(1), [0, 1], [2, 3]);
error('mdemat:test:InjectedRenderFailure', 'Fail after allocating an owned artist.');
end
function h = drawSamples(axes, ~)
h = line(axes(1), [0, 1], [2, 3]);
end
function db = samples()
view = struct('Values', [20, 21, 22, 23; 40, 41, 42, 43], 'ComponentIds', [2, 4]);
view.Axes = {[0, 1], [0, 1]};
db.NUMERIC = struct('Type', 'numeric', 'Data', struct('u', view));
end
