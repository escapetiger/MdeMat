classdef TestFigureManager < matlab.unittest.TestCase
    % Figure ownership belongs to the session rather than a second manager.
    methods (Test)
        function borrowedFigureAndRootStateSurvive(test)
            previousUnits = get(groot, 'Units');
            unitsGuard = onCleanup(@() set(groot, 'Units', previousUnits));
            set(groot, 'Units', 'inches');
            external = figure('Visible', 'off'); unrelated = figure('Visible', 'off');
            guard = onCleanup(@() delete([external, unrelated]));
            ax = axes(external); retained = line(ax, [0, 1], [3, 4]);
            panel = struct('Mode', 'custom', 'RenderFn', @drawSamples);
            layout = struct('Name', 'panels', 'Size', [1, 1], 'Panels', panel, 'Figure', external);
            vis = mdemat.tools.visual.Visualizer(layout, visible='off');
            cleanup = onCleanup(@() vis.finalize(force=true));
            set(groot, 'CurrentFigure', unrelated);
            vis.render(struct(), title='Owned panels');
            test.verifyEqual(get(groot, 'CurrentFigure'), unrelated);
            test.verifyEqual(get(groot, 'Units'), 'inches');
            test.verifyTrue(isgraphics(retained));
            figure(vis.getFigure('panels'));
            test.verifyEqual(get(groot, 'CurrentFigure'), external);
            vis.finalize(force=true);
            test.verifyTrue(isgraphics(external)); test.verifyTrue(isgraphics(retained));
            test.verifyTrue(isgraphics(unrelated));
        end
        function allOwnedFiguresAreRegisteredAndClosed(test)
            p = struct('Mode', 'custom', 'RenderFn', @drawSamples);
            specs = repmat(struct('Name', '', 'Size', [1, 1], 'Panels', p), 1, 3);
            for k = 1:3, specs(k).Name = sprintf('plot%d', k); end
            vis = mdemat.tools.visual.Visualizer(specs, visible='off');
            cleanup = onCleanup(@() vis.finalize(force=true)); vis.initialize();
            figures = struct2cell(vis.FigureHandles);
            test.verifyNumElements(figures, 3);
            for k = 1:3
                test.verifyTrue(isgraphics(figures{k}, 'figure'));
                figure(figures{k});
                test.verifyEqual(get(groot, 'CurrentFigure'), figures{k});
            end
            vis.finalize();
            for k = 1:3, test.verifyFalse(isgraphics(figures{k})); end
            test.verifyEmpty(fieldnames(vis.FigureHandles));
        end
    end
end
function h = drawSamples(axes, ~)
h = line(axes(1), [0, 1], [2, 3]);
end
