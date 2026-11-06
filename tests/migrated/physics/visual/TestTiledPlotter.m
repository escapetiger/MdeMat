classdef TestTiledPlotter < matlab.unittest.TestCase
    % Layout specifications are values; the session owns all allocated axes.
    methods (Test)
        function spanningPanelsUseRowMajorTiles(test)
            p = struct('Mode', 'custom', 'RenderFn', @drawSamples, 'Spans', [1, 2]);
            panels = repmat(p, 1, 3);
            panels(2).Spans = [1, 1]; panels(3).Spans = [1, 1];
            spec = struct('Name', 'layout', 'Size', [2, 2], 'Panels', panels, ...
                'TileSpacing', 'compact', 'Padding', 'tight');
            vis = mdemat.tools.visual.Visualizer(spec, visible='off');
            cleanup = onCleanup(@() vis.finalize(force=true));
            vis.render(struct(), title='Three panels');
            test.verifyEqual([vis.Layouts.Panels.Tiles], [1, 3, 4]);
            test.verifyEqual(vis.Axes.layout{1}.Layout.TileSpan, [1, 2]);
            test.verifyNumElements(vis.Axes.layout, 3);
            test.verifyEqual(vis.Axes.layout{1}.Parent.TileSpacing, 'compact');
            test.verifyEqual(vis.Axes.layout{1}.Parent.Padding, 'tight');
        end
        function invalidCompleteLayoutCreatesNothing(test)
            before = findall(groot, 'Type', 'figure');
            p = struct('Mode', 'custom', 'RenderFn', @drawSamples, 'Tiles', 1);
            spec = struct('Name', 'bad', 'Size', [1, 1], 'Panels', [p, p]);
            test.verifyError(@() mdemat.tools.visual.Visualizer(spec), ...
                'mdemat:tools:visual:InvalidLayout');
            test.verifyEqual(findall(groot, 'Type', 'figure'), before);
        end
        function deletedAxesAreRecreatedWithoutDeletingBorrowedContent(test)
            f = figure('Visible', 'off'); cleanupF = onCleanup(@() delete(f));
            externalAxes = axes(f); outside = line(externalAxes, [0, 1], [4, 5]);
            p = struct('Mode', 'custom', 'RenderFn', @drawSamples);
            spec = struct('Name', 'layout', 'Size', [1, 1], 'Panels', p, 'Figure', f);
            vis = mdemat.tools.visual.Visualizer(spec, visible='off');
            cleanup = onCleanup(@() vis.finalize(force=true)); vis.render(struct());
            old = vis.Axes.layout{1}; delete(old);
            vis.render(struct());
            test.verifyTrue(all(isgraphics(vis.Axes.layout{1}, 'axes')));
            test.verifyTrue(isgraphics(outside)); test.verifyFalse(isgraphics(old));
        end
    end
end
function h = drawSamples(axes, ~)
h = line(axes(1), [0, 1], [2, 3]);
end
