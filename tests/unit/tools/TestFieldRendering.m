classdef TestFieldRendering < matlab.unittest.TestCase
    methods (Test)
        function lineLabelsLogScaleAndComponentIdentity(test)
            f = figure('Visible', 'off'); guard = onCleanup(@() delete(f)); ax = axes(f);
            field = struct('Values', [2, 4, 9; 20, 40, 90], 'ComponentIds', [2, 4]);
            field.Axes = {[1, 2, 3]}; datasets.EXACT = struct('Type', 'exact', 'Data', struct('u', field));
            style = mdemat.tools.visual.plotStyle(logX=true, logY=true);
            h = mdemat.tools.visual.plotField(ax, datasets, 'u', 2, style=style);
            curve = h(isgraphics(h, 'line')); test.verifyEqual(curve.YData, [2, 4, 9]);
            test.verifyEqual(ax.XScale, 'log'); test.verifyEqual(ax.YScale, 'log');
            delete(h(isgraphics(h)));
            mdemat.tools.visual.plotField(ax, datasets, 'u', 4);
            test.verifyEqual(ax.XScale, 'linear'); test.verifyEqual(ax.YScale, 'linear');
            curve = findall(ax, 'Type', 'line'); test.verifyEqual(curve.YData, [20, 40, 90]);
            test.verifyFalse(ishold(ax));
        end
        function nonuniformNodesRetainPhysicalPositionsAndLastRow(test)
            f = figure('Visible', 'off'); guard = onCleanup(@() delete(f));
            ax = [subplot(1, 2, 1, 'Parent', f), subplot(1, 2, 2, 'Parent', f)];
            x = [-2, 0, 1]; y = [0, 1, 3, 5]; [X, Y] = ndgrid(x, y); U = X + 3 * Y;
            db = database({x, y}, U);
            style = mdemat.tools.visual.plotStyle(symmetricColorLimits=true);
            h = mdemat.tools.visual.plotField(ax, db, 'u', 1, mode='2d', style=style);
            h = h(isgraphics(h, 'surface'));
            test.verifyNumElements(h, 2);
            test.verifyEqual(h(1).XData, X); test.verifyEqual(h(1).YData, Y);
            test.verifyEqual(h(1).CData, U); test.verifyEqual(h(2).CData, 2 * U);
            test.verifyEqual(ax(1).CLim, [-32, 32]); test.verifyEqual(ax(2).CLim, [-32, 32]);
        end
        function volumeSlicesMatchAnIndependentAffineField(test)
            f = figure('Visible', 'off'); guard = onCleanup(@() delete(f));
            ax = [subplot(1, 2, 1, 'Parent', f), subplot(1, 2, 2, 'Parent', f)];
            x = [-1, 0, 1]; y = [-2, 0, 2]; z = [-3, 0, 3];
            [X, Y, Z] = ndgrid(x, y, z); U = X + 2 * Y + 3 * Z;
            mdemat.tools.visual.plotField(ax, database({x, y, z}, U), 'u', 1, mode='3d');
            surfaces = findall(ax(1), 'Type', 'surface'); test.verifyNumElements(surfaces, 3);
            for item = surfaces'
                test.verifyEqual(item.CData, item.XData + 2 * item.YData + 3 * item.ZData, AbsTol=2e-14);
            end
            test.verifyEqual(ax(1).CLim, ax(2).CLim);
        end
        function diagonalCoversRectangularPhysicalDomain(test)
            f = figure('Visible', 'off'); guard = onCleanup(@() delete(f)); ax = axes(f);
            for nd = 2:3
                cla(ax); coordinates = {[-1, 0, 2], [-2, 1, 4], [-3, 0, 6]};
                coordinates = coordinates(1:nd); grids = cell(1, nd);
                [grids{:}] = ndgrid(coordinates{:}); U = zeros(size(grids{1}));
                for d = 1:nd, U = U + d * grids{d}; end
                db = database(coordinates, U); db = rmfield(db, 'DATA2');
                mdemat.tools.visual.plotField(ax, db, 'u', 1, mode='slice1d');
                curve = findall(ax, 'Type', 'line'); test.verifyNumElements(curve, 1);
                lower = cellfun(@min, coordinates); upper = cellfun(@max, coordinates);
                length = norm(upper - lower); r = curve.XData;
                expected = (1:nd) * ((lower + upper)' / 2) + ...
                    r * ((1:nd) * (upper - lower)') / length;
                test.verifyEqual(curve.YData, expected, AbsTol=4e-13);
                test.verifyEqual(r([1, end]), [-length / 2, length / 2], AbsTol=2e-14);
            end
        end
        function palettesAndExplicitComplexRenderingRemainAvailable(test)
            for name = ["bright", "muted", "vibrant", "light", "high_contrast", "filbet"]
                style = mdemat.tools.visual.plotStyle(discreteColorScheme=name);
                test.verifyNotEmpty(style.Palette);
            end
            for name = ["sunset", "burd", "iridescent"]
                style = mdemat.tools.visual.plotStyle(continuousColorScheme=name);
                test.verifySize(style.ColorMap, [size(style.ColorMap, 1), 3]);
                test.verifyGreaterThanOrEqual(min(style.ColorMap(:)), 0);
                test.verifyLessThanOrEqual(max(style.ColorMap(:)), 1);
            end
            cyclic = mdemat.tools.visual.plotStyle(lineScheme={'-', '--'}, scatterScheme={'o', 's'});
            line = mdemat.tools.visual.curveStyle(cyclic, 10, 4);
            points = mdemat.tools.visual.curveStyle(cyclic, 10, 4, kind='scatter');
            test.verifyEqual(line.Color, '#EE6677'); test.verifyEqual(line.LineStyle, '--');
            test.verifyEqual(points.Marker, 's'); test.verifyEqual(points.LineStyle, 'none');
            f = figure('Visible', 'off'); guard = onCleanup(@() delete(f)); ax = axes(f);
            db = database({[0, 1]}, [1+2i, 3+4i]); db = rmfield(db, 'DATA2');
            test.verifyError(@() mdemat.tools.visual.plotField(ax, db, 'u', 1), ...
                'mdemat:tools:visual:ComplexDisplayRequired');
            style = mdemat.tools.visual.plotStyle(complexDisplay='imag');
            h = mdemat.tools.visual.plotField(ax, db, 'u', 1, style=style);
            line = h(isgraphics(h, 'line')); test.verifyEqual(line.YData, [2, 4]);
        end
        function dgConnectivityDoesNotConnectSeparateTraces(test)
            f = figure('Visible', 'off'); guard = onCleanup(@() delete(f)); ax = axes(f);
            field = struct('Values', [2, 3, 8, 9], 'ComponentIds', 1, ...
                'Coordinates', [0, 1, 1, 2], 'Connectivity', [1, 2; 3, 4]);
            db.DG = struct('Type', 'fixed', 'Data', struct('u', field));
            h = mdemat.tools.visual.plotField(ax, db, 'u', 1);
            curve = h(isgraphics(h, 'line'));
            test.verifyEqual(curve.XData, [0, 1, NaN, 1, 2, NaN]);
            test.verifyEqual(curve.YData, [2, 3, NaN, 8, 9, NaN]);
        end
        function sparseDisplayDensifiesOnlyTheRequestedComponent(test)
            f = figure('Visible', 'off'); guard = onCleanup(@() delete(f)); ax = axes(f);
            [X, Y] = ndgrid([0, 0.2, 1], [-2, 2]); U = X + 2 * Y;
            field = struct('Values', sparse([U(:).'; 1e6 * ones(1, numel(U))]), ...
                'ComponentIds', [2, 4]);
            field.Axes = {[0, 0.2, 1], [-2, 2]};
            db.RUN = struct('Type', 'numeric', 'Data', struct('u', field));
            h = mdemat.tools.visual.plotField(ax, db, 'u', 2, mode='2d');
            surface = h(isgraphics(h, 'surface'));
            test.verifyEqual(surface.CData, U);
            test.verifyTrue(issparse(db.RUN.Data.u.Values));
            test.verifyEqual(db.RUN.Data.u.Values(2, :), sparse(1e6 * ones(1, numel(U))));
        end
        function failedRenderingRollsBackOnlyNewArtists(test)
            f = figure('Visible', 'off'); guard = onCleanup(@() delete(f)); ax = axes(f);
            original = line(ax, [0, 1], [4, 5]); hold(ax, 'on');
            db = database({[0, 1]}, [1, 2]); db.DATA2.Data.u.Values = [1, Inf];
            test.verifyError(@() mdemat.tools.visual.plotField(ax, db, 'u', 1), ...
                'mdemat:tools:visual:InfiniteValues');
            test.verifyTrue(ishold(ax)); test.verifyTrue(isgraphics(original));
            test.verifyNumElements(findall(ax, 'Type', 'line'), 1);
        end
    end
end
function db = database(coordinates, values)
db = struct();
for factor = 1:2
    field = struct('Values', factor * values(:).', 'ComponentIds', 1);
    field.Axes = coordinates;
    db.(sprintf('DATA%d', factor)) = struct('Type', 'exact', 'Data', struct('u', field));
end
end
