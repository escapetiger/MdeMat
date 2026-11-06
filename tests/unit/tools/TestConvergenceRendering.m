classdef TestConvergenceRendering < matlab.unittest.TestCase
    methods (Test)
        function fittedAndPrescribedOrdersMatchIndependentPowerLaws(test)
            f = figure('Visible', 'off'); guard = onCleanup(@() delete(f)); ax = axes(f);
            h = 2 .^ (0:-1:-5); errors = 3 * h .^ 2.5;
            retained = line(ax, [1, 2], [4, 5]); hold(ax, 'off');
            [artists, fit] = mdemat.tools.visual.plotConvergence(ax, h, errors, ...
                fitted=true, fitRange=[2, 6], label='L2');
            test.verifyEqual(fit.Order, 2.5, AbsTol=3e-14);
            test.verifyEqual(fit.LogAmplitude, log(3), AbsTol=3e-14);
            test.verifyEqual(artists.Curve.XData(:), h(:));
            test.verifyEqual(artists.Curve.YData(:), errors(:));
            test.verifyEqual(artists.Reference.YData(:), errors(2:6).', AbsTol=3e-14);
            test.verifyEqual(ax.XScale, 'log'); test.verifyEqual(ax.YScale, 'log');
            test.verifyFalse(ishold(ax)); test.verifyTrue(isgraphics(retained));
            [specified, fit] = mdemat.tools.visual.plotConvergence(ax, h, errors, ...
                referenceOrder=1, fitRange=[2, 4], colorIndex=2);
            expected = errors(2) * h(2:4) / h(2);
            test.verifyEqual(specified.Reference.XData(:), h(2:4).');
            test.verifyEqual(specified.Reference.YData(:), expected.', AbsTol=3e-14);
            test.verifyEqual(fit.Anchor, [h(2), errors(2)], AbsTol=3e-14);
            test.verifyEqual(fit.ReferenceOrder, 1);
        end
        function zeroAndUnassessedPointsDoNotInventAFit(test)
            f = figure('Visible', 'off'); guard = onCleanup(@() delete(f)); ax = axes(f);
            [artists, fit] = mdemat.tools.visual.plotConvergence(ax, [1, .5, .25], [NaN, 0, 2]);
            test.verifyEqual(artists.Curve.XData, .25); test.verifyEqual(artists.Curve.YData, 2);
            test.verifyTrue(isnan(fit.Order)); test.verifyEmpty(artists.Reference);
            test.verifyError(@() mdemat.tools.visual.plotConvergence(ax, [1, .5], [0, NaN], ...
                fitted=true), 'mdemat:tools:visual:InsufficientFitPoints');
            test.verifyError(@() mdemat.tools.visual.plotConvergence(ax, [1, .5], [1, Inf]), ...
                'mdemat:tools:visual:InvalidConvergenceData');
        end
    end
end
