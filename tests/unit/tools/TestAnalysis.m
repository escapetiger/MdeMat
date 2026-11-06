classdef TestAnalysis < matlab.unittest.TestCase
    methods (Test)
        function physicalWeightedNormsMatchHandCalculation(testCase)
            reference = [1, 2, 100; 0, 0, -100];
            values = reference + [1, -2, 999; 2, 1, -999];
            report = mdemat.tools.analyze.weightedError(values, reference, [0.25; 0.75; 0]);
            testCase.verifyEqual(report.L1, 3);
            testCase.verifyEqual(report.L2, sqrt(5), 'AbsTol', 8*eps);
            testCase.verifyEqual(report.LInf, 2);
            testCase.verifyEqual(report.ComponentL2, sqrt([3.25; 1.75]), 'AbsTol', 8*eps);
            testCase.verifyEqual(report.ReferenceL2, sqrt(3.25), 'AbsTol', 8*eps);
            testCase.verifyEqual(report.RelativeL2, sqrt(5/3.25), 'AbsTol', 8*eps);
            testCase.verifyTrue(report.RelativeDefined);
        end

        function complexFieldsUseMagnitudeAndZeroReferenceIsExplicit(testCase)
            report = mdemat.tools.analyze.weightedError([1i, 2], [0, 0], [0.5, 0.5]);
            testCase.verifyEqual(report.L2, sqrt(2.5), 'AbsTol', 8*eps);
            testCase.verifyFalse(report.RelativeDefined);
            testCase.verifyTrue(isnan(report.RelativeL2));
            small = mdemat.tools.analyze.weightedError([1, 1], [1e-12, 1e-12], ...
                [0.5, 0.5], referenceFloor=1e-10);
            testCase.verifyFalse(small.RelativeDefined);
            testCase.verifyError(@() mdemat.tools.analyze.weightedError(ones(2), ...
                ones(2), [0, 0]), 'mdemat:analysis:InvalidShape');
            masked = mdemat.tools.analyze.weightedError([1, realmax], [0, -realmax], [1, 0]);
            testCase.verifyEqual(masked.L2, 1);
            overflow = mdemat.tools.analyze.weightedError(realmax, -realmax, 1);
            testCase.verifyFalse(overflow.FiniteNorms);
            testCase.verifyFalse(overflow.RelativeDefined);
        end

        function rateStudyPreservesZeroAndDivergingCases(testCase)
            scales = [0.5; 0.25; 0.125];
            errors = [scales.^3, [1; 0; 0], scales.^-1];
            report = mdemat.tools.analyze.convergenceRates(scales, errors);
            testCase.verifyEqual(report.LocalRates(:, 1), [3; 3], 'AbsTol', 8*eps);
            testCase.verifyEqual(report.LocalRates(:, 3), [-1; -1], 'AbsTol', 8*eps);
            testCase.verifyFalse(any(report.Defined(:, 2)));
            testCase.verifyTrue(all(isnan(report.LocalRates(:, 2))));
            testCase.verifyError(@() mdemat.tools.analyze.convergenceRates( ...
                [0.5; 0.5], [1; 0.1]), 'mdemat:analysis:InvalidRefinement');
        end

        function nearbyScalesAndWeightedExtremeDifferencesRemainRepresentable(testCase)
            h = [1e200;1e200*(1-eps)];
            rates = mdemat.tools.analyze.convergenceRates(h,[2;1]);
            expected = log(2)/log1p((h(1)-h(2))/h(2));
            testCase.verifyEqual(rates.LocalRates,expected,RelTol=8*eps);
            testCase.verifyTrue(rates.Finite);
            report = mdemat.tools.analyze.weightedError(realmax,-realmax,1e-300);
            testCase.verifyTrue(isfinite(report.L2));
            testCase.verifyEqual(report.L2,2*(realmax*1e-150),RelTol=3e-13);
            testCase.verifyTrue(isinf(report.LInf));
            complexReport = mdemat.tools.analyze.weightedError(complex(realmax,realmax),0,1e-300);
            testCase.verifyEqual(complexReport.L1,sqrt(2)*(realmax*1e-300),RelTol=3e-13);
            testCase.verifyEqual(complexReport.L2,sqrt(2)*(realmax*1e-150),RelTol=3e-13);
            testCase.verifyTrue(isinf(complexReport.LInf));
            ratio = mdemat.tools.analyze.weightedError(1,1e-310,1);
            testCase.verifyFalse(ratio.RelativeDefined);
            testCase.verifyFalse(ratio.RelativeFinite);
            testCase.verifyTrue(isnan(ratio.RelativeL2));
        end

        function unassessedSeriesNeverInventsAnErrorOrRate(testCase)
            first = struct('Scale',0.5,'Field','U','ComponentIds',1,'SpatialNorm','L2', ...
                'ComponentNorm','none','Comparison','unassessed','Reference','','Error',NaN);
            second = first;
            second.Scale = 0.25;
            summary = mdemat.tools.analyze.summarizeConvergence([first,second]);
            testCase.verifyTrue(all(isnan(summary.Errors)));
            testCase.verifyTrue(isnan(summary.LocalRates));
            testCase.verifyFalse(summary.Defined);
            first.Comparison = 'exact';
            testCase.verifyError(@() mdemat.tools.analyze.summarizeConvergence(first), ...
                'mdemat:analysis:InvalidObservation');
        end

        function allModeSequenceRetainsOneCurveAndEveryActualCutoff(testCase)
            scales = [0.5;0.25;0.125];
            resolved = {1:2;1:4;1:8};
            observations = repmat(struct('Scale',1,'Field','D','ComponentIds','all', ...
                'ResolvedComponentIds',[],'SpatialNorm','L2','ComponentNorm','sequence', ...
                'Comparison','interlevel','Reference','coarse:D','Error',0),1,3);
            for level = 1:3
                % Constants on a domain of measure 2: every component has
                % L2 error sqrt(2)*3*h^2/sqrt(N), independently of cutoff N.
                modes = numel(resolved{level});
                values = repmat(3*scales(level)^2/sqrt(modes),modes,2);
                report = mdemat.tools.analyze.weightedError(values,zeros(size(values)),[0.5,1.5]);
                observations(level).Scale = scales(level);
                observations(level).ResolvedComponentIds = resolved{level};
                observations(level).Error = norm(report.ComponentL2,2);
            end
            summary = mdemat.tools.analyze.summarizeConvergence(observations);
            testCase.verifyNumElements(summary,1);
            testCase.verifyEqual(summary.Identity.ComponentIds,'all');
            testCase.verifyEqual(summary.ResolvedComponentIds,resolved);
            testCase.verifyEqual(summary.Errors,3*sqrt(2)*scales.^2,AbsTol=2e-15);
            testCase.verifyEqual(summary.LocalRates,[2;2],AbsTol=8*eps);
        end

        function fixedSelectionOrderAndAllSelectionHaveDistinctIdentities(testCase)
            first = struct('Scale',0.5,'Field','D','ComponentIds',[2,4], ...
                'ResolvedComponentIds',[2,4],'SpatialNorm','L2','ComponentNorm','sequence', ...
                'Comparison','exact','Reference','supplied:D','Error',0.25);
            reversed = first;
            reversed.ComponentIds = [4,2];
            reversed.ResolvedComponentIds = [4,2];
            allSelection = first;
            allSelection.ComponentIds = "all";
            allSelection.ResolvedComponentIds = 1:4;
            summary = mdemat.tools.analyze.summarizeConvergence([first,reversed,allSelection]);
            testCase.verifyNumElements(summary,3);
            testCase.verifyEqual(summary(1).Identity.ComponentIds,[2,4]);
            testCase.verifyEqual(summary(2).Identity.ComponentIds,[4,2]);
            testCase.verifyEqual(summary(3).Identity.ComponentIds,'all');
            numericOnly = rmfield(first,'ResolvedComponentIds');
            legacy = mdemat.tools.analyze.summarizeConvergence(numericOnly);
            testCase.verifyEqual(legacy.ResolvedComponentIds,{[2,4]});
        end

        function allSelectionRequiresAnExplicitCompleteOrderedResolution(testCase)
            observation = struct('Scale',0.5,'Field','D','ComponentIds','all', ...
                'ResolvedComponentIds',1:4,'SpatialNorm','L2','ComponentNorm','sequence', ...
                'Comparison','exact','Reference','supplied:D','Error',0.25);
            testCase.verifyError(@() mdemat.tools.analyze.summarizeConvergence( ...
                rmfield(observation,'ResolvedComponentIds')),'mdemat:analysis:InvalidObservation');
            invalid = observation;
            invalid.ComponentNorm = 'L2';
            testCase.verifyError(@() mdemat.tools.analyze.summarizeConvergence(invalid), ...
                'mdemat:analysis:InvalidObservation');
            for ids = {[],[1,3],[2,1],[1,1]}
                invalid = observation;
                invalid.ResolvedComponentIds = ids{1};
                testCase.verifyError(@() mdemat.tools.analyze.summarizeConvergence(invalid), ...
                    'mdemat:analysis:InvalidObservation');
            end
            invalid = observation;
            invalid.ComponentIds = [2,4];
            invalid.ResolvedComponentIds = [4,2];
            testCase.verifyError(@() mdemat.tools.analyze.summarizeConvergence(invalid), ...
                'mdemat:analysis:InvalidObservation');
        end

        function timingDeclaresWorkCountAndNativeScope(testCase)
            count = 0;
            rngBefore = rng;
            report = mdemat.tools.analyze.benchmark(@operation, ...
                warmup=2, repetitions=4, scope='test assembly');
            testCase.verifyEqual(count, 6);
            testCase.verifySize(report.Seconds, [4, 1]);
            testCase.verifyGreaterThanOrEqual(report.Seconds, zeros(4, 1));
            testCase.verifyEqual(report.MedianSeconds, median(report.Seconds));
            testCase.verifyEqual(report.Scope, 'test assembly');
            testCase.verifyEqual(report.Process, 'current MATLAB');
            testCase.verifyEqual(rng, rngBefore);
            function operation()
                count = count + 1;
            end
        end

        function memoryIsSnapshotOfMatlabAndDoesNotClaimPeak(testCase)
            report = mdemat.tools.analyze.memorySnapshot();
            testCase.verifyEqual(report.Kind, 'snapshot');
            testCase.verifyFalse(report.PeakMeasured);
            testCase.verifyTrue(isnan(report.PeakBytes));
            testCase.verifyEqual(report.Process, 'current MATLAB');
            if ispc
                testCase.verifyTrue(report.Available, report.Reason);
                testCase.verifyGreaterThan(report.Bytes, 0);
            else
                testCase.verifyFalse(report.Available);
                testCase.verifyTrue(isnan(report.Bytes));
            end
        end
    end
end
