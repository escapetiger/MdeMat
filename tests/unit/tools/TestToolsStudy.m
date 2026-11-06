classdef TestToolsStudy < matlab.unittest.TestCase
    methods (Test)
        function boundedStudyValidatesEveryCaseAndRestoresCallerState(testCase)
            root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
            originalPath = path;
            pathGuard = onCleanup(@() path(originalPath)); %#ok<NASGU>
            addpath(fullfile(root,'benchmarks'));
            originalDirectory = pwd;
            originalRng = rng;
            originalFigure = get(groot,'CurrentFigure');
            originalFigures = findall(groot,'Type','figure');
            report = toolsStudy(warmup=0,repetitions=1,denseSize=[8,8], ...
                spatialCells=3,graphEdges=5,seed=17);
            names = string({report.Cases.Name});
            required = ["collection-selection","collection-construction","collection-copy-update", ...
                "fem-plan","fem-samples","weighted-error","dense-mat","dense-json", ...
                "sparse-mat","sparse-json","complex-mat","complex-json", ...
                "small-leaves-mat","small-leaves-json", ...
                "linear-solve-recordHash-0","linear-solve-recordHash-1", ...
                "timer-intervals","metadata-version","metadata-environment", ...
                "graph-batch-0","graph-batch-1"];
            testCase.verifyTrue(all(ismember(required,names)));
            testCase.verifyEqual(numel(unique(names)),numel(names));
            testCase.verifyTrue(all([report.Cases.Validated]));
            testCase.verifyEqual([report.Cases.UntimedExecutions],ones(size(report.Cases)));
            for index = 1:numel(report.Cases)
                item = report.Cases(index);
                testCase.verifySize(item.Timing.Seconds,[1,1]);
                testCase.verifyGreaterThanOrEqual(item.Timing.Seconds,0);
                testCase.verifyFalse(item.MemoryBefore.PeakMeasured);
                testCase.verifyFalse(item.MemoryAfter.PeakMeasured);
            end
            batch = report.Cases(strcmp({report.Cases.Name},'graph-batch-1'));
            separate = report.Cases(strcmp({report.Cases.Name},'graph-batch-0'));
            testCase.verifyEqual(batch.ArtistCount,1);
            testCase.verifyEqual(separate.ArtistCount,5);
            testCase.verifyEqual(sum(isfinite([report.Cases.SerializedBytes])),8);
            small = report.Cases(startsWith(names,'small-leaves-'));
            testCase.verifyEqual([small.LeafCount],[64,64]);
            testCase.verifyEqual(report.ExecutionStatus,'COMPLETED');
            testCase.verifyEqual(rng,originalRng);
            testCase.verifyEqual(pwd,originalDirectory);
            testCase.verifyEqual(get(groot,'CurrentFigure'),originalFigure);
            testCase.verifyEqual(findall(groot,'Type','figure'),originalFigures);
        end
    end
end
