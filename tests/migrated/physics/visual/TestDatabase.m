classdef TestDatabase < matlab.unittest.TestCase
    methods (Test)
        function explicitComponentIdsAndNumericFieldNames(test)
            field = struct('Values', [20, 21; 40, 41], 'ComponentIds', [2, 4], ...
                'Coordinates', [1, 0]);
            db.A = struct('Type', 'numeric', 'Data', struct('u2', field));
            out = mdemat.tools.visual.selectComponent(db, 'u2', 4);
            test.verifyEqual(out.A.Data.u2.Values, [40, 41]);
            test.verifyEqual(out.A.Data.u2.ComponentIds, 4);
            test.verifyEqual(out.A.Data.u2.Coordinates, [1, 0]);
        end
        function repeatedDgPointsKeepTheirSeparateValues(test)
            field = struct('Values', [2, 3, 8, 9], 'ComponentIds', 1, ...
                'Coordinates', [0, 1, 1, 2], 'ElementIds', [1, 1, 2, 2], ...
                'Connectivity', [1, 2; 3, 4]);
            db.A = struct('Type', 'numeric', 'Data', struct('u', field));
            out = mdemat.tools.visual.selectComponent(db, 'u', 1);
            test.verifyEqual(out.A.Data.u, field);
        end
        function missingAndEmptyAreExplicit(test)
            [group, missing] = mdemat.tools.visual.selectComponent(struct(), 'u', 1, datasetIds="A");
            test.verifyEmpty(fieldnames(group)); test.verifyEqual(missing.DatasetId, "A");
            test.verifyError(@() mdemat.tools.visual.selectComponent(struct(), 'u', 1, ...
                datasetIds="A", required="A"), 'mdemat:tools:visual:MissingSamples');
            test.verifyError(@() mdemat.tools.visual.selectComponent(struct(), 'u', 1, ...
                required="A"), 'mdemat:tools:visual:MissingSamples');
            empty.A = struct('Type', 'numeric', 'Data', struct('u', ...
                struct('Values', zeros(1, 0), 'ComponentIds', 1, 'Coordinates', zeros(2, 0))));
            out = mdemat.tools.visual.selectComponent(empty, 'u', 1);
            test.verifySize(out.A.Data.u.Values, [1, 0]);
            test.verifyEqual(numel(fieldnames(out)), 1);
        end
    end
end
