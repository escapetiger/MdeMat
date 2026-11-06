classdef TestDataset < matlab.unittest.TestCase
    % Value samples replace the former handle; mutation must be explicit.
    methods (Test)
        function coordinatesCannotCollideWithPhysicalFieldNames(test)
            view = struct('Values', [2, 4], 'ComponentIds', 1, 'Coordinates', [0, 1]);
            for field = ["x", "v", "Values"]
                fields = struct(); fields.(field) = view;
                db.RUN = struct('Type', 'numeric', 'Data', fields);
                selected = mdemat.tools.visual.selectComponent(db, field, 1);
                test.verifyEqual(selected.RUN.Data.(field).Coordinates, [0, 1]);
                test.verifyEqual(selected.RUN.Data.(field).Values, [2, 4]);
            end
        end
        function numericSnapshotsAreIndependent(test)
            db.RUN = struct('Type', 'numeric', 'Data', struct('u', ...
                struct('Values', [2, 4], 'ComponentIds', 1, 'Coordinates', [0, 1])));
            copy = db; copy.RUN.Data.u.Values(1) = 20;
            test.verifyEqual(db.RUN.Data.u.Values, [2, 4]);
            test.verifyEqual(copy.RUN.Data.u.Values, [20, 4]);
            test.verifyTrue(isfield(copy.RUN.Data, 'u'));
            test.verifyEqual(fieldnames(copy.RUN.Data), {'u'});
            copy.RUN.Data = struct();
            test.verifyEmpty(fieldnames(copy.RUN.Data));
            test.verifyEqual(copy.RUN.Type, 'numeric');
        end
        function scalarAndGeometryValidation(test)
            bad.RUN = struct('Type', 'numeric', 'Data', repmat(struct('u', []), 1, 2));
            test.verifyError(@() mdemat.tools.visual.selectComponent(bad, 'u', 1), ...
                'mdemat:tools:visual:InvalidDataset');
            view = struct('Values', [1, 2], 'ComponentIds', 1, 'Coordinates', [0, 1]);
            view.Axes = {[0, 1]}; bad.RUN.Data = struct('u', view);
            test.verifyError(@() mdemat.tools.visual.selectComponent(bad, 'u', 1), ...
                'mdemat:tools:visual:InvalidCoordinates');
        end
    end
end
