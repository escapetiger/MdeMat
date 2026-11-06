classdef TestRunRecord < matlab.unittest.TestCase
    methods (Test)
        function successfulExecutionWithoutAssessmentRemainsUnverified(testCase)
            rngBefore = rng;
            record = mdemat.tools.io.createRunRecord(struct('U', [1, 2]), ...
                struct('Converged', true, 'SolveSeconds', 0.01));
            testCase.verifyEqual(record.ValidationReport.Status, 'UNVERIFIED');
            testCase.verifyEmpty(record.Provenance.UTC);
            testCase.verifyEmpty(record.Provenance.Seed);
            testCase.verifyEmpty(record.Provenance.CodeVersion);
            testCase.verifyFalse(record.Environment.InstalledProductsIncluded);
            testCase.verifyEmpty(record.Environment.InstalledProducts);
            testCase.verifyEqual(rng, rngBefore);
            inventory = mdemat.tools.environment.describe(installedProducts=true);
            testCase.verifyTrue(inventory.InstalledProductsIncluded);
            testCase.verifyTrue(any(strcmp({inventory.InstalledProducts.Name}, 'MATLAB')));
            testCase.verifyFalse(isfield(inventory, 'License'));
            testCase.verifyFalse(isfield(inventory, 'Hostname'));
        end

        function explicitReportsAndPackageVersionsSurviveNumericRecordIO(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            file = fullfile(fixture.Folder, 'run.json');
            validation = struct('Status', "PASS", 'Criterion', "weighted L2 <= 1e-8", ...
                'Measurements', struct('L2', 1e-10), 'Reference', "analytic quadratic", ...
                'Bounds', struct('L2', 1e-8), 'Reason', "");
            record = mdemat.tools.io.createRunRecord(struct('U', [1, 2i]), ...
                struct('Status', 'COMPLETED'), validation=validation, ...
                utc='2026-11-05T10:00:00Z', seed=uint32(7), codeVersion='explicit-commit');
            mdemat.tools.io.writeRecord(file, record);
            [restored, info] = mdemat.tools.io.readRecord(file);
            identity = mdemat.tools.environment.versionInfo();
            testCase.verifyEqual(restored, record);
            testCase.verifyEqual(info.ToolboxVersion, identity.ToolboxVersion);
            testCase.verifyEqual(record.ToolboxVersion, info.ToolboxVersion);
            testCase.verifyEqual(record.NumericRecordSchemaVersion, info.SchemaVersion);
            testCase.verifyEqual(record.SchemaVersion, identity.RunRecordSchemaVersion);
            invalid = validation; invalid.Bounds = struct();
            testCase.verifyError(@() mdemat.tools.io.createRunRecord(struct(), struct(), ...
                validation=invalid), 'mdemat:io:InvalidValidationReport');
            testCase.verifyError(@() mdemat.tools.io.createRunRecord(struct(), struct(), ...
                utc='local morning'), 'mdemat:io:InvalidProvenance');
            testCase.verifyError(@() mdemat.tools.io.createRunRecord(struct(), struct(), ...
                seed=-1), 'mdemat:io:InvalidProvenance');
            testCase.verifyError(@() mdemat.tools.io.createRunRecord(struct('f', @sin), struct()), ...
                'mdemat:io:UnsupportedValue');
            record.Solution.callback = @sin;
            testCase.verifyError(@() mdemat.tools.io.writeRecord(file,record,overwrite=true), ...
                'mdemat:io:UnsupportedValue');
            testCase.verifyEqual(mdemat.tools.io.readRecord(file),restored);
        end

        function checkoutDiscoveryRejectsInstalledLikeLayout(testCase)
            root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
            actual = mdemat.tools.environment.projectRoot();
            testCase.verifyTrue(strcmpi(actual, root));
            testCase.verifyTrue(isfile(fullfile(actual, 'MdeMat.prj')));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            installed = fullfile(fixture.Folder, 'Toolboxes', 'MdeMat', '+mdemat', ...
                '+tools', '+environment');
            mkdir(installed);
            testCase.verifyError(@() mdemat.tools.environment.projectRoot(installed), ...
                'mdemat:environment:NotCheckout');
        end

        function benchmarkKeepsStudyAPIAndPersistsASeparateRunRecord(testCase)
            root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
            original = path;
            guard = onCleanup(@() path(original)); %#ok<NASGU>
            addpath(fullfile(root, 'benchmarks'));
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            file = fullfile(fixture.Folder, 'study.json');
            study = lowRankTransportStudy(SpatialCells=6, Velocities=3, ...
                Tolerances=1e-4, Repetitions=1, OutputFile=file, CodeVersion='test-input');
            record = mdemat.tools.io.readRecord(file);
            testCase.verifyEqual(record.Solution, study);
            testCase.verifyEqual(record.Format, 'mdemat.runRecord');
            testCase.verifyEqual(record.SolveReport.Status, 'COMPLETED');
            testCase.verifyEqual(record.ValidationReport.Status, 'UNVERIFIED');
            testCase.verifyEqual(record.ValidationReport.Reference, ...
                'Independent dense semidiscrete matrix exponential');
            testCase.verifyEqual(record.Provenance.CodeVersion, 'test-input');
        end
    end
end
