classdef TestRecordIO < matlab.unittest.TestCase
    properties
        Folder
    end
    methods (TestMethodSetup)
        function temporaryOutput(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.Folder = fixture.Folder;
        end
    end
    methods (Test)
        function fixedV1ContainersRemainReadableWithoutProducerMetadata(testCase)
            % Literal v1 fixtures are independent of the current writeRecord.
            file = fullfile(testCase.Folder,'v1.json');
            descriptor = ['{"format":"mdemat.numericRecord","schemaVersion":1,', ...
                '"tree":{"kind":"array","variable":"n000001"}}'];
            literal = [descriptor(1:end-1),',"arrays":{"n000001":', ...
                '{"type":"double","shape":[1,2],"sparse":false,"complex":false,', ...
                '"real":"AAAAAAAA8D8AAAAAAAAAQA=="}}}'];
            writeText(file,literal);
            [value,info] = mdemat.tools.io.readRecord(file);
            testCase.verifyEqual(value,[1,2]);
            testCase.verifyEqual(info.ToolboxVersion,'');
            manifest = unicode2native(descriptor,'UTF-8'); %#ok<NASGU>
            n000001 = [1,2]; %#ok<NASGU>
            mat = fullfile(testCase.Folder,'v1.mat');
            save(mat,'manifest','n000001','-v7');
            testCase.verifyEqual(mdemat.tools.io.readRecord(mat),[1,2]);
            mask = ['{"format":"mdemat.numericRecord","schemaVersion":1,', ...
                '"tree":{"kind":"string","shape":[1,2],"text":["alpha",""],', ...
                '"missing":[0,1]},"arrays":{}}'];
            writeText(file,mask);
            testCase.verifyEqual(mdemat.tools.io.readRecord(file),["alpha",missing]);
        end

        function matAndJsonPreserveDataTypesShapesAndValues(testCase)
            data = struct();
            data.coefficients = reshape(single(1:24), [2, 3, 4]);
            data.integer = bitshift(uint64(1), 60) + uint64(3);
            data.signed = int8([-128, 0, 127]);
            data.special = [NaN, Inf, -Inf, -0.0];
            data.complex = [1 + 2i, 0; -3i, 4];
            data.sparse = sparse([1, 0, 3i; 0, -2, 0]);
            data.mask = sparse(logical([1, 0; 0, 1]));
            data.empty = {zeros(0, 3, 'single'), cell(2, 0), repmat('x', 1, 0), ''};
            data.text = ["alpha", ""; "β", missing];
            data.emptyText = strings(0, 3);
            data.otherEmptyText = strings(2,0);
            data.scalarEmptyText = "";
            data.scalarMissingText = string(missing);
            data.characterMatrix = ['αβ';'γδ'];
            data.zeroImaginary = complex([0,1],[0,0]);
            data.emptyStruct = repmat(struct(),0,3);
            data.fieldlessArray = repmat(struct(),2,3);
            data.otherPrimitives = {uint8([0,255]),int16([-32768,32767]), ...
                uint16([0,65535]),int32([intmin('int32'),intmax('int32')]), ...
                uint32([0,intmax('uint32')]),int64([intmin('int64'),intmax('int64')]), ...
                single([1+2i,-3i]),logical([true,false])};
            data.records = [struct('x', 2, 'label', 'a'), struct('x', 3, 'label', 'b')];
            data.noFields = struct();
            for extension = {'.mat', '.json'}
                file = fullfile(testCase.Folder, ['record', extension{1}]);
                mdemat.tools.io.writeRecord(file, data);
                [restored, info] = mdemat.tools.io.readRecord(file);
                testCase.verifyEqual(restored, data);
                testCase.verifyEqual(info.SchemaVersion, 1);
                testCase.verifyTrue(issparse(restored.sparse));
                testCase.verifyTrue(issparse(restored.mask));
                testCase.verifyClass(restored.integer, 'uint64');
                testCase.verifySize(restored.empty{3}, [1, 0]);
                testCase.verifyEqual(isreal(restored.zeroImaginary),isreal(data.zeroImaginary));
            end
        end

        function refusesOverwriteAndRequiresExplicitDestinationCreation(testCase)
            file = fullfile(testCase.Folder, 'nested', 'record.json');
            testCase.verifyError(@() mdemat.tools.io.writeRecord(file, 1), ...
                'mdemat:io:DirectoryNotFound');
            testCase.verifyFalse(isfolder(fileparts(file)));
            mdemat.tools.io.writeRecord(file, 1, createDirectory=true);
            testCase.verifyError(@() mdemat.tools.io.writeRecord(file, 2), ...
                'mdemat:io:FileExists');
            testCase.verifyEqual(mdemat.tools.io.readRecord(file), 1);
            mdemat.tools.io.writeRecord(file, 2, overwrite=true);
            testCase.verifyEqual(mdemat.tools.io.readRecord(file), 2);
        end

        function rejectsOpaqueValuesBeforeCreatingOutput(testCase)
            file = fullfile(testCase.Folder, 'invalid.mat');
            testCase.verifyError(@() mdemat.tools.io.writeRecord(file, {@sin}), ...
                'mdemat:io:UnsupportedValue');
            testCase.verifyFalse(isfile(file));
            testCase.verifyError(@() mdemat.tools.io.writeRecord(file, containers.Map), ...
                'mdemat:io:UnsupportedValue');
            testCase.verifyFalse(isfile(file));
        end

        function rejectsArbitraryMatSerialization(testCase)
            file = fullfile(testCase.Folder, 'opaque.mat');
            opaque = @sin; %#ok<NASGU>
            save(file, 'opaque');
            testCase.verifyError(@() mdemat.tools.io.readRecord(file), ...
                'mdemat:io:InvalidRecord');
        end

        function reportsInvalidSchemaAndTruncatedData(testCase)
            file = fullfile(testCase.Folder, 'record.json');
            mdemat.tools.io.writeRecord(file, [1, 2]);
            envelope = jsondecode(fileread(file));
            envelope.schemaVersion = 99;
            writeText(file, jsonencode(envelope));
            testCase.verifyError(@() mdemat.tools.io.readRecord(file), ...
                'mdemat:io:UnsupportedSchema');
            envelope.schemaVersion = 1;
            envelope.arrays.n000001.real = '';
            writeText(file, jsonencode(envelope));
            testCase.verifyError(@() mdemat.tools.io.readRecord(file), ...
                'mdemat:io:InvalidRecord');
            testCase.verifyNumElements(dir(fullfile(testCase.Folder, '*.json')), 1);
        end

        function rejectsBroadcastImaginaryAndMalformedFlagsBeforeReconstruction(testCase)
            file = fullfile(testCase.Folder,'complex.json');
            mdemat.tools.io.writeRecord(file,[1+2i,3+4i]);
            original = jsondecode(fileread(file));
            envelope = original;
            envelope.arrays.n000001.imaginary = matlab.net.base64encode(typecast(9,'uint8'));
            writeText(file,jsonencode(envelope));
            testCase.verifyError(@() mdemat.tools.io.readRecord(file),'mdemat:io:InvalidRecord');
            envelope = original;
            envelope.arrays.n000001.complex = [true,false];
            writeText(file,jsonencode(envelope));
            testCase.verifyError(@() mdemat.tools.io.readRecord(file),'mdemat:io:InvalidRecord');
        end

        function resourceBudgetsRejectDescriptionsWithoutLargeAllocations(testCase)
            file = fullfile(testCase.Folder,'shape.json');
            mdemat.tools.io.writeRecord(file,struct());
            envelope = jsondecode(fileread(file));
            envelope.tree.shape = [1e6,1e6];
            writeText(file,jsonencode(envelope));
            testCase.verifyError(@() mdemat.tools.io.readRecord(file), ...
                'mdemat:io:ResourceLimit');
            mdemat.tools.io.writeRecord(file,ones(1,8),overwrite=true);
            testCase.verifyError(@() mdemat.tools.io.readRecord(file,limits=struct('MaxArrayBytes',8)), ...
                'mdemat:io:ResourceLimit');
            testCase.verifyError(@() mdemat.tools.io.readRecord(file,limits=struct('MaxFileBytes',8)), ...
                'mdemat:io:ResourceLimit');
            testCase.verifyError(@() mdemat.tools.io.writeRecord(file,ones(1,9),overwrite=true, ...
                limits=struct('MaxArrayBytes',8)),'mdemat:io:ResourceLimit');
            testCase.verifyEqual(mdemat.tools.io.readRecord(file),ones(1,8));
            testCase.verifyNumElements(dir(fullfile(testCase.Folder,'*')),3);
        end

        function sparseDuplicatesAndMissingMasksAreNotSilentlyCoerced(testCase)
            file = fullfile(testCase.Folder,'sparse.json');
            mdemat.tools.io.writeRecord(file,sparse([1,2],[1,1],[3,4],2,2));
            envelope = jsondecode(fileread(file));
            envelope.arrays.n000001.rows = [1,1];
            writeText(file,jsonencode(envelope));
            testCase.verifyError(@() mdemat.tools.io.readRecord(file),'mdemat:io:InvalidRecord');
            mdemat.tools.io.writeRecord(file,["a",missing],overwrite=true);
            envelope = jsondecode(fileread(file));
            envelope.tree.missing = [false,2];
            writeText(file,jsonencode(envelope));
            testCase.verifyError(@() mdemat.tools.io.readRecord(file),'mdemat:io:InvalidRecord');
        end

        function invalidDescriptionsAndLimitPoliciesFailAtTheBoundary(testCase)
            file = fullfile(testCase.Folder,'invalid.json');
            mdemat.tools.io.writeRecord(file,logical([1,0]));
            original = jsondecode(fileread(file));
            bad = original;
            bad.arrays.n000001.real = matlab.net.base64encode(uint8([2,0]));
            writeText(file,jsonencode(bad));
            testCase.verifyError(@() mdemat.tools.io.readRecord(file),'mdemat:io:InvalidRecord');
            bad = original;
            bad.tree.variable = 'n999999';
            writeText(file,jsonencode(bad));
            testCase.verifyError(@() mdemat.tools.io.readRecord(file),'mdemat:io:InvalidRecord');
            bad = original;
            bad.tree.extra = 1;
            writeText(file,jsonencode(bad));
            testCase.verifyError(@() mdemat.tools.io.readRecord(file),'mdemat:io:InvalidRecord');
            bad = original;
            bad.arrays.n000001.real = '????';
            writeText(file,jsonencode(bad));
            testCase.verifyError(@() mdemat.tools.io.readRecord(file),'mdemat:io:InvalidRecord');
            testCase.verifyError(@() mdemat.tools.io.readRecord(file,limits=struct('Unknown',1)), ...
                'mdemat:io:InvalidLimits');
            testCase.verifyError(@() mdemat.tools.io.readRecord(file,limits=struct('MaxDepth',0)), ...
                'mdemat:io:InvalidLimits');
            testCase.verifyError(@() mdemat.tools.io.writeRecord(file,{{1}},overwrite=true, ...
                limits=struct('MaxDepth',2)),'mdemat:io:ResourceLimit');
            testCase.verifyError(@() mdemat.tools.io.writeRecord(file,struct('a',1,'b',2), ...
                overwrite=true,limits=struct('MaxTotalBytes',8)),'mdemat:io:ResourceLimit');
        end

        function failedCommitRetainsTheExistingDestinationAndRemovesTemporaryFile(testCase)
            destination = fullfile(testCase.Folder,'occupied.mat');
            mkdir(destination);
            testCase.verifyError(@() mdemat.tools.io.writeRecord(destination,1), ...
                'mdemat:io:WriteFailed');
            testCase.verifyTrue(isfolder(destination));
            listing = dir(testCase.Folder);
            testCase.verifyEmpty(listing(~[listing.isdir]));
        end
    end
end

function writeText(file, text)
stream = fopen(file, 'w', 'n', 'UTF-8');
guard = onCleanup(@() fclose(stream)); %#ok<NASGU>
fwrite(stream, unicode2native(text, 'UTF-8'), 'uint8');
end
