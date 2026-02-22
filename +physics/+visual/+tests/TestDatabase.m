classdef TestDatabase < matlab.unittest.TestCase
    % TESTDATABASE Unit tests for the simplified Database class.
    %
    %   TestDatabase provides comprehensive test coverage for the simplified 
    %   Database class functionality including dataset management and 
    %   groupBy processing.
    %
    % See also:
    %   physics.visual.Database

    properties (Constant)
        Tolerance = 1e-12
    end

    methods (Test)
        function testConstructor(testCase)
            db = physics.visual.Database();
            testCase.verifyClass(db, 'physics.visual.Database');
            testCase.verifyClass(db.Datasets, 'struct');
            testCase.verifyEmpty(fieldnames(db.Datasets));
        end

        function testSetAndGetDataset(testCase)
            db = physics.visual.Database();
            dataset = physics.visual.Dataset('type', 'numeric');
            
            db = db.setDataset('test', dataset);
            testCase.verifyTrue(isfield(db.Datasets, 'test'));
            
            retrieved = db.getDataset('test');
            testCase.verifyEqual(retrieved, dataset);
            
            nonExistent = db.getDataset('nonexistent');
            testCase.verifyEmpty(nonExistent);
        end

        function testGroupByBasic(testCase)
            db = physics.visual.Database();
            
            testData = struct('u', [1, 2; 3, 4; 5, 6], 'x', {[0, 1], [0, 1]});
            dataset = physics.visual.Dataset('type', 'numeric', 'data', testData);
            db = db.setDataset('test', dataset);
            
            processedDb = db.groupBy('u1');
            testCase.verifyClass(processedDb, 'physics.visual.Database');
            
            processedDataset = processedDb.getDataset('test');
            testCase.verifyClass(processedDataset, 'physics.visual.Dataset');
            testCase.verifyEqual(processedDataset.Type, 'numeric');
            testCase.verifyTrue(isfield(processedDataset.Data, 'x'));
            testCase.verifyTrue(isfield(processedDataset.Data, 'u'));
        end

        function testGroupByMultipleDatasets(testCase)
            db = physics.visual.Database();
            
            testData1 = struct('pressure', [10, 20; 30, 40], 'x', {[0, 1], [0, 1]});
            dataset1 = physics.visual.Dataset('type', 'exact', 'data', testData1);
            db = db.setDataset('exact_data', dataset1);
            
            testData2 = struct('pressure', [15, 25; 35, 45], 'x', {[0, 1], [0, 1]});
            dataset2 = physics.visual.Dataset('type', 'numeric', 'data', testData2);
            db = db.setDataset('numeric_data', dataset2);
            
            processedDb = db.groupBy('pressure1');
            testCase.verifyLength(fieldnames(processedDb.Datasets), 2);
            
            exact = processedDb.getDataset('exact_data');
            numeric = processedDb.getDataset('numeric_data');
            
            testCase.verifyEqual(exact.Type, 'exact');
            testCase.verifyEqual(numeric.Type, 'numeric');
        end

        function testGroupByInvalidVariable(testCase)
            db = physics.visual.Database();
            testCase.verifyError(@() db.groupBy('invalidformat'), ...
                'physics:visual:Database:InvalidInput');
        end

        function testGroupByMissingField(testCase)
            db = physics.visual.Database();
            
            testData = struct('temperature', [1, 2, 3], 'x', {[0, 1, 2]});
            dataset = physics.visual.Dataset('type', 'numeric', 'data', testData);
            db = db.setDataset('test', dataset);
            
            processedDb = db.groupBy('pressure1');
            testCase.verifyEmpty(fieldnames(processedDb.Datasets));
        end

        function testGroupByComplexWorkflow(testCase)
            db = physics.visual.Database();
            
            data1 = struct('velocity', [1, 2, 3; 4, 5, 6], 'x', {[0, 0.5, 1]});
            dataset1 = physics.visual.Dataset('type', 'fixed', 'data', data1);
            db = db.setDataset('reference', dataset1);
            
            data2 = struct('velocity', [1.1, 2.1, 3.1; 4.1, 5.1, 6.1], 'x', {[0, 0.5, 1]});
            dataset2 = physics.visual.Dataset('type', 'numeric', 'data', data2);
            db = db.setDataset('simulation', dataset2);
            
            processedDb = db.groupBy('velocity2');
            
            refDataset = processedDb.getDataset('reference');
            simDataset = processedDb.getDataset('simulation');
            
            testCase.verifyEqual(refDataset.Type, 'fixed');
            testCase.verifyEqual(simDataset.Type, 'numeric');
            
            testCase.verifyEqual(refDataset.Data.u, [2; 5]);
            testCase.verifyEqual(simDataset.Data.u, [2.1; 5.1]);
            
            testCase.verifyEqual(refDataset.Data.x{1}, [0, 0.5, 1]);
            testCase.verifyEqual(simDataset.Data.x{1}, [0, 0.5, 1]);
        end

        function testNDatasets(testCase)
            % Test nDatasets property after groupBy
            db = physics.visual.Database();
            
            data1 = struct('temp', [1, 2], 'x', {[0, 1]});
            dataset1 = physics.visual.Dataset('type', 'fixed', 'data', data1);
            db = db.setDataset('test1', dataset1);
            
            data2 = struct('temp', [3, 4], 'x', {[0, 1]});
            dataset2 = physics.visual.Dataset('type', 'numeric', 'data', data2);
            db = db.setDataset('test2', dataset2);
            
            processedDb = db.groupBy('temp1');
            
            % Check if nDatasets property exists and has correct value
            if isprop(processedDb, 'nDatasets') || isfield(processedDb, 'nDatasets')
                testCase.verifyEqual(processedDb.nDatasets, 2);
            end
        end
    end
end