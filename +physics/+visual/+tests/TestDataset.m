classdef TestDataset < matlab.unittest.TestCase
    % TESTDATASET Unit tests for the Dataset class.
    %
    %   TestDataset provides comprehensive test coverage for the Dataset 
    %   class functionality including constructor validation, data 
    %   manipulation methods, and type handling.
    %
    % See also:
    %   physics.visual.Dataset

    properties (Constant)
        Tolerance = 1e-12
    end

    methods (Test)
        function testEmptyConstructor(testCase)
            dataset = physics.visual.Dataset();
            testCase.verifyEqual(dataset.Type, 'numeric');
            testCase.verifyClass(dataset.Data, 'struct');
            testCase.verifyEmpty(fieldnames(dataset.Data));
        end

        function testConstructorWithType(testCase)
            dataset = physics.visual.Dataset('type', 'fixed');
            testCase.verifyEqual(dataset.Type, 'fixed');
            testCase.verifyClass(dataset.Data, 'struct');
            
            dataset = physics.visual.Dataset('type', 'exact');
            testCase.verifyEqual(dataset.Type, 'exact');
            
            dataset = physics.visual.Dataset('type', 'numeric');
            testCase.verifyEqual(dataset.Type, 'numeric');
        end

        function testConstructorWithTypeAndData(testCase)
            data = struct('u', [1, 2, 3], 'v', [4, 5, 6]);
            dataset = physics.visual.Dataset('type', 'fixed', 'data', data);
            testCase.verifyEqual(dataset.Type, 'fixed');
            testCase.verifyEqual(dataset.Data, data);
        end

        function testInvalidTypeValidation(testCase)
            testCase.verifyError(@() physics.visual.Dataset('type', 'invalid'), ...
                'MATLAB:validators:mustBeMember');
            testCase.verifyError(@() physics.visual.Dataset('type', 123), ...
                'MATLAB:ISMEMBER:InputClass');
        end
        
        function testSetAndGetData(testCase)
            dataset = physics.visual.Dataset('type', 'fixed');
            
            testData = [1, 2, 3, 4];
            dataset = dataset.setData('velocity', testData);
            
            retrievedData = dataset.getData('velocity');
            testCase.verifyEqual(retrievedData, testData);
            
            emptyData = dataset.getData('pressure');
            testCase.verifyEmpty(emptyData);
        end
        
        function testHasField(testCase)
            dataset = physics.visual.Dataset('type', 'exact');
            dataset = dataset.setData('temperature', [10, 20, 30]);
            
            testCase.verifyTrue(dataset.hasField('temperature'));
            testCase.verifyFalse(dataset.hasField('nonexistent'));
        end
        
        function testGetFieldNames(testCase)
            dataset = physics.visual.Dataset('type', 'numeric');
            dataset = dataset.setData('u', [1, 2, 3]);
            dataset = dataset.setData('v', [4, 5, 6]);
            
            fieldNames = dataset.getFieldNames();
            testCase.verifyLength(fieldNames, 2);
            testCase.verifyTrue(ismember('u', fieldNames));
            testCase.verifyTrue(ismember('v', fieldNames));
        end
        
        function testClearData(testCase)
            dataset = physics.visual.Dataset('type', 'fixed');
            dataset = dataset.setData('field1', [1, 2, 3]);
            dataset = dataset.setData('field2', [4, 5, 6]);
            
            testCase.verifyTrue(dataset.hasField('field1'));
            testCase.verifyTrue(dataset.hasField('field2'));
            
            dataset = dataset.clearData();
            
            testCase.verifyFalse(dataset.hasField('field1'));
            testCase.verifyFalse(dataset.hasField('field2'));
            testCase.verifyEqual(dataset.Type, 'fixed');
        end
        
        function testDataValidation(testCase)
            dataset = physics.visual.Dataset();
            
            testCase.verifyError(@() dataset.setData(123, [1, 2, 3]), ...
                'MATLAB:validators:mustBeTextScalar');
            testCase.verifyError(@() dataset.setData('field', 'invalid'), ...
                'MATLAB:validators:mustBeNumeric');
        end
        
        function testComplexWorkflow(testCase)
            data = struct('pressure', [1, 4, 9], 'velocity', [2, 3, 5]);
            dataset = physics.visual.Dataset('type', 'exact', 'data', data);
            
            testCase.verifyEqual(dataset.Type, 'exact');
            testCase.verifyLength(dataset.getFieldNames(), 2);
            testCase.verifyEqual(dataset.getData('pressure'), [1, 4, 9]);
            
            dataset = dataset.setData('density', [7, 8, 9]);
            testCase.verifyLength(dataset.getFieldNames(), 3);
            testCase.verifyTrue(dataset.hasField('density'));
            
            dataset = dataset.clearData();
            testCase.verifyEmpty(dataset.getFieldNames());
            testCase.verifyEqual(dataset.Type, 'exact');
        end
    end
end