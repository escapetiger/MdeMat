classdef TestOrderedMap < matlab.unittest.TestCase
    % TESTORDEREDMAP Unit tests for OrderedMap class.
    %
    %   TestOrderedMap provides test coverage for OrderedMap functionality
    %   including constructor modes, key-value storage, order preservation,
    %   and manipulation methods.
    %
    % Examples:
    %   % Run all tests
    %   results = runtests(core.data.tests.TestOrderedMap);
    %
    % See Also:
    %   core.data.OrderedMap

    properties
        emptyMap     % Empty OrderedMap for testing
        populatedMap % Pre-populated OrderedMap for testing
        testKeys     % Standard test keys
        testValues   % Standard test values
    end

    properties (TestParameter)
        dataType = struct(...
            'numeric', [1, 2, 3], ...
            'string', 'test string', ...
            'cell', {{'a', 'b', 'c'}}, ...
            'logical', [true, false, true], ...
            'complex', [1+2i, 3+4i], ...
            'struct', struct('field1', 10, 'field2', 'value'));
    end

    properties (Access = private, Constant)
        EXPECTED_ERROR_PREFIX = 'core:data:OrderedMap:';
    end

    methods(TestMethodSetup)
        function setupMaps(testCase)
            % SETUPMAPS Initialize test fixtures for each test method.
            
            testCase.emptyMap = core.data.OrderedMap();

            testCase.testKeys = {'first', 'second', 'third'};
            testCase.testValues = {42, 'test value', [1, 2, 3]};

            testCase.populatedMap = core.data.OrderedMap(...
                testCase.testKeys, testCase.testValues);
        end
    end

    methods(Test)
        function testEmptyConstructor(testCase)
            % TESTEMPTYCONSTRUCTOR Test empty constructor creates empty map.

            map = core.data.OrderedMap();
            testCase.verifyEqual(length(map), 0);
            testCase.verifyEqual(map.keys, {});
            testCase.verifyEqual(map.values, {});
        end

        function testConstructorWithKeysValues(testCase)
            % TESTCONSTRUCTORWITHKEYSVALUES Test array constructor mode.

            keys = {'a', 'b', 'c'};
            values = {1, 2, 3};
            map = core.data.OrderedMap(keys, values);

            testCase.verifyEqual(length(map), 3);
            testCase.verifyEqual(map.keys, keys);

            for iKey = 1:length(keys)
                testCase.verifyEqual(map(keys{iKey}), values{iKey});
            end
        end

        function testConstructorWithPairs(testCase)
            % TESTCONSTRUCTORWITHPAIRS Test key-value pairs constructor mode.

            map = core.data.OrderedMap('key1', 'value1', 'key2', 'value2', 'key3', 'value3');

            expectedKeys = {'key1', 'key2', 'key3'};
            expectedValues = {'value1', 'value2', 'value3'};

            testCase.verifyEqual(length(map), 3);
            testCase.verifyEqual(map.keys, expectedKeys);

            for iKey = 1:length(expectedKeys)
                testCase.verifyEqual(map(expectedKeys{iKey}), expectedValues{iKey});
            end
        end

        function testConstructorWithEmptyArrays(testCase)
            % TESTCONSTRUCTORWITHEMPTYARRAYS Test empty arrays constructor.

            map = core.data.OrderedMap({}, {});
            testCase.verifyEqual(length(map), 0);
            testCase.verifyEqual(map.keys, {});
            testCase.verifyEqual(map.values, {});
        end

        function testConstructorErrorInvalidArgumentCount(testCase)
            % TESTCONSTRUCTORERRORINVALIDARGUMENTCOUNT Test error on odd number of arguments.

            testCase.verifyError(...
                @() core.data.OrderedMap('key1', 'value1', 'key2'), ...
                [testCase.EXPECTED_ERROR_PREFIX, 'InvalidInput']);
        end

        function testConstructorErrorLengthMismatch(testCase)
            % TESTCONSTRUCTORERRORLENGTHMISMATCH Test error on mismatched array lengths.

            testCase.verifyError(...
                @() core.data.OrderedMap({'a', 'b'}, {'x', 'y', 'z'}), ...
                [testCase.EXPECTED_ERROR_PREFIX, 'LengthMismatch']);
        end

        function testManualPopulation(testCase)
            % TESTMANUALPOPULATION Test adding key-value pairs individually.

            map = core.data.OrderedMap();

            for iKey = 1:length(testCase.testKeys)
                map(testCase.testKeys{iKey}) = testCase.testValues{iKey};
            end

            testCase.verifyEqual(length(map), 3);
            testCase.verifyEqual(map.keys, testCase.testKeys);

            for iKey = 1:length(testCase.testKeys)
                testCase.verifyEqual(map(testCase.testKeys{iKey}), testCase.testValues{iKey});
            end
        end

        function testLength(testCase)
            % TESTLENGTH Test length method returns correct count.

            testCase.verifyEqual(length(testCase.emptyMap), 0);
            testCase.verifyEqual(length(testCase.populatedMap), 3);
        end

        function testCount(testCase)
            % TESTCOUNT Test count method compatibility.

            testCase.verifyEqual(testCase.emptyMap.count(), 0);
            testCase.verifyEqual(testCase.populatedMap.count(), 3);
        end

        function testSize(testCase)
            % TESTSIZE Test size method with and without dimension parameter.

            testCase.verifyEqual(size(testCase.emptyMap), [1 0]);
            testCase.verifyEqual(size(testCase.populatedMap), [1 3]);

            testCase.verifyEqual(size(testCase.populatedMap, 1), 1);
            testCase.verifyEqual(size(testCase.populatedMap, 2), 3);
            testCase.verifyEqual(size(testCase.populatedMap, 3), 1);
        end

        function testDataStorage(testCase)
            % TESTDATASTORAGE Test storing and overwriting values.

            map = testCase.emptyMap;

            map('temperature') = [20, 22, 25];
            testCase.verifyEqual(length(map), 1);
            testCase.verifyEqual(map('temperature'), [20, 22, 25]);

            % Test overwriting existing key
            map('temperature') = [30, 32, 35];
            testCase.verifyEqual(length(map), 1);
            testCase.verifyEqual(map('temperature'), [30, 32, 35]);

            map('pressure') = [1013, 1014, 1015];
            testCase.verifyEqual(length(map), 2);
            testCase.verifyEqual(map('pressure'), [1013, 1014, 1015]);
        end

        function testDataRetrieval(testCase)
            % TESTDATARETRIEVAL Test retrieving values and error on missing keys.

            map = testCase.populatedMap;

            for iKey = 1:length(testCase.testKeys)
                key = testCase.testKeys{iKey};
                expectedValue = testCase.testValues{iKey};
                testCase.verifyEqual(map(key), expectedValue);
            end

            % Test key not found
            testCase.verifyError(@() map('nonexistent'), ...
                [testCase.EXPECTED_ERROR_PREFIX, 'KeyNotFound']);
        end

        function testMultipleDataTypes(testCase, dataType)
            % TESTMULTIPLEDATATYPES Test storage and retrieval of various data types.

            map = testCase.emptyMap;
            map('testVar') = dataType;
            retrievedValue = map('testVar');
            testCase.verifyEqual(retrievedValue, dataType);
        end

        function testIsKey(testCase)
            % TESTISKEY Test key existence checking.

            map = testCase.populatedMap;

            for iKey = 1:length(testCase.testKeys)
                testCase.verifyTrue(map.isKey(testCase.testKeys{iKey}));
            end

            testCase.verifyFalse(map.isKey('nonexistent'));
        end

        function testRemove(testCase)
            % TESTREMOVE Test key removal and error handling.

            map = testCase.populatedMap;
            initialLength = length(map);

            keyToRemove = testCase.testKeys{1};
            map.remove(keyToRemove);

            testCase.verifyEqual(length(map), initialLength - 1);
            testCase.verifyFalse(map.isKey(keyToRemove));

            % Test removing non-existent key
            testCase.verifyError(@() map.remove('nonexistent'), ...
                [testCase.EXPECTED_ERROR_PREFIX, 'KeyNotFound']);
        end

        function testKeysMethod(testCase)
            % TESTKEYSMETHOD Test keys property returns ordered keys.

            map = testCase.populatedMap;
            keySet = map.keys;
            testCase.verifyEqual(keySet, testCase.testKeys);
        end

        function testValuesMethod(testCase)
            % TESTVALUESMETHOD Test values property returns ordered values.

            map = testCase.populatedMap;
            valueSet = map.values;

            testCase.verifyEqual(length(valueSet), length(testCase.testValues));

            for iKey = 1:length(testCase.testKeys)
                testCase.verifyEqual(valueSet{iKey}, testCase.testValues{iKey});
            end
        end

        function testKeysOrdering(testCase)
            % TESTKEYSORDERING Test keys maintain insertion order.

            map = core.data.OrderedMap();
            keys = {'d', 'a', 'c', 'b'};

            for iKey = 1:length(keys)
                map(keys{iKey}) = iKey;
            end

            returnedKeys = map.keys;
            testCase.verifyEqual(returnedKeys, keys);
        end

        function testItems(testCase)
            % TESTITEMS Test items property returns key-value pairs matrix.

            map = testCase.populatedMap;
            items = map.items;

            testCase.verifyEqual(size(items), [3, 2]);

            for iItem = 1:3
                testCase.verifyEqual(items{iItem, 1}, testCase.testKeys{iItem});
                testCase.verifyEqual(items{iItem, 2}, testCase.testValues{iItem});
            end
        end

        function testClear(testCase)
            % TESTCLEAR Test clearing all key-value pairs.

            map = testCase.populatedMap;
            map.clear();

            testCase.verifyEqual(length(map), 0);
            testCase.verifyEqual(map.keys, {});
            testCase.verifyEqual(map.values, {});
        end

        function testMove(testCase)
            % TESTMOVE Test moving keys to different positions.

            map = testCase.populatedMap;

            % Move 'first' to position 2
            map.move('first', 2);
            expectedKeys = {'second', 'first', 'third'};
            testCase.verifyEqual(map.keys, expectedKeys);

            % Move 'third' to position 1
            map.move('third', 1);
            expectedKeys = {'third', 'second', 'first'};
            testCase.verifyEqual(map.keys, expectedKeys);

            % Test invalid position
            testCase.verifyError(@() map.move('first', 4), ...
                [testCase.EXPECTED_ERROR_PREFIX, 'InvalidPosition']);

            % Test invalid key
            testCase.verifyError(@() map.move('nonexistent', 1), ...
                [testCase.EXPECTED_ERROR_PREFIX, 'KeyNotFound']);
        end

        function testMoveToFront(testCase)
            % TESTMOVETOFRONT Test moving key to beginning of order.

            map = testCase.populatedMap;

            map.moveToFront('third');
            expectedKeys = {'third', 'first', 'second'};
            testCase.verifyEqual(map.keys, expectedKeys);

            % Test with non-existent key
            map.moveToFront('nonexistent');
            testCase.verifyEqual(map.keys, expectedKeys);
        end

        function testMoveToEnd(testCase)
            % TESTMOVETOEND Test moving key to end of order.

            map = testCase.populatedMap;

            map.moveToEnd('first');
            expectedKeys = {'second', 'third', 'first'};
            testCase.verifyEqual(map.keys, expectedKeys);

            % Test with non-existent key
            map.moveToEnd('nonexistent');
            testCase.verifyEqual(map.keys, expectedKeys);
        end

        function testKeyType(testCase)
            % TESTKEYTYPE Test keyType property returns correct type.

            map1 = core.data.OrderedMap();
            testCase.verifyEqual(map1.keyType, 'char');

            map2 = core.data.OrderedMap({1}, {'k'});
            testCase.verifyEqual(map2.keyType, 'double');
        end

        function testValueType(testCase)
            % TESTVALUETYPE Test valueType property returns correct type.

            map1 = core.data.OrderedMap();
            testCase.verifyEqual(map1.valueType, 'any');

            map2 = core.data.OrderedMap({'char'}, {'char'});
            testCase.verifyEqual(map2.valueType, 'char');
        end

        function testToTable(testCase)
            % TESTTOTABLE Test conversion to MATLAB table.

            map = testCase.populatedMap;
            tbl = map.toTable();

            testCase.verifyEqual(size(tbl), [3, 2]);
            testCase.verifyEqual(tbl.Properties.VariableNames, {'Key', 'Value'});

            testCase.verifyEqual(tbl.Key, testCase.testKeys(:));

            for iValue = 1:length(testCase.testValues)
                testCase.verifyEqual(tbl.Value{iValue}, testCase.testValues{iValue});
            end

            % Test empty map
            emptyTbl = testCase.emptyMap.toTable();
            testCase.verifyEqual(size(emptyTbl), [0, 2]);
            testCase.verifyEqual(emptyTbl.Properties.VariableNames, {'Key', 'Value'});
        end

        function testDisplay(testCase)
            % TESTDISPLAY Test custom display output.

            map = testCase.populatedMap;

            output = evalc('disp(map)');

            testCase.verifyTrue(contains(output, 'OrderedMap'));
            testCase.verifyTrue(contains(output, '3 key/value pairs'));

            for iKey = 1:length(testCase.testKeys)
                testCase.verifyTrue(contains(output, testCase.testKeys{iKey}));
            end

            % Test empty map display
            emptyOutput = evalc('disp(testCase.emptyMap)');
            testCase.verifyTrue(contains(emptyOutput, '0 key/value pairs'));
        end

        function testDuplicateKeysInConstructor(testCase)
            % TESTDUPLICATEKEYSINCONSTRUCTOR Test behavior with duplicate keys.

            keys = {'a', 'b', 'a'};
            values = {1, 2, 3};
            map = core.data.OrderedMap(keys, values);

            testCase.verifyEqual(map('a'), 3);
            testCase.verifyEqual(length(map), 2);

            % Test with pairs constructor
            map2 = core.data.OrderedMap('key1', 'value1', 'key1', 'value2');
            testCase.verifyEqual(map2('key1'), 'value2');
            testCase.verifyEqual(length(map2), 1);
        end

        function testNumericKeys(testCase)
            % TESTNUMERICKEYS Test using numeric keys.

            keys = {1, 2, 3};
            values = {'one', 'two', 'three'};
            map = core.data.OrderedMap(keys, values);

            testCase.verifyEqual(length(map), 3);
            testCase.verifyEqual(map(1), 'one');
            testCase.verifyEqual(map(2), 'two');
            testCase.verifyEqual(map(3), 'three');
        end

        function testOverwriteExistingKey(testCase)
            % TESTOVERWRITEEXISTINGKEY Test overwriting values preserves order.

            map = testCase.populatedMap;
            originalKeys = map.keys;

            map('second') = 'new value';

            testCase.verifyEqual(map.keys, originalKeys);
            testCase.verifyEqual(map('second'), 'new value');
            testCase.verifyEqual(length(map), 3);
        end
    end
end