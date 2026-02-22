classdef TestL1MultiIndexer < matlab.unittest.TestCase
    % TESTL1MULTIINDEXER Unit tests for L1MultiIndexer class.
    %
    %   TestL1MultiIndexer provides test coverage for L1MultiIndexer
    %   functionality including index generation, caching, and conversion
    %   between multi-indices and linear indices.
    %
    % Examples:
    %   % Run all tests
    %   results = runtests(core.linalg.tests.TestL1MultiIndexer);
    %
    % See Also:
    %   core.linalg.L1MultiIndexer

    properties (TestParameter)
        style = {'F', 'C'}     % Storage ordering styles
        dimension = {1, 2, 3, 4} % Test dimensions
    end

    methods (Test)
        function test1DIndices(testCase, style)
            % TEST1DINDICES Test 1D index generation.

            indexer = core.linalg.L1MultiIndexer(style);

            indexer.setCache(1, 3);
            M = indexer.cache;

            expected = (1:3)';
            testCase.verifyEqual(M, expected);
        end

        function test2DIndicesUpper2(testCase)
            % TEST2DINDICESUPPER2 Test 2D indices with upper bound 2.

            indexer = core.linalg.L1MultiIndexer('F');

            indexer.setCache(2, 2);
            M = indexer.cache;

            expected = [1, 1];
            testCase.verifyEqual(M, expected);
        end

        function test2DIndicesUpper3F(testCase)
            % TEST2DINDICESUPPER3F Test 2D indices with F-style ordering.

            indexer = core.linalg.L1MultiIndexer('F');

            indexer.setCache(2, 3);
            M = indexer.cache;

            expected = [; ...
                1, 1; ...
                1, 2; ...
                2, 1; ...
                ];

            M = sortrows(M);
            expected = sortrows(expected);
            testCase.verifyEqual(M, expected);
        end

        function test2DIndicesUpper3C(testCase)
            % TEST2DINDICESUPPER3C Test 2D indices with C-style ordering.

            indexer = core.linalg.L1MultiIndexer('C');

            indexer.setCache(2, 3);
            M = indexer.cache;

            expected = [; ...
                1, 1; ...
                2, 1; ...
                1, 2; ...
                ];

            M = sortrows(M);
            expected = sortrows(expected);
            testCase.verifyEqual(M, expected);
        end

        function test3DIndices(testCase)
            % TEST3DINDICES Test 3D index generation.

            indexer = core.linalg.L1MultiIndexer('F');

            indexer.setCache(3, 3);
            M = indexer.cache;

            expected = [1, 1, 1];
            testCase.verifyEqual(M, expected);
        end

        function test3DIndicesUpper4(testCase)
            % TEST3DINDICESUPPER4 Test 3D indices with upper bound 4.

            indexer = core.linalg.L1MultiIndexer('F');

            indexer.setCache(3, 4);
            M = indexer.cache;

            expected = [; ...
                1, 1, 1; ...
                1, 1, 2; ...
                1, 2, 1; ...
                2, 1, 1; ...
                ];

            M = sortrows(M);
            expected = sortrows(expected);
            testCase.verifyEqual(M, expected);
        end

        function testInvalidInput(testCase)
            % TESTINVALIDINPUT Test error handling for invalid inputs.

            indexer = core.linalg.L1MultiIndexer();

            testCase.verifyError(@() indexer.setCache(2, 1), ...
                'core:linalg:L1MultiIndexer:InvalidInput');
            testCase.verifyError(@() indexer.setCache(2, -1), ...
                'core:linalg:L1MultiIndexer:InvalidInput');
            testCase.verifyError(@() indexer.setCache(), ...
                'core:linalg:L1MultiIndexer:InvalidInput');
        end

        function testLargeUpper(testCase)
            % TESTLARGEUPPER Test generation with large upper bound.

            indexer = core.linalg.L1MultiIndexer('F');

            indexer.setCache(2, 10);
            M = indexer.cache;

            testCase.verifyTrue(all(M(:) >= 1));
            testCase.verifyEqual(size(unique(M, 'rows'), 1), size(M, 1));

            rowSums = sum(M, 2);
            testCase.verifyTrue(all(rowSums <= 10));
            testCase.verifyTrue(all(rowSums >= 2));
        end

        function testDifferentDimensions(testCase, dimension)
            % TESTDIFFERENTDIMENSIONS Test generation with various dimensions.

            if dimension < 5
                indexer = core.linalg.L1MultiIndexer('F');

                upper = dimension + 1;
                indexer.setCache(dimension, upper);
                M = indexer.cache;

                testCase.verifyEqual(size(M, 2), dimension);

                rowSums = sum(M, 2);
                testCase.verifyTrue(all(rowSums <= upper));
                testCase.verifyTrue(all(rowSums >= dimension));

                testCase.verifyTrue(all(M(:) >= 1));
                testCase.verifyTrue(all(floor(M(:)) == M(:)));
            end
        end

        function testResultConsistency(testCase)
            % TESTRESULTCONSISTENCY Test consistent results across calls.

            indexer = core.linalg.L1MultiIndexer('F');

            indexer.setCache(2, 4);
            M1 = indexer.cache;

            indexer.setCache(2, 4);
            M2 = indexer.cache;

            testCase.verifyEqual(M1, M2);
        end

        function testStorageOrderDifference(testCase)
            % TESTSTORAGEORDERDIFFERENCE Test F vs C style ordering.

            indexerF = core.linalg.L1MultiIndexer('F');
            indexerC = core.linalg.L1MultiIndexer('C');

            indexerF.setCache(3, 4);
            indexerC.setCache(3, 4);

            MF = indexerF.cache;
            MC = indexerC.cache;

            testCase.verifyEqual(size(MF), size(MC));

            MF_sorted = sortrows(MF);
            MC_sorted = sortrows(MC);

            testCase.verifyEqual(MF_sorted, MC_sorted);
        end

        function testMinimalDimensionSum(testCase)
            % TESTMINIMALDIMENSIONSUM Test minimal sum equals dimension.

            for d = 1:4
                indexer = core.linalg.L1MultiIndexer('F');
                indexer.setCache(d, d);
                M = indexer.cache;

                expected = ones(1, d);

                testCase.verifyEqual(M, expected);
                testCase.verifyEqual(sum(M), d);
            end
        end

        function testLowerUpperBounds(testCase)
            % TESTLOWERUPPERBOUNDS Test with both lower and upper bounds.

            indexer = core.linalg.L1MultiIndexer('F');

            indexer.setCache(2, 3, 4);
            M = indexer.cache;

            expected = [; ...
                1, 2; ...
                2, 1; ...
                1, 3; ...
                2, 2; ...
                3, 1; ...
                ];

            M = sortrows(M);
            expected = sortrows(expected);

            testCase.verifyEqual(M, expected);
        end

        function testMultiToLinearF(testCase)
            % TESTMULTITOLINEARF Test multi-index to linear conversion.

            indexer = core.linalg.L1MultiIndexer('F');

            indexer.setCache(3, 4);

            testIndices = [; ...
                1, 1, 1; ...
                2, 1, 1; ...
                1, 2, 1; ...
                1, 1, 2; ...
                ];

            linearIndices = indexer.multiToLinear(testIndices);
            expected = (1:size(testIndices, 1))';

            testCase.verifyEqual(linearIndices, expected);
        end

        function testLinearToMulti(testCase)
            % TESTLINEARTOMULTI Test linear to multi-index conversion.

            indexer = core.linalg.L1MultiIndexer('F');

            indexer.setCache(3, 4);
            originalCache = indexer.cache;

            linearIndices = (1:size(originalCache, 1))';

            multiIndices = indexer.linearToMulti(linearIndices);

            testCase.verifyEqual(multiIndices, originalCache);
        end

        function testRoundTripConversion(testCase)
            % TESTROUNDTRIPCONVERSION Test round-trip index conversion.

            indexer = core.linalg.L1MultiIndexer('F');

            indexer.setCache(3, 5);
            original = indexer.cache;

            linear = indexer.multiToLinear(original);
            reconstructed = indexer.linearToMulti(linear);

            testCase.verifyEqual(reconstructed, original);
        end

        function testInvalidLinearIndex(testCase)
            % TESTINVALIDLINEARINDEX Test error handling for invalid linear indices.

            indexer = core.linalg.L1MultiIndexer('F');

            indexer.setCache(2, 3);

            testCase.verifyError(@() indexer.linearToMulti(0), ...
                'core:linalg:CachedMultiIndexer:InvalidInput');
            testCase.verifyError(@() indexer.linearToMulti(-1), ...
                'core:linalg:CachedMultiIndexer:InvalidInput');
            testCase.verifyError(@() indexer.linearToMulti(10000), ...
                'core:linalg:CachedMultiIndexer:InvalidInput');
        end

        function testNIndices(testCase)
            % TESTNINDICES Test nIndices property.

            indexer = core.linalg.L1MultiIndexer('F');

            indexer.setCache(2, 5);

            testCase.verifyEqual(indexer.nIndices, size(indexer.cache, 1));
        end
    end
end