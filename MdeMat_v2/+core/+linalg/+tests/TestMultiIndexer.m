classdef TestMultiIndexer < matlab.unittest.TestCase
    % TESTMULTIINDEXER Unit tests for MultiIndexer class.
    %
    %   TestMultiIndexer provides test coverage for MultiIndexer
    %   functionality including index generation and conversion between
    %   multi-indices and linear indices.
    %
    % Examples:
    %   % Run all tests
    %   results = runtests(core.linalg.tests.TestMultiIndexer);
    %
    % See Also:
    %   core.linalg.MultiIndexer

    properties (TestParameter)
        style = {'F', 'C'}   % Storage ordering styles
        dimension = {1, 2, 3} % Test dimensions
    end

    methods (Test)
        function test1DGenerate(testCase, style)
            % TEST1DGENERATE Test 1D index generation.

            indexer = core.linalg.MultiIndexer(3, style);

            M = indexer.generate();

            expected = (1:3)';
            testCase.verifyEqual(M, expected);
        end

        function test2DGenerateF(testCase)
            % TEST2DGENERATEF Test 2D index generation with F-style.

            indexer = core.linalg.MultiIndexer([2, 3], 'F');

            M = indexer.generate();

            expected = [; ...
                1, 1; ...
                2, 1; ...
                1, 2; ...
                2, 2; ...
                1, 3; ...
                2, 3; ...
                ];

            testCase.verifyEqual(M, expected);
        end

        function test2DGenerateC(testCase)
            % TEST2DGENERATEC Test 2D index generation with C-style.

            indexer = core.linalg.MultiIndexer([2, 3], 'C');

            M = indexer.generate();

            expected = [; ...
                1, 1; ...
                1, 2; ...
                1, 3; ...
                2, 1; ...
                2, 2; ...
                2, 3; ...
                ];

            testCase.verifyEqual(M, expected);
        end

        function testMultiToLinear2DF(testCase)
            % TESTMULTITOLINEAR2DF Test 2D multi-index to linear conversion.

            indexer = core.linalg.MultiIndexer([2, 2], 'F');

            M = [1, 1; 2, 1; 1, 2; 2, 2];
            L = indexer.multiToLinear(M);

            expected = (1:4)';
            testCase.verifyEqual(L, expected);
        end

        function testMultiToLinear2DC(testCase)
            % TESTMULTITOLINEAR2DC Test 2D multi-index to linear conversion C-style.

            indexer = core.linalg.MultiIndexer([2, 2], 'C');

            M = [1, 1; 1, 2; 2, 1; 2, 2];
            L = indexer.multiToLinear(M);

            expected = (1:4)';
            testCase.verifyEqual(L, expected);
        end

        function testMultiToLinearCellInput(testCase)
            % TESTMULTITOLINEARCELLIMPUT Test cell array input for conversion.

            indexer = core.linalg.MultiIndexer([2, 2], 'F');

            M = {1:2, [1; 1; 2; 2]};
            L = indexer.multiToLinear(M);

            expected = [1; 2; 1; 2; 3; 4; 3; 4];
            testCase.verifyEqual(L, expected);
        end

        function testMultiToLinearPeriodicBoundary(testCase)
            % TESTMULTITOLINEARPERIODICBOUNDARY Test periodic boundary handling.

            indexer = core.linalg.MultiIndexer([2, 2], 'F');

            M = [0, 1; 3, 1; 1, 0; 1, 3];
            L = indexer.multiToLinear(M, 0);

            expected = [2; 1; 3; 1];
            testCase.verifyEqual(L, expected);
        end

        function testMultiToLinearStrictBoundary(testCase)
            % TESTMULTITOLINEARSTRICTBOUNDARY Test strict boundary handling.

            indexer = core.linalg.MultiIndexer([2, 2], 'F');

            M = [1, 1; 3, 1; 1, 3; 0, 1];
            L = indexer.multiToLinear(M, 1);

            expected = [1; 0; 0; 0];
            testCase.verifyEqual(L, expected);
        end

        function testLinearToMulti2DF(testCase)
            % TESTLINEARTOMULTI2DF Test 2D linear to multi-index conversion.

            indexer = core.linalg.MultiIndexer([2, 2], 'F');

            L = [1; 2; 3; 4];
            M = indexer.linearToMulti(L);

            expected = [; ...
                1, 1; ...
                2, 1; ...
                1, 2; ...
                2, 2; ...
                ];

            testCase.verifyEqual(M, expected);
        end

        function testLinearToMulti2DC(testCase)
            % TESTLINEARTOMULTI2DC Test 2D linear to multi-index conversion C-style.

            indexer = core.linalg.MultiIndexer([2, 2], 'C');

            L = [1; 2; 3; 4];
            M = indexer.linearToMulti(L);

            expected = [; ...
                1, 1; ...
                1, 2; ...
                2, 1; ...
                2, 2; ...
                ];

            testCase.verifyEqual(M, expected);
        end

        function testInvalidInput(testCase)
            % TESTINVALIDINPUT Test error handling for invalid inputs.

            indexer = core.linalg.MultiIndexer();
            testCase.verifyError(@() indexer.generate(), ...
                'core:linalg:MultiIndexer:MissingShape');

            testCase.verifyError(@() core.linalg.MultiIndexer([2, -1]), ...
                'core:linalg:MultiIndexer:InvalidShape');
            testCase.verifyError(@() core.linalg.MultiIndexer([2, 0]), ...
                'core:linalg:MultiIndexer:InvalidShape');

            indexer = core.linalg.MultiIndexer([2, 2]);

            testCase.verifyError(@() indexer.linearToMulti(5), ...
                'core:linalg:MultiIndexer:IndexOutOfBounds');
            testCase.verifyError(@() indexer.linearToMulti(0), ...
                'core:linalg:MultiIndexer:IndexOutOfBounds');

            testCase.verifyError(@() indexer.multiToLinear({1, 2, 3}), ...
                'core:linalg:MultiIndexer:DimensionMismatch');
        end

        function testRoundTrip(testCase, dimension, style)
            % TESTROUNDTRIP Test round-trip conversion consistency.

            if dimension < 4
                shape = 2:(dimension + 1);

                indexer = core.linalg.MultiIndexer(shape, style);

                M = indexer.generate();
                L = indexer.multiToLinear(M);
                M2 = indexer.linearToMulti(L);

                testCase.verifyEqual(M, M2);

                L2 = indexer.multiToLinear(M2);
                testCase.verifyEqual(L, L2);
            end
        end

        function testLargeShape(testCase)
            % TESTLARGESHAPE Test indexer with larger tensor shapes.

            indexer = core.linalg.MultiIndexer([10, 5], 'F');

            M = indexer.generate();

            testCase.verifyEqual(size(M), [prod([10, 5]), 2]);
            testCase.verifyTrue(all(M(:, 1) >= 1 & M(:, 1) <= 10));
            testCase.verifyTrue(all(M(:, 2) >= 1 & M(:, 2) <= 5));

            L = indexer.multiToLinear(M);
            testCase.verifyEqual(length(L), prod([10, 5]));
            testCase.verifyTrue(all(L >= 1 & L <= prod([10, 5])));

            M2 = indexer.linearToMulti(L);
            testCase.verifyEqual(M, M2);
        end

        function testNonSquareTensor(testCase)
            % TESTNONSQUARETENSOR Test non-square tensor indexing.

            indexer = core.linalg.MultiIndexer([2, 3, 4], 'F');

            M = indexer.generate();

            testCase.verifyEqual(size(M), [prod([2, 3, 4]), 3]);

            testCase.verifyTrue(all(M(:, 1) >= 1 & M(:, 1) <= 2));
            testCase.verifyTrue(all(M(:, 2) >= 1 & M(:, 2) <= 3));
            testCase.verifyTrue(all(M(:, 3) >= 1 & M(:, 3) <= 4));

            L = indexer.multiToLinear(M);
            testCase.verifyEqual(L, (1:prod([2, 3, 4]))');

            M2 = indexer.linearToMulti(L);
            testCase.verifyEqual(M, M2);
        end
    end
end