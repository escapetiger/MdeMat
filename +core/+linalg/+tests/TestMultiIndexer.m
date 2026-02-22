classdef TestMultiIndexer < matlab.unittest.TestCase
    % TESTMULTIINDEXER Unit tests for MultiIndexer class.
    %
    %   TestMultiIndexer provides test coverage for MultiIndexer
    %   functionality including index generation and conversion between
    %   multi-indices and linear indices with various padding modes.
    %
    % Examples:
    %   % Run all tests
    %   results = runtests(core.linalg.tests.TestMultiIndexer);
    %
    % See Also:
    %   core.linalg.MultiIndexer

    properties (TestParameter)
        Style = {'F', 'C'}   % Storage ordering styles
        Dimension = {1, 2, 3} % Test dimensions
        Pad = {'empty', 'wrap', 'edge', 'reflect', 'symmetric'} % Padding modes
    end

    methods (Test)
        function test1DGenerate(testCase, Style)
            % TEST1DGENERATE Test 1D index generation.

            indexer = core.linalg.MultiIndexer(shape=3, style=Style);

            M = indexer.generate();

            expected = (1:3)';
            testCase.verifyEqual(M, expected);
        end

        function test2DGenerateF(testCase)
            % TEST2DGENERATEF Test 2D index generation with F-style.

            indexer = core.linalg.MultiIndexer(shape=[2, 3], style='F');

            M = indexer.generate();

            expected = [1, 1; ...
                        2, 1; ...
                        1, 2; ...
                        2, 2; ...
                        1, 3; ...
                        2, 3];

            testCase.verifyEqual(M, expected);
        end

        function test2DGenerateC(testCase)
            % TEST2DGENERATEC Test 2D index generation with C-style.

            indexer = core.linalg.MultiIndexer(shape=[2, 3], style='C');

            M = indexer.generate();

            expected = [1, 1; ...
                        1, 2; ...
                        1, 3; ...
                        2, 1; ...
                        2, 2; ...
                        2, 3];

            testCase.verifyEqual(M, expected);
        end

        function testMultiToLinear2DF(testCase)
            % TESTMULTITOLINEAR2DF Test 2D multi-index to linear conversion.

            indexer = core.linalg.MultiIndexer(shape=[2, 2], style='F');

            M = [1, 1; 2, 1; 1, 2; 2, 2];
            L = indexer.multiToLinear(M);

            expected = (1:4)';
            testCase.verifyEqual(L, expected);
        end

        function testMultiToLinear2DC(testCase)
            % TESTMULTITOLINEAR2DC Test 2D multi-index to linear conversion C-style.

            indexer = core.linalg.MultiIndexer(shape=[2, 2], style='C');

            M = [1, 1; 1, 2; 2, 1; 2, 2];
            L = indexer.multiToLinear(M);

            expected = (1:4)';
            testCase.verifyEqual(L, expected);
        end

        function testMultiToLinearCellInput(testCase)
            % TESTMULTITOLINEARCELLIMPUT Test cell array input for conversion.

            indexer = core.linalg.MultiIndexer(shape=[2, 2], style='F');

            F = {1:2, [1; 1; 2; 2]};
            M = indexer.factorToMulti(F);
            L = indexer.multiToLinear(M);

            expected = [1; 2; 1; 2; 3; 4; 3; 4];
            testCase.verifyEqual(L, expected);
        end

        function testMultiToLinearPaddingEmpty(testCase)
            % TESTMULTITOLINEARPPADDINGEMPTY Test empty padding mode.

            indexer = core.linalg.MultiIndexer(shape=4, style='F');

            % Test with out-of-bounds indices
            M = [-1; 0; 1; 2; 3; 4; 5; 6];
            L = indexer.multiToLinear(M, Pad='empty');

            expected = [0; 0; 1; 2; 3; 4; 0; 0];
            testCase.verifyEqual(L, expected);
        end

        function testMultiToLinearPaddingWrap(testCase)
            % TESTMULTITOLINEARPADDINGWRAP Test wrap padding mode.

            indexer = core.linalg.MultiIndexer(shape=4, style='F');

            % Test with out-of-bounds indices
            M = [-1; 0; 1; 2; 3; 4; 5; 6];
            L = indexer.multiToLinear(M, Pad='wrap');

            expected = [3; 4; 1; 2; 3; 4; 1; 2];
            testCase.verifyEqual(L, expected);
        end

        function testMultiToLinearPaddingEdge(testCase)
            % TESTMULTITOLINEARPPADDINGEDGE Test edge padding mode.

            indexer = core.linalg.MultiIndexer(shape=4, style='F');

            % Test with out-of-bounds indices
            M = [-1; 0; 1; 2; 3; 4; 5; 6];
            L = indexer.multiToLinear(M, Pad='edge');

            expected = [1; 1; 1; 2; 3; 4; 4; 4];
            testCase.verifyEqual(L, expected);
        end

        function testMultiToLinearPaddingReflect(testCase)
            % TESTMULTITOLINEARPPADDINGREFLECT Test reflect padding mode.

            indexer = core.linalg.MultiIndexer(shape=4, style='F');

            % Test with out-of-bounds indices
            M = [-1; 0; 1; 2; 3; 4; 5; 6];
            L = indexer.multiToLinear(M, Pad='reflect');

            expected = [3; 2; 1; 2; 3; 4; 3; 2];
            testCase.verifyEqual(L, expected);
        end

        function testMultiToLinearPaddingSymmetric(testCase)
            % TESTMULTITOLINEARPPADDDINGSYMMETRIC Test symmetric padding mode.

            indexer = core.linalg.MultiIndexer(shape=4, style='F');

            % Test with out-of-bounds indices
            M = [-1; 0; 1; 2; 3; 4; 5; 6];
            L = indexer.multiToLinear(M, Pad='symmetric');

            expected = [2; 1; 1; 2; 3; 4; 4; 3];
            testCase.verifyEqual(L, expected);
        end

        function testMultiToLinear2DPads(testCase, Pad)
            % TESTMULTITOLINEAR2DPADDINGMODES Test 2D padding modes.

            indexer = core.linalg.MultiIndexer(shape=[3, 3], style='F');

            % Test with multi-dimensional out-of-bounds indices
            M = [0, 0; 0, 1; 1, 0; 4, 4; 4, 1; 1, 4];
            L = indexer.multiToLinear(M, pad=Pad);

            % Verify no NaN or Inf values
            testCase.verifyTrue(all(isfinite(L)));
            
            % Verify output is numeric and has correct size
            testCase.verifyTrue(isnumeric(L));
            testCase.verifyEqual(size(L), [size(M, 1), 1]);

            % For empty padding, out-of-bounds should be 0
            if strcmp(Pad, 'empty')
                testCase.verifyTrue(all(L == 0 | (L >= 1 & L <= 9)));
            else
                % For other padding modes, all values should be valid indices
                testCase.verifyTrue(all(L >= 1 & L <= 9));
            end
        end

        function testMultiToLinearDefaultPadding(testCase)
            % TESTMULTITOLINEARDEFAULTPADDING Test default padding behavior.

            indexer = core.linalg.MultiIndexer(shape=[2, 2], style='F');

            M = [0, 1; 3, 1; 1, 0; 1, 3];
            
            % Test with no padding argument (should default to 'empty')
            L1 = indexer.multiToLinear(M);
            L2 = indexer.multiToLinear(M, Pad='empty');

            testCase.verifyEqual(L1, L2);
            testCase.verifyEqual(L1, [0; 0; 0; 0]);
        end

        function testLinearToMulti2DF(testCase)
            % TESTLINEARTOMULTI2DF Test 2D linear to multi-index conversion.

            indexer = core.linalg.MultiIndexer(shape=[2, 2], style='F');

            L = [1; 2; 3; 4];
            M = indexer.linearToMulti(L);

            expected = [1, 1; ...
                        2, 1; ...
                        1, 2; ...
                        2, 2];

            testCase.verifyEqual(M, expected);
        end

        function testLinearToMulti2DC(testCase)
            % TESTLINEARTOMULTI2DC Test 2D linear to multi-index conversion C-style.

            indexer = core.linalg.MultiIndexer(shape=[2, 2], style='C');

            L = [1; 2; 3; 4];
            M = indexer.linearToMulti(L);

            expected = [1, 1; ...
                        1, 2; ...
                        2, 1; ...
                        2, 2];

            testCase.verifyEqual(M, expected);
        end

        function testInvalidInput(testCase)
            % TESTINVALIDINPUT Test error handling for invalid inputs.

            indexer = core.linalg.MultiIndexer();
            testCase.verifyError(@() indexer.generate(), ...
                'core:linalg:MultiIndexer:MissingShape');

            testCase.verifyError(@() core.linalg.MultiIndexer(shape=[2, -1]), ...
                'MATLAB:validators:mustBePositive');
            testCase.verifyError(@() core.linalg.MultiIndexer(shape=[2, 0]), ...
                'MATLAB:validators:mustBePositive');

            indexer = core.linalg.MultiIndexer(shape=[2, 2]);

            testCase.verifyError(@() indexer.linearToMulti(5), ...
                'core:linalg:MultiIndexer:IndexOutOfBounds');
            testCase.verifyError(@() indexer.linearToMulti(0), ...
                'MATLAB:validators:mustBePositive');
        end

        function testRoundTrip(testCase, Dimension, Style)
            % TESTROUNDTRIP Test round-trip conversion consistency.

            if Dimension < 4
                shape = 2:(Dimension + 1);

                indexer = core.linalg.MultiIndexer(shape=shape, style=Style);

                M = indexer.generate();
                L = indexer.multiToLinear(M);
                M2 = indexer.linearToMulti(L);

                testCase.verifyEqual(M, M2);

                L2 = indexer.multiToLinear(M2);
                testCase.verifyEqual(L, L2);
            end
        end

        function testRoundTripWithPadding(testCase, Pad)
            % TESTROUNDTRIPWITHPADDING Test round-trip with padding modes.

            indexer = core.linalg.MultiIndexer(shape=[3, 3], style='F');

            % Generate valid indices
            M = indexer.generate();
            L = indexer.multiToLinear(M, pad=Pad);
            
            % Only test round-trip for non-zero indices (valid for empty padding)
            if strcmp(Pad, 'empty')
                validMask = L > 0;
                M_valid = M(validMask, :);
                L_valid = L(validMask);
                M2 = indexer.linearToMulti(L_valid);
                testCase.verifyEqual(M_valid, M2);
            else
                M2 = indexer.linearToMulti(L);
                testCase.verifyEqual(M, M2);
            end
        end

        function testLargeShape(testCase)
            % TESTLARGESHAPE Test indexer with larger tensor shapes.

            indexer = core.linalg.MultiIndexer(shape=[10, 5], style='F');

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

            indexer = core.linalg.MultiIndexer(shape=[2, 3, 4], style='F');

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

        function testPaddingEdgeCases(testCase)
            % TESTPADDINGEDGECASES Test edge cases for padding modes.

            indexer = core.linalg.MultiIndexer(shape=[2, 2], style='F');

            % Test extremely out-of-bounds values
            M_extreme = [-100, 1; 100, 1; 1, -100; 1, 100];
            
            % Test each padding mode handles extreme values
            paddingModes = {'empty', 'wrap', 'edge', 'reflect', 'symmetric'};
            for iMode = 1:length(paddingModes)
                mode = paddingModes{iMode};
                L = indexer.multiToLinear(M_extreme, Pad=mode);
                
                % Verify output is finite and has correct size
                testCase.verifyTrue(all(isfinite(L)), ...
                    sprintf('Padding mode %s produced non-finite values', mode));
                testCase.verifyEqual(size(L), [size(M_extreme, 1), 1], ...
                    sprintf('Padding mode %s produced wrong output size', mode));
                
                % For non-empty modes, all values should be valid indices
                if ~strcmp(mode, 'empty')
                    testCase.verifyTrue(all(L >= 1 & L <= 4), ...
                        sprintf('Padding mode %s produced invalid indices', mode));
                end
            end
        end
    end
end