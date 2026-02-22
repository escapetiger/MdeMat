classdef TestLxMultiIndexer < matlab.unittest.TestCase
    % TESTLXMULTIINDEXER Unit tests for LxMultiIndexer class.
    %
    %   TestLxMultiIndexer provides test coverage for LxMultiIndexer
    %   functionality including index generation, caching, and conversion
    %   between multi-indices and linear indices.
    %
    % Examples:
    %   % Run all tests
    %   results = runtests(core.linalg.tests.TestLxMultiIndexer);
    %
    % See Also:
    %   core.linalg.LxMultiIndexer
    
    properties (TestParameter)
        Style = {'F', 'C'}   % Storage ordering styles
        Dimension = {1, 2, 3} % Test dimensions
        Upper = {1, 2, 3}    % Test upper bounds
    end
    
    methods (Test)
        function test1DIndices(testCase, Style)
            % TEST1DINDICES Test 1D index generation.

            indexer = core.linalg.LxMultiIndexer(style=Style);
            
            indexer.setCache(1, 3);
            M = indexer.Cache;
            
            expected = (1:3)';
            testCase.verifyEqual(M, expected);
        end
        
        function test2DIndicesF(testCase)
            % TEST2DINDICESF Test 2D indices with F-style ordering.

            indexer = core.linalg.LxMultiIndexer(style='F');
            
            indexer.setCache(2, 2);
            M = indexer.Cache;
            
            expected = [
                1, 1;
                2, 1;
                1, 2;
                2, 2
            ];
            
            M = sortrows(M);
            expected = sortrows(expected);
            testCase.verifyEqual(M, expected);
        end
        
        function test2DIndicesC(testCase)
            % TEST2DINDICESC Test 2D indices with C-style ordering.

            indexer = core.linalg.LxMultiIndexer(style='C');
            
            indexer.setCache(2, 2);
            M = indexer.Cache;
            
            expected = [
                1, 1;
                1, 2;
                2, 1;
                2, 2
            ];
            
            M = sortrows(M);
            expected = sortrows(expected);
            testCase.verifyEqual(M, expected);
        end
        
        function test3DIndices(testCase)
            % TEST3DINDICES Test 3D index generation.

            indexer = core.linalg.LxMultiIndexer(style='F');
            
            indexer.setCache(3, 2);
            M = indexer.Cache;
            
            testCase.verifyEqual(size(M), [8, 3]);
            testCase.verifyTrue(all(M(:) >= 1 & M(:) <= 2));
            
            expectedIndices = [1, 2];
            [X, Y, Z] = ndgrid(expectedIndices, expectedIndices, expectedIndices);
            expected = [X(:), Y(:), Z(:)];
            
            M = sortrows(M);
            expected = sortrows(expected);
            testCase.verifyEqual(M, expected);
        end
        
        function testInvalidInput(testCase)
            % TESTINVALIDINPUT Test error handling for invalid inputs.

            indexer = core.linalg.LxMultiIndexer();
            
            testCase.verifyError(@() indexer.setCache(2, -1), ...
                'MATLAB:validators:mustBePositive');
            testCase.verifyError(@() indexer.setCache(2, 0), ...
                'MATLAB:validators:mustBePositive');
        end
        
        function testMemoryAdaptation(testCase)
            % TESTMEMORYADAPTATION Test memory scaling with different sizes.

            indexer = core.linalg.LxMultiIndexer(style='F');
            
            indexer.setCache(3, 2);
            M_small = indexer.Cache;
            testCase.verifyEqual(size(M_small), [2^3, 3]);
            
            indexer.setCache(3, 3);
            M_medium = indexer.Cache;
            testCase.verifyEqual(size(M_medium), [3^3, 3]);
            testCase.verifyTrue(all(M_medium(:) >= 1 & M_medium(:) <= 3));
        end
        
        function testAllSizeCombinations(testCase, Dimension, Upper, Style)
            % TESTALLSIZECOMBINATIONS Test various dimension/upper combinations.

            if Dimension < 4 && Upper < 4
                indexer = core.linalg.LxMultiIndexer(style=Style);
                indexer.setCache(Dimension, Upper);
                M = indexer.Cache;
                
                testCase.verifyEqual(size(M, 2), Dimension);
                testCase.verifyEqual(size(M, 1), Upper^Dimension);
                testCase.verifyTrue(all(M(:) >= 1 & M(:) <= Upper));
            end
        end
        
        function testResultConsistency(testCase)
            % TESTRESULTCONSISTENCY Test consistent results across calls.

            indexer = core.linalg.LxMultiIndexer(style='F');
            
            indexer.setCache(2, 3);
            M1 = indexer.Cache;
            
            indexer.setCache(2, 3);
            M2 = indexer.Cache;
            
            testCase.verifyEqual(M1, M2);
        end
        
        function testStorageOrderDifference(testCase)
            % TESTSTORAGEORDERDIFFERENCE Test F vs C style ordering.

            indexerF = core.linalg.LxMultiIndexer(style='F');
            indexerC = core.linalg.LxMultiIndexer(style='C');
            
            indexerF.setCache(2, 3);
            indexerC.setCache(2, 3);
            
            MF = indexerF.Cache;
            MC = indexerC.Cache;
            
            testCase.verifyEqual(size(MF), size(MC));
            
            MF_sorted = sortrows(MF);
            MC_sorted = sortrows(MC);
            
            testCase.verifyEqual(MF_sorted, MC_sorted);
        end
        
        function testLowerUpperBounds(testCase)
            % TESTLOWERUPPERBOUNDS Test with both lower and upper bounds.

            indexer = core.linalg.LxMultiIndexer(style='F');
            
            indexer.setCache(2, 3, minIdx=2);
            M = indexer.Cache;
            
            expected = [
                2, 2;
                3, 2;
                2, 3;
                3, 3
            ];
            
            M = sortrows(M);
            expected = sortrows(expected);
            
            testCase.verifyEqual(M, expected);
        end
        
        function testVectorLowerUpper(testCase)
            % TESTVECTORLOWERUPPER Test with vector lower and upper bounds.

            indexer = core.linalg.LxMultiIndexer(style='F');
            
            lower = [1, 2];
            upper = [2, 3];
            indexer.setCache(2, upper, minIdx = lower);
            M = indexer.Cache;
            
            expected = [
                1, 2;
                2, 2;
                1, 3;
                2, 3
            ];
            
            M = sortrows(M);
            expected = sortrows(expected);
            
            testCase.verifyEqual(M, expected);
        end
        
        function testCacheAccess(testCase)
            % TESTCACHEACCESS Test cache property access.

            indexer = core.linalg.LxMultiIndexer(style='F');
            
            indexer.setCache(2, 3);
            
            testCase.verifyEqual(indexer.NDims, 2);
            testCase.verifyEqual(indexer.NIndices, 9);
            testCase.verifyEqual(size(indexer.Cache), [9, 2]);
        end
        
        function testRoundTripConversion(testCase)
            % TESTROUNDTRIPCONVERSION Test round-trip index conversion.

            indexer = core.linalg.LxMultiIndexer(style='F');
            
            indexer.setCache(2, 3);
            
            original = indexer.Cache;
            linear = indexer.multiToLinear(original);
            reconstructed = indexer.linearToMulti(linear);
            
            testCase.verifyEqual(reconstructed, original);
        end
    end
end