classdef TestUniformGrid < matlab.unittest.TestCase
    
    properties (TestParameter)
        resolution = {[5], [5, 4], [3, 4, 5]}
    end
    
    methods (Test)
        function testConstructor(testCase, resolution)
            nDims = length(resolution);
            n = 5 * ones(1, nDims);
            bbox = repmat([0, 1], 1, nDims);
            
            grid = approx.mesh.UniformGrid(n, bbox);
            
            testCase.verifyEqual(grid.nDims, nDims);
            testCase.verifyEqual(grid.resolution, n);
            testCase.verifyEqual(length(grid.centroids), nDims);
            testCase.verifyEqual(length(grid.spacings), nDims);
            testCase.verifyEqual(length(grid.elements), nDims);
            testCase.verifyEqual(length(grid.boundary), 2*nDims);
        end
        
        function testConstructorValidation(testCase)
            % Test dimension mismatch between n and bbox
            testCase.verifyError(...
                @() approx.mesh.UniformGrid([1, 2], [0, 1]), ...
                'approx:mesh:UniformGrid:DimensionMismatch');
            
            testCase.verifyError(...
                @() approx.mesh.UniformGrid([1, 2], [0, 1, 0]), ...
                'approx:mesh:UniformGrid:DimensionMismatch');
            
            % Test invalid bounds (upper <= lower)
            testCase.verifyError(...
                @() approx.mesh.UniformGrid([5], [1, 0]), ...
                'approx:mesh:UniformGrid:InvalidInput');
            
            testCase.verifyError(...
                @() approx.mesh.UniformGrid([3, 4], [0, 1, 2, 1]), ...
                'approx:mesh:UniformGrid:InvalidInput');
            
            % Test invalid n values (non-positive or non-integer)
            testCase.verifyError(...
                @() approx.mesh.UniformGrid([0], [0, 1]), ...
                'approx:mesh:UniformGrid:InvalidInput');
            
            testCase.verifyError(...
                @() approx.mesh.UniformGrid([2.5], [0, 1]), ...
                'approx:mesh:UniformGrid:InvalidInput');
        end
        
        function testMeasure(testCase, resolution)
            nDims = length(resolution);
            n = 5 * ones(1, nDims);
            bbox = repmat([0, 1], 1, nDims);
            
            grid = approx.mesh.UniformGrid(n, bbox);
            
            expectedMeasure = 0.2;
            
            testCase.verifyEqual(grid.measure, expectedMeasure, 'AbsTol', 1e-10);
        end
        
        function testNTotalElements(testCase, resolution)
            nDims = length(resolution);
            bbox = repmat([0, 1], 1, nDims);
            
            grid = approx.mesh.UniformGrid(resolution, bbox);
            
            expectedTotal = prod(resolution);
            
            testCase.verifyEqual(grid.nTotalElements, expectedTotal);
        end
        
        function testAllElementIndices(testCase, resolution)
            nDims = length(resolution);
            bbox = repmat([0, 1], 1, nDims);
            
            grid = approx.mesh.UniformGrid(resolution, bbox);
            
            indices = grid.allElementMultiIndices;
            
            testCase.verifyEqual(size(indices, 1), prod(resolution));
            testCase.verifyEqual(size(indices, 2), nDims);
            
            for d = 1:nDims
                testCase.verifyTrue(all(indices(:,d) >= 1 & indices(:,d) <= resolution(d)));
            end
            
            uniqueIndices = unique(indices, 'rows');
            testCase.verifyEqual(size(uniqueIndices, 1), size(indices, 1));
        end
        
        function testBoundaryElementIndices(testCase, resolution)
            nDims = length(resolution);
            bbox = repmat([0, 1], 1, nDims);
            
            grid = approx.mesh.UniformGrid(resolution, bbox);
            
            boundaryIndices = grid.boundaryElementMultiIndices;
            for iBoundary = 1:2*nDims
                indices = boundaryIndices{iBoundary};
                
                dim = ceil(iBoundary/2);
                isLower = mod(iBoundary, 2) == 1;
                
                expectedCount = prod(resolution) / resolution(dim);
                
                testCase.verifyEqual(size(indices, 1), expectedCount);
                
                if isLower
                    testCase.verifyTrue(all(indices(:, dim) == 1));
                else
                    testCase.verifyTrue(all(indices(:, dim) == resolution(dim)));
                end
            end
        end
        
        function testGraph(testCase, resolution)
            nDims = length(resolution);
            n = 4 * ones(1, nDims);
            bbox = repmat([0, 1], 1, nDims);
            
            grid = approx.mesh.UniformGrid(n, bbox);
            graph = grid.graph;
            
            testCase.verifyEqual(graph.nDims, nDims);
            
            nVertices = graph.nVertices;
            for d = 1:nDims
                testCase.verifyEqual(nVertices(d), n(d) + 1);
            end
            
            nEdges = graph.nEdges;
            for d = 1:nDims
                testCase.verifyEqual(nEdges(d), n(d));
            end
        end
        
        function testMagnitudes(testCase, resolution)
            nDims = length(resolution);
            n = 5 * ones(1, nDims);
            bbox = repmat([0, 1], 1, nDims);
            
            grid = approx.mesh.UniformGrid(n, bbox);
            
            expectedMagnitude = (0.2)^nDims;
            
            testCase.verifyEqual(grid.magnitudes, expectedMagnitude, 'AbsTol', 1e-10);
        end
        
        function testCollocate(testCase, resolution)
            nDims = length(resolution);
            n = 3 * ones(1, nDims);
            bbox = repmat([0, 1], 1, nDims);
            
            grid = approx.mesh.UniformGrid(n, bbox);
            
            X = cell(1, nDims);
            for d = 1:nDims
                X{d} = 0;
            end
            
            Y = grid.collocate(X);
            
            testCase.verifyTrue(iscell(Y));
            testCase.verifyEqual(length(Y), nDims);
            
            if nDims == 1
                indices = [2];
            elseif nDims == 2
                indices = [2, 2];
            else
                indices = [2, 2, 2];
            end
            
            Y2 = grid.collocate(X, indices);
            
            testCase.verifyTrue(ismatrix(Y2));
            testCase.verifyEqual(size(Y2, 2), 1);
            testCase.verifyEqual(size(Y2, 1), nDims);
            
            for d = 1:nDims
                expectedCoord = grid.centroids{d}(indices(d));
                testCase.verifyEqual(Y2(d), expectedCoord, 'AbsTol', 1e-10);
            end
        end
        
        function testCollocateEdgeCases(testCase)
            n = [2, 3];
            bbox = [0, 2, 0, 3];
            
            grid = approx.mesh.UniformGrid(n, bbox);
            
            X = {[0.5, -0.5], [0, 0.5]};
            I = [1, 1; 2, 3];
            Y = grid.collocate(X, I);
            
            testCase.verifyEqual(size(Y), [2, 4]);
            
            expectedY = zeros(2, 4);
            for i = 1:2
                expectedY(1, 2*i-1:2*i) = grid.centroids{1}(I(i, 1)) + grid.spacings{1} * X{1};
                expectedY(2, 2*i-1:2*i) = grid.centroids{2}(I(i, 2)) + grid.spacings{2} * X{2};
            end
            
            testCase.verifyEqual(Y, expectedY, 'AbsTol', 1e-10);
        end
        
        function testDifferentBoundingBoxes(testCase)
            % Test with non-unit bounding boxes
            n = [4, 3];
            bbox = [-1, 2, 0.5, 3.5];
            
            grid = approx.mesh.UniformGrid(n, bbox);
            
            testCase.verifyEqual(grid.nDims, 2);
            testCase.verifyEqual(grid.resolution, n);
            
            % Check spacing calculations
            expectedSpacing1 = (2 - (-1)) / 4;  % 0.75
            expectedSpacing2 = (3.5 - 0.5) / 3;  % 1.0
            
            testCase.verifyEqual(grid.spacings{1}, expectedSpacing1, 'AbsTol', 1e-10);
            testCase.verifyEqual(grid.spacings{2}, expectedSpacing2, 'AbsTol', 1e-10);
            
            % Check centroid calculations
            expectedCentroids1 = -1 + expectedSpacing1/2 + (0:3) * expectedSpacing1;
            expectedCentroids2 = 0.5 + expectedSpacing2/2 + (0:2) * expectedSpacing2;
            
            testCase.verifyEqual(grid.centroids{1}, expectedCentroids1, 'AbsTol', 1e-10);
            testCase.verifyEqual(grid.centroids{2}, expectedCentroids2, 'AbsTol', 1e-10);
        end
        
        function testConstructorWithDifferentDimensions(testCase)
            % Test 1D case
            grid1D = approx.mesh.UniformGrid([10], [-2, 3]);
            testCase.verifyEqual(grid1D.nDims, 1);
            testCase.verifyEqual(grid1D.resolution, [10]);
            testCase.verifyEqual(grid1D.spacings{1}, 0.5, 'AbsTol', 1e-10);
            
            % Test 3D case
            grid3D = approx.mesh.UniformGrid([5, 4, 6], [0, 1, -1, 2, 0.5, 1.5]);
            testCase.verifyEqual(grid3D.nDims, 3);
            testCase.verifyEqual(grid3D.resolution, [5, 4, 6]);
            testCase.verifyEqual(grid3D.spacings{1}, 0.2, 'AbsTol', 1e-10);
            testCase.verifyEqual(grid3D.spacings{2}, 0.75, 'AbsTol', 1e-10);
            testCase.verifyEqual(grid3D.spacings{3}, 1/6, 'AbsTol', 1e-10);
        end
    end
end