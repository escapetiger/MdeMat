classdef TestNonuniformGrid < matlab.unittest.TestCase
    
    properties (TestParameter)
        gridType = struct(...
            'linear1D',     {{linspace(0, 1, 6)}}, ...
            'nonlinear1D',  {{[0, 0.1, 0.3, 0.6, 1.0]}}, ...
            'linear2D',     {{linspace(0, 1, 5), linspace(0, 1, 4)}}, ...
            'nonlinear2D',  {{[0, 0.2, 0.5, 0.9, 1], [0, 0.3, 0.7, 1]}}, ...
            'linear3D',     {{linspace(0, 1, 4), linspace(0, 1, 3), linspace(0, 1, 5)}}, ...
            'nonlinear3D',  {{[0, 0.25, 0.75, 1], [0, 0.4, 1], [0, 0.2, 0.6, 0.9, 1]}})
    end
    
    methods (Test)
        function testConstructor(testCase, gridType)
            grid = approx.mesh.NonuniformGrid(gridType);
            
            nDims = length(gridType);
            resolution = zeros(1, nDims);
            for i = 1:nDims
                resolution(i) = length(gridType{i}) - 1;
            end
            
            testCase.verifyEqual(grid.nDims, nDims);
            testCase.verifyEqual(grid.resolution, resolution);
            testCase.verifyEqual(length(grid.centroids), nDims);
            testCase.verifyEqual(length(grid.spacings), nDims);
            testCase.verifyEqual(length(grid.elements), nDims);
            testCase.verifyEqual(length(grid.boundary), 2*nDims);
            
            for i = 1:nDims
                vertices = gridType{i};
                expectedCentroids = (vertices(1:end-1) + vertices(2:end)) / 2;
                testCase.verifyEqual(grid.centroids{i}, expectedCentroids, 'AbsTol', 1e-10);
            end
            
            for i = 1:nDims
                vertices = gridType{i};
                expectedSpacings = diff(vertices);
                testCase.verifyEqual(grid.spacings{i}, expectedSpacings, 'AbsTol', 1e-10);
            end
        end
        
        function testMeasure(testCase, gridType)
            grid = approx.mesh.NonuniformGrid(gridType);
            
            minSpacings = zeros(1, length(gridType));
            for i = 1:length(gridType)
                minSpacings(i) = min(diff(gridType{i}));
            end
            expectedMeasure = min(minSpacings);
            
            testCase.verifyEqual(grid.measure, expectedMeasure, 'AbsTol', 1e-10);
        end
        
        function testNTotalElements(testCase, gridType)
            grid = approx.mesh.NonuniformGrid(gridType);
            
            resolution = zeros(1, length(gridType));
            for i = 1:length(gridType)
                resolution(i) = length(gridType{i}) - 1;
            end
            expectedTotal = prod(resolution);
            
            testCase.verifyEqual(grid.nTotalElements, expectedTotal);
        end
        
        function testAllElementIndices(testCase, gridType)
            grid = approx.mesh.NonuniformGrid(gridType);
            
            indices = grid.allElementMultiIndices;
            
            resolution = zeros(1, length(gridType));
            for i = 1:length(gridType)
                resolution(i) = length(gridType{i}) - 1;
            end
            expectedCount = prod(resolution);
            
            testCase.verifyEqual(size(indices, 1), expectedCount);
            testCase.verifyEqual(size(indices, 2), grid.nDims);
                
            for d = 1:grid.nDims
                testCase.verifyTrue(all(indices(:,d) >= 1 & indices(:,d) <= grid.resolution(d)));
            end
            
            uniqueIndices = unique(indices, 'rows');
            testCase.verifyEqual(size(uniqueIndices, 1), size(indices, 1));
        end
        
        function testBoundaryElementIndices(testCase, gridType)
            grid = approx.mesh.NonuniformGrid(gridType);
            
            boundaryIndices = grid.boundaryElementMultiIndices;
            for iBoundary = 1:2*grid.nDims
                indices = boundaryIndices{iBoundary};
                
                dim = ceil(iBoundary/2);
                isLower = mod(iBoundary, 2) == 1;
                
                expectedCount = prod(grid.resolution) / grid.resolution(dim);
                
                testCase.verifyEqual(size(indices, 1), expectedCount);
                
                if isLower
                    testCase.verifyTrue(all(indices(:, dim) == 1));
                else
                    testCase.verifyTrue(all(indices(:, dim) == grid.resolution(dim)));
                end
            end
        end
        
        function testGraph(testCase, gridType)
            grid = approx.mesh.NonuniformGrid(gridType);
            
            graph = grid.graph;
            
            testCase.verifyEqual(graph.nDims, grid.nDims);
                
            nVertices = graph.nVertices;
            for d = 1:grid.nDims
                testCase.verifyEqual(nVertices(d), grid.resolution(d) + 1);
            end
            
            nEdges = graph.nEdges;
            for d = 1:grid.nDims
                testCase.verifyEqual(nEdges(d), grid.resolution(d));
            end
        end
        
        function testMagnitudes(testCase, gridType)
            grid = approx.mesh.NonuniformGrid(gridType);
            
            magnitudes = grid.magnitudes;
            
            testCase.verifyTrue(isnumeric(magnitudes));
            
            expectedCount = prod(grid.resolution);
            testCase.verifyEqual(length(magnitudes), expectedCount);
            
            testCase.verifyTrue(all(magnitudes > 0));
            
            if all(cellfun(@(x) length(unique(diff(x))) == 1, gridType))
                expectedMagnitude = prod(cellfun(@(x) diff(x(1:2)), gridType));
                testCase.verifyEqual(magnitudes, expectedMagnitude * ones(size(magnitudes)), ...
                    'AbsTol', 1e-10);
            end
        end
        
        function testCollocate(testCase, gridType)
            grid = approx.mesh.NonuniformGrid(gridType);
            
            X = cell(1, grid.nDims);
            for d = 1:grid.nDims
                X{d} = 0;
            end
            
            Y = grid.collocate(X);
            
            testCase.verifyTrue(iscell(Y));
            testCase.verifyEqual(length(Y), grid.nDims);
            
            if grid.nDims == 1
                indices = [1];
            elseif grid.nDims == 2
                indices = [1, 1];
            else
                indices = [1, 1, 1];
            end
            
            Y2 = grid.collocate(X, indices);
            
            testCase.verifyTrue(ismatrix(Y2));
            testCase.verifyEqual(size(Y2, 2), 1);
            testCase.verifyEqual(size(Y2, 1), grid.nDims);
            
            for d = 1:grid.nDims
                expectedCoord = grid.centroids{d}(indices(d));
                testCase.verifyEqual(Y2(d), expectedCoord, 'AbsTol', 1e-10);
            end
        end
        
        function testNonuniformCollocate(testCase)
            gridCoords = {[0, 0.25, 0.75, 1.0]};
            grid = approx.mesh.NonuniformGrid(gridCoords);
            
            refNodes = {[-0.5, 0, 0.5]};
            
            indices = [2];
            
            Y = grid.collocate(refNodes, indices);
            
            h = grid.spacings{1}(indices);
            c = grid.centroids{1}(indices);
            expectedCoords = c + h * refNodes{1};
            
            testCase.verifyEqual(Y(1,:), expectedCoords, 'AbsTol', 1e-10);
        end
    end
end