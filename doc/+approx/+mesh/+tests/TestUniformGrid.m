classdef TestUniformGrid < matlab.unittest.TestCase

    properties (TestParameter)
        resolution = {[5], [5, 4], [3, 4, 5]}
    end

    methods (Test)
        function testConstructor(testCase, resolution)
            nDims = length(resolution);
            bbox = repmat([0, 1], 1, nDims);

            grid = approx.mesh.UniformGrid(bbox, resolution);

            testCase.verifyEqual(grid.nDims, nDims);
            testCase.verifyEqual(grid.resolution, resolution);
            testCase.verifyEqual(length(grid.centroids), nDims);
            testCase.verifyEqual(length(grid.spacings), nDims);
            testCase.verifyEqual(length(grid.elements), nDims);
            testCase.verifyEqual(length(grid.boundary), 2*nDims);
            
            %< Verify bbox is stored correctly
            testCase.verifyEqual(grid.bbox, bbox);
            
            %< Verify spacing calculations
            for d = 1:nDims
                expectedSpacing = 1 / resolution(d);
                testCase.verifyEqual(grid.spacings{d}, expectedSpacing, 'AbsTol', 1e-10);
            end
            
            %< Verify centroid calculations
            for d = 1:nDims
                h = 1 / resolution(d);
                expectedCentroids = h/2 + (0:resolution(d)-1) * h;
                testCase.verifyEqual(grid.centroids{d}, expectedCentroids, 'AbsTol', 1e-10);
            end
        end

        function testConstructorValidation(testCase)
            %< Test invalid bounds (upper <= lower)
            testCase.verifyError( ...
                @() approx.mesh.UniformGrid([1, 0], [5]), ...
                'approx:mesh:UniformGrid:InvalidInput');

            testCase.verifyError( ...
                @() approx.mesh.UniformGrid([0, 1, 2, 1], [3, 4]), ...
                'approx:mesh:UniformGrid:InvalidInput');
        end

        function testMeasure(testCase, resolution)
            nDims = length(resolution);
            bbox = repmat([0, 1], 1, nDims);

            grid = approx.mesh.UniformGrid(bbox, resolution);

            expectedMeasure = 1 / max(resolution);

            testCase.verifyEqual(grid.computeMeasure(), expectedMeasure, 'AbsTol', 1e-10);
        end

        function testNElements(testCase, resolution)
            nDims = length(resolution);
            bbox = repmat([0, 1], 1, nDims);

            grid = approx.mesh.UniformGrid(bbox, resolution);

            expectedTotal = prod(resolution);

            testCase.verifyEqual(grid.nElements, expectedTotal);
        end

        function testLinearToMultiIndices(testCase, resolution)
            nDims = length(resolution);
            bbox = repmat([0, 1], 1, nDims);

            grid = approx.mesh.UniformGrid(bbox, resolution);

            %< Test conversion for a few linear indices
            nElements = grid.nElements;
            testIndices = [1, min(5, nElements), nElements];

            for linearIdx = testIndices
                multiIdx = grid.linearToMulti(linearIdx);
                backToLinear = grid.multiToLinear(multiIdx);
                
                testCase.verifyEqual(backToLinear, linearIdx);
                testCase.verifyEqual(size(multiIdx, 1), 1);
                testCase.verifyEqual(size(multiIdx, 2), nDims);

                %< Verify indices are within bounds
                for d = 1:nDims
                    testCase.verifyTrue(multiIdx(d) >= 1 && multiIdx(d) <= resolution(d));
                end
            end
        end

        function testFindInteriorElements(testCase, resolution)
            nDims = length(resolution);
            bbox = repmat([0, 1], 1, nDims);

            grid = approx.mesh.UniformGrid(bbox, resolution);

            %< Skip test for grids with resolution <= 2 in any dimension
            if any(resolution <= 2)
                return;
            end

            interiorElements = grid.findInteriorElements();
            interiorMultiIndices = grid.linearToMulti(interiorElements);

            %< Verify interior elements don't touch boundaries
            for d = 1:nDims
                testCase.verifyTrue(all(interiorMultiIndices(:, d) > 1));
                testCase.verifyTrue(all(interiorMultiIndices(:, d) < resolution(d)));
            end

            %< Verify expected count
            expectedCount = prod(max(0, resolution - 2));
            testCase.verifyEqual(length(interiorElements), expectedCount);
        end

        function testFindBoundaryElements(testCase, resolution)
            nDims = length(resolution);
            bbox = repmat([0, 1], 1, nDims);

            grid = approx.mesh.UniformGrid(bbox, resolution);

            for iBoundary = 1:2*nDims
                boundaryElements = grid.findBoundaryElements(iBoundary);
                boundaryMultiIndices = grid.linearToMulti(boundaryElements);

                dim = ceil(iBoundary/2);
                isLower = mod(iBoundary, 2) == 1;

                expectedCount = prod(resolution) / resolution(dim);
                testCase.verifyEqual(length(boundaryElements), expectedCount);

                if isLower
                    testCase.verifyTrue(all(boundaryMultiIndices(:, dim) == 1));
                else
                    testCase.verifyTrue(all(boundaryMultiIndices(:, dim) == resolution(dim)));
                end
            end
        end

        function testFindNeighborElements(testCase, resolution)
            nDims = length(resolution);
            bbox = repmat([0, 1], 1, nDims);

            grid = approx.mesh.UniformGrid(bbox, resolution);

            %< Test neighbor finding for element in the middle (if exists)
            if all(resolution >= 3)
                centerMultiIdx = ceil(resolution / 2);
                centerLinearIdx = grid.multiToLinear(centerMultiIdx);

                for iBoundary = 1:2*nDims
                    %< Use periodic boundary condition (bc = 0) for interior elements
                    neighborLinearIdx = grid.findNeighborElements(iBoundary, centerLinearIdx, 'periodic');
                    
                    if ~isempty(neighborLinearIdx) && neighborLinearIdx > 0
                        neighborMultiIdx = grid.linearToMulti(neighborLinearIdx);

                        dim = ceil(iBoundary/2);
                        expectedDiff = -(mod(iBoundary, 2) * 2 - 1); %< -1 for odd, +1 for even

                        expectedNeighborMultiIdx = centerMultiIdx;
                        expectedNeighborMultiIdx(dim) = expectedNeighborMultiIdx(dim) + expectedDiff;

                        testCase.verifyEqual(neighborMultiIdx, expectedNeighborMultiIdx);
                    end
                end
            end
        end

        function testGraphify(testCase, resolution)
            nDims = length(resolution);
            bbox = repmat([0, 1], 1, nDims);

            grid = approx.mesh.UniformGrid(bbox, resolution);
            graph = grid.graphify();

            testCase.verifyEqual(graph.nDims, nDims);

            nVertices = graph.nVertices;
            for d = 1:nDims
                testCase.verifyEqual(nVertices(d), resolution(d)+1);
            end

            nEdges = graph.nEdges;
            for d = 1:nDims
                testCase.verifyEqual(nEdges(d), resolution(d));
            end
        end
        
        function testElementJacobianDeterminants(testCase, resolution)
            nDims = length(resolution);
            bbox = repmat([0, 1], 1, nDims);

            grid = approx.mesh.UniformGrid(bbox, resolution);

            detJ = grid.computeElementJacobianDeterminants();
            expectedDeterminant = prod(1 ./ resolution);

            %< For uniform grid, should return scalar
            testCase.verifyTrue(isscalar(detJ));
            testCase.verifyEqual(detJ, expectedDeterminant, 'AbsTol', 1e-10);
            testCase.verifyTrue(detJ > 0);
        end
        
        function testElementInverseJacobians(testCase, resolution)
            nDims = length(resolution);
            bbox = repmat([0, 1], 1, nDims);

            grid = approx.mesh.UniformGrid(bbox, resolution);

            invJ = grid.computeElementInverseJacobians();
            
            %< For uniform grid, should return nDims x nDims matrix
            testCase.verifyEqual(size(invJ), [nDims, nDims]);
            
            %< Should be diagonal matrix
            testCase.verifyTrue(isdiag(invJ));
            
            %< Check diagonal entries
            expectedDiag = resolution;
            actualDiag = diag(invJ);
            testCase.verifyEqual(actualDiag(:), expectedDiag(:), 'AbsTol', 1e-10);
            
            %< All diagonal entries should be positive
            testCase.verifyTrue(all(diag(invJ) > 0));
        end
        
        function testFaceJacobianDeterminants(testCase, resolution)
            nDims = length(resolution);
            bbox = repmat([0, 1], 1, nDims);

            grid = approx.mesh.UniformGrid(bbox, resolution);
            
            for faceIndex = 1:2*nDims
                detJFace = grid.computeFaceJacobianDeterminants(faceIndex);
                
                %< For uniform grid, should return scalar
                testCase.verifyTrue(isscalar(detJFace));
                testCase.verifyTrue(detJFace > 0);
                
                %< Compute expected value
                dim = ceil(faceIndex / 2);
                spacings = 1 ./ resolution;
                otherDims = 1:nDims;
                otherDims(dim) = [];
                
                if isempty(otherDims)
                    %< 1D case: face is a point
                    expectedDetJFace = 1;
                else
                    expectedDetJFace = prod(spacings(otherDims));
                end
                
                testCase.verifyEqual(detJFace, expectedDetJFace, 'AbsTol', 1e-10);
            end
        end
        
        function testOutwardNormals(testCase, resolution)
            nDims = length(resolution);
            bbox = repmat([0, 1], 1, nDims);

            grid = approx.mesh.UniformGrid(bbox, resolution);
            
            for faceIndex = 1:2*nDims
                normals = grid.computeOutwardNormals(faceIndex);
                
                %< Should return nDims x 1 vector
                testCase.verifyEqual(size(normals), [nDims, 1]);
                
                %< Should be unit vector
                testCase.verifyEqual(norm(normals), 1, 'AbsTol', 1e-10);
                
                %< Should be axis-aligned
                dim = ceil(faceIndex / 2);
                isPositive = (mod(faceIndex, 2) == 0);
                
                expectedNormal = zeros(nDims, 1);
                if isPositive
                    expectedNormal(dim) = 1;
                else
                    expectedNormal(dim) = -1;
                end
                
                testCase.verifyEqual(normals, expectedNormal, 'AbsTol', 1e-10);
            end
        end
        
        function testInvalidFaceIndex(testCase)
            bbox = [0, 1];
            resolution = [5];
            grid = approx.mesh.UniformGrid(bbox, resolution);
            
            %< Test invalid face indices
            testCase.verifyError(...
                @() grid.computeFaceJacobianDeterminants(0), ...
                'approx:mesh:UniformGrid:InvalidInput');
            testCase.verifyError(...
                @() grid.computeFaceJacobianDeterminants(3), ...
                'approx:mesh:UniformGrid:InvalidInput');
            testCase.verifyError(...
                @() grid.computeOutwardNormals(0), ...
                'approx:mesh:UniformGrid:InvalidInput');
            testCase.verifyError(...
                @() grid.computeOutwardNormals(3), ...
                'approx:mesh:UniformGrid:InvalidInput');
        end

        function testCollocate(testCase, resolution)
            nDims = length(resolution);
            bbox = repmat([0, 1], 1, nDims);

            grid = approx.mesh.UniformGrid(bbox, resolution);

            %< Test collocation at element centers (reference coordinate = 0)
            X = cell(1, nDims);
            for d = 1:nDims
                X{d} = 0;
            end

            Y = grid.collocate(X);

            testCase.verifyTrue(iscell(Y));
            testCase.verifyEqual(length(Y), nDims);

            %< Test collocation for specific element
            if all(resolution >= 2)
                testElement = grid.multiToLinear(repmat(2, 1, nDims));
                Y2 = grid.collocate(X, testElement);

                testCase.verifyTrue(ismatrix(Y2));
                testCase.verifyEqual(size(Y2, 2), 1);
                testCase.verifyEqual(size(Y2, 1), nDims);

                for d = 1:nDims
                    expectedCoord = grid.centroids{d}(2);
                    testCase.verifyEqual(Y2(d), expectedCoord, 'AbsTol', 1e-10);
                end
            end
        end

        function testCollocateEdgeCases(testCase)
            bbox = [0, 2, 0, 3];
            n = [2, 3];

            grid = approx.mesh.UniformGrid(bbox, n);

            X = {[0.5, -0.5], [0, 0.5]};
            testMultiIndices = [1, 1; 2, 3];
            testLinearIndices = grid.multiToLinear(testMultiIndices);
            Y = grid.collocate(X, testLinearIndices);

            testCase.verifyEqual(size(Y), [2, 4]);

            expectedY = zeros(2, 4);
            for i = 1:2
                expectedY(1, 2*i-1:2*i) = grid.centroids{1}(testMultiIndices(i, 1)) + grid.spacings{1} * X{1};
                expectedY(2, 2*i-1:2*i) = grid.centroids{2}(testMultiIndices(i, 2)) + grid.spacings{2} * X{2};
            end

            testCase.verifyEqual(Y, expectedY, 'AbsTol', 1e-10);
        end

        function testDifferentBoundingBoxes(testCase)
            %< Test with non-unit bounding boxes
            bbox = [-1, 2, 0.5, 3.5];
            n = [4, 3];

            grid = approx.mesh.UniformGrid(bbox, n);

            testCase.verifyEqual(grid.nDims, 2);
            testCase.verifyEqual(grid.resolution, n);

            %< Check spacing calculations
            expectedSpacing1 = (2 - (-1)) / 4; %< 0.75
            expectedSpacing2 = (3.5 - 0.5) / 3; %< 1.0

            testCase.verifyEqual(grid.spacings{1}, expectedSpacing1, 'AbsTol', 1e-10);
            testCase.verifyEqual(grid.spacings{2}, expectedSpacing2, 'AbsTol', 1e-10);

            %< Check centroid calculations
            expectedCentroids1 = -1 + expectedSpacing1 / 2 + (0:3) * expectedSpacing1;
            expectedCentroids2 = 0.5 + expectedSpacing2 / 2 + (0:2) * expectedSpacing2;

            testCase.verifyEqual(grid.centroids{1}, expectedCentroids1, 'AbsTol', 1e-10);
            testCase.verifyEqual(grid.centroids{2}, expectedCentroids2, 'AbsTol', 1e-10);
        end

        function testConstructorWithDifferentDimensions(testCase)
            %< Test 1D case
            grid1D = approx.mesh.UniformGrid([-2, 3], [10]);
            testCase.verifyEqual(grid1D.nDims, 1);
            testCase.verifyEqual(grid1D.resolution, [10]);
            testCase.verifyEqual(grid1D.spacings{1}, 0.5, 'AbsTol', 1e-10);

            %< Test 3D case
            grid3D = approx.mesh.UniformGrid([0, 1, -1, 2, 0.5, 1.5], [5, 4, 6]);
            testCase.verifyEqual(grid3D.nDims, 3);
            testCase.verifyEqual(grid3D.resolution, [5, 4, 6]);
            testCase.verifyEqual(grid3D.spacings{1}, 0.2, 'AbsTol', 1e-10);
            testCase.verifyEqual(grid3D.spacings{2}, 0.75, 'AbsTol', 1e-10);
            testCase.verifyEqual(grid3D.spacings{3}, 1/6, 'AbsTol', 1e-10);
        end

        function testRefine(testCase)
            %< Test refinement for simple 2D case
            bbox = [0, 1, 0, 1];
            n = [2, 2];

            grid = approx.mesh.UniformGrid(bbox, n);

            %< Test no refinement
            refinedGrid0 = grid.refine(0);
            testCase.verifyEqual(refinedGrid0.resolution, n);

            %< Test single level refinement
            refinedGrid1 = grid.refine(1);
            testCase.verifyEqual(refinedGrid1.resolution, 2 * n);
            testCase.verifyEqual(refinedGrid1.nDims, grid.nDims);
            testCase.verifyEqual(refinedGrid1.bbox, grid.bbox);
        end
    end
end