classdef TestGrid < matlab.unittest.TestCase
    
    properties (TestParameter)
        GridType = struct(...
            'linear1D',     {{linspace(0, 1, 6)}}, ...
            'nonlinear1D',  {{[0, 0.1, 0.3, 0.6, 1.0]}}, ...
            'linear2D',     {{linspace(0, 1, 5), linspace(0, 1, 4)}}, ...
            'nonlinear2D',  {{[0, 0.2, 0.5, 0.9, 1], [0, 0.3, 0.7, 1]}}, ...
            'linear3D',     {{linspace(0, 1, 4), linspace(0, 1, 3), linspace(0, 1, 5)}}, ...
            'nonlinear3D',  {{[0, 0.25, 0.75, 1], [0, 0.4, 1], [0, 0.2, 0.6, 0.9, 1]}})
    end
    
    methods (Test)
        function testConstructor(testCase, GridType)
            grid = approx.mesh.Grid(GridType);
            
            nDims = length(GridType);
            resolution = zeros(1, nDims);
            for i = 1:nDims
                resolution(i) = length(GridType{i}) - 1;
            end
            
            testCase.verifyEqual(grid.NDims, nDims);
            testCase.verifyEqual(grid.Resolution, resolution);
            testCase.verifyEqual(length(grid.Centroids), nDims);
            testCase.verifyEqual(length(grid.Spacings), nDims);
            testCase.verifyEqual(length(grid.Elements), nDims);
            testCase.verifyEqual(length(grid.Boundary), 2*nDims);
            
            % Verify centroids computation
            for i = 1:nDims
                vertices = GridType{i};
                expectedCentroids = (vertices(1:end-1) + vertices(2:end)) / 2;
                testCase.verifyEqual(grid.Centroids{i}, expectedCentroids, 'AbsTol', 1e-10);
            end
            
            % Verify spacings computation
            for i = 1:nDims
                vertices = GridType{i};
                expectedSpacings = diff(vertices);
                testCase.verifyEqual(grid.Spacings{i}, expectedSpacings, 'AbsTol', 1e-10);
            end
        end
        
        function testMeasure(testCase, GridType)
            grid = approx.mesh.Grid(GridType);
            
            % Compute expected minimum spacing across all dimensions
            minSpacings = zeros(1, length(GridType));
            for i = 1:length(GridType)
                minSpacings(i) = min(diff(GridType{i}));
            end
            expectedMeasure = min(minSpacings);
            
            testCase.verifyEqual(grid.computeMeasure(), expectedMeasure, 'AbsTol', 1e-10);
        end
        
        function testNElements(testCase, GridType)
            grid = approx.mesh.Grid(GridType);
            
            resolution = zeros(1, length(GridType));
            for i = 1:length(GridType)
                resolution(i) = length(GridType{i}) - 1;
            end
            expectedTotal = prod(resolution);
            
            testCase.verifyEqual(grid.NElements, expectedTotal);
        end
        
        function testLinearToMultiIndices(testCase, GridType)
            grid = approx.mesh.Grid(GridType);
            
            % Test conversion for a few linear indices
            nElements = grid.NElements;
            testIndices = [1, min(5, nElements), nElements];
            
            for linearIdx = testIndices
                multiIdx = grid.Indexer.linearToMulti(linearIdx);
                backToLinear = grid.Indexer.multiToLinear(multiIdx);
                
                testCase.verifyEqual(backToLinear, linearIdx);
                testCase.verifyEqual(size(multiIdx, 1), 1);
                testCase.verifyEqual(size(multiIdx, 2), grid.NDims);
                
                % Verify indices are within bounds
                for d = 1:grid.NDims
                    testCase.verifyTrue(multiIdx(d) >= 1 && multiIdx(d) <= grid.Resolution(d));
                end
            end
        end
        
        function testGetInternalElements(testCase, GridType)
            grid = approx.mesh.Grid(GridType);
            
            % Skip test for grids with resolution <= 2 in any dimension
            if any(grid.Resolution <= 2)
                return;
            end
            
            internalElements = grid.getInternalElements();
            internalMultiIndices = grid.Indexer.linearToMulti(internalElements);
            
            % Verify internal elements don't touch boundaries
            for d = 1:grid.NDims
                testCase.verifyTrue(all(internalMultiIndices(:, d) > 1));
                testCase.verifyTrue(all(internalMultiIndices(:, d) < grid.Resolution(d)));
            end
            
            % Verify expected count
            expectedCount = prod(max(0, grid.Resolution - 2));
            testCase.verifyEqual(length(internalElements), expectedCount);
        end
        
        function testGetBoundaryElements(testCase, GridType)
            grid = approx.mesh.Grid(GridType);
            
            for iBoundary = 1:2*grid.NDims
                boundaryElements = grid.getBoundaryElements(iBoundary);
                boundaryMultiIndices = grid.Indexer.linearToMulti(boundaryElements);
                
                dim = ceil(iBoundary/2);
                isLower = mod(iBoundary, 2) == 1;
                
                expectedCount = prod(grid.Resolution) / grid.Resolution(dim);
                testCase.verifyEqual(length(boundaryElements), expectedCount);
                
                if isLower
                    testCase.verifyTrue(all(boundaryMultiIndices(:, dim) == 1));
                else
                    testCase.verifyTrue(all(boundaryMultiIndices(:, dim) == grid.Resolution(dim)));
                end
            end
        end
        
        function testFindNeighborElements(testCase, GridType)
            grid = approx.mesh.Grid(GridType);
            
            % Test neighbor finding for element in the middle (if exists)
            if all(grid.Resolution >= 3)
                centerMultiIdx = ceil(grid.Resolution / 2);
                centerLinearIdx = grid.Indexer.multiToLinear(centerMultiIdx);
                
                for iBoundary = 1:2*grid.NDims
                    % Use periodic boundary condition for internal elements
                    [neighborLinearIdx, ~] = grid.findNeighborElements(centerLinearIdx, iBoundary, 'periodic');
                    
                    if ~isempty(neighborLinearIdx) && neighborLinearIdx > 0
                        neighborMultiIdx = grid.Indexer.linearToMulti(neighborLinearIdx);
                        
                        dim = ceil(iBoundary/2);
                        expectedDiff = -(mod(iBoundary, 2) * 2 - 1); % -1 for odd, +1 for even
                        
                        expectedNeighborMultiIdx = centerMultiIdx;
                        expectedNeighborMultiIdx(dim) = expectedNeighborMultiIdx(dim) + expectedDiff;
                        
                        testCase.verifyEqual(neighborMultiIdx, expectedNeighborMultiIdx);
                    end
                end
            end
        end
        
        function testGraphify(testCase, GridType)
            grid = approx.mesh.Grid(GridType);
            
            graph = grid.graphify();
            
            testCase.verifyEqual(graph.NDims, grid.NDims);
                
            nVertices = graph.NVertices;
            for d = 1:grid.NDims
                testCase.verifyEqual(nVertices(d), grid.Resolution(d) + 1);
            end
            
            nEdges = graph.NEdges;
            for d = 1:grid.NDims
                testCase.verifyEqual(nEdges(d), grid.Resolution(d));
            end
        end
        
        function testElementJacobianDeterminants(testCase, GridType)
            grid = approx.mesh.Grid(GridType);
            
            detJ = grid.computeElementJacobianDeterminants();
            
            % Should return vector of length nElements
            testCase.verifyTrue(isnumeric(detJ));
            testCase.verifyEqual(length(detJ), grid.NElements);
            testCase.verifyTrue(all(detJ > 0));
            
            % For uniform spacing in each dimension, verify calculation
            if all(cellfun(@(x) length(unique(diff(x))) == 1, GridType))
                expectedDeterminant = prod(cellfun(@(x) diff(x(1:2)), GridType));
                testCase.verifyEqual(detJ, expectedDeterminant * ones(size(detJ)), ...
                    'AbsTol', 1e-10);
            end
        end
        
        function testElementInverseJacobians(testCase, GridType)
            grid = approx.mesh.Grid(GridType);
            
            invJ = grid.computeElementInverseJacobians();
            
            % Should return nDims x nDims x nElements array
            expectedSize = [grid.NDims, grid.NDims, grid.NElements];
            testCase.verifyEqual(size(invJ), expectedSize);
            
            % Each matrix should be diagonal
            for iElement = 1:grid.NElements
                testCase.verifyTrue(isdiag(invJ(:, :, iElement)));
                
                % Diagonal entries should be positive
                testCase.verifyTrue(all(diag(invJ(:, :, iElement)) > 0));
            end
            
            % Verify against expected values for specific elements
            if grid.NElements >= 1
                % Test first element
                M = grid.Indexer.linearToMulti(1);
                expectedInvJ = zeros(grid.NDims);
                for d = 1:grid.NDims
                    expectedInvJ(d, d) = 1 / grid.Spacings{d}(M(d));
                end
                testCase.verifyEqual(invJ(:, :, 1), expectedInvJ, 'AbsTol', 1e-10);
            end
        end
        
        function testFaceJacobianDeterminants(testCase, GridType)
            grid = approx.mesh.Grid(GridType);
            
            for faceIndex = 1:2*grid.NDims
                detJFace = grid.computeFaceJacobianDeterminants([], faceIndex);
                
                % Should return vector of length nElements
                testCase.verifyEqual(size(detJFace), [grid.NElements, 1]);
                testCase.verifyTrue(all(detJFace > 0));
                
                % Verify against expected calculation
                dim = ceil(faceIndex / 2);
                otherDims = 1:grid.NDims;
                otherDims(dim) = [];
                
                if isempty(otherDims)
                    % 1D case: should all be 1
                    testCase.verifyEqual(detJFace, ones(grid.NElements, 1), 'AbsTol', 1e-10);
                else
                    % Check a few elements manually
                    for testElement = [1, min(3, grid.NElements)]
                        M = grid.Indexer.linearToMulti(testElement);
                        expectedDetJ = 1;
                        for d = otherDims
                            expectedDetJ = expectedDetJ * grid.Spacings{d}(M(d));
                        end
                        testCase.verifyEqual(detJFace(testElement), expectedDetJ, 'AbsTol', 1e-10);
                    end
                end
            end
        end
        
        function testOutwardNormals(testCase, GridType)
            grid = approx.mesh.Grid(GridType);
            
            for faceIndex = 1:2*grid.NDims
                normals = grid.computeOutwardNormals([], faceIndex);
                
                % Should return nDims x 1 vector
                testCase.verifyEqual(size(normals), [grid.NDims, 1]);
                
                % Should be unit vector
                testCase.verifyEqual(norm(normals), 1, 'AbsTol', 1e-10);
                
                % Should be axis-aligned
                dim = ceil(faceIndex / 2);
                isPositive = (mod(faceIndex, 2) == 0);
                
                expectedNormal = zeros(grid.NDims, 1);
                if isPositive
                    expectedNormal(dim) = 1;
                else
                    expectedNormal(dim) = -1;
                end
                
                testCase.verifyEqual(normals, expectedNormal, 'AbsTol', 1e-10);
            end
        end

        function testCollocate(testCase, GridType)
            grid = approx.mesh.Grid(GridType);

            % Test collocation at element centers (reference coordinate = 0)
            X = cell(1, grid.NDims);
            for d = 1:grid.NDims
                X{d} = 0;
            end

            Y = grid.collocate(X);

            testCase.verifyTrue(ismatrix(Y));
            testCase.verifyEqual(size(Y, 1), grid.NDims);

            % Test collocation for specific element
            if all(grid.Resolution >= 2)
                testElement = grid.Indexer.multiToLinear(repmat(2, 1, grid.NDims));
                Y2 = grid.collocate(X, testElement);

                testCase.verifyTrue(ismatrix(Y2));
                testCase.verifyEqual(size(Y2, 2), 1);
                testCase.verifyEqual(size(Y2, 1), grid.NDims);

                for d = 1:grid.NDims
                    expectedCoord = grid.Centroids{d}(2);
                    testCase.verifyEqual(Y2{d}, expectedCoord, 'AbsTol', 1e-10);
                end
            end
        end

        function testNonuniformCollocate(testCase)
            % Test specific nonuniform case
            gridCoords = {[0, 0.25, 0.75, 1.0]};
            grid = approx.mesh.Grid(gridCoords);

            refNodes = {[-0.5, 0, 0.5]};
            elementIdx = 2;

            Y = grid.collocate(refNodes, elementIdx);

            h = grid.Spacings{1}(elementIdx);
            c = grid.Centroids{1}(elementIdx);
            expectedCoords = c + h * refNodes{1};

            testCase.verifyEqual(Y{1}, expectedCoords, 'AbsTol', 1e-10);
        end

        function testRefine(testCase)
            % Test refinement for simple 1D case
            gridCoords = {[0, 0.5, 1.0]};
            grid = approx.mesh.Grid(gridCoords);

            % Test no refinement
            refinedGrid0 = grid.refine(0);
            testCase.verifyEqual(refinedGrid0.Resolution, grid.Resolution);

            % Test single level refinement
            refinedGrid1 = grid.refine(1);
            testCase.verifyEqual(refinedGrid1.Resolution, 2 * grid.Resolution);
            testCase.verifyEqual(refinedGrid1.NDims, grid.NDims);
        end
        
        function testSpecificNonuniformCase(testCase)
            % Test a specific nonuniform case with known values
            gridCoords = {[0, 0.2, 0.7, 1.0], [0, 0.3, 1.0]};
            grid = approx.mesh.Grid(gridCoords);
            
            % Test element Jacobian determinants
            detJ = grid.computeElementJacobianDeterminants();
            
            % Expected spacings: [0.2, 0.5, 0.3] x [0.3, 0.7]
            expectedSpacings1 = [0.2, 0.5, 0.3];
            expectedSpacings2 = [0.3, 0.7];
            
            % Expected determinants for each element (in linear order)
            expectedDetJ = [];
            for j = 1:2
                for i = 1:3
                    expectedDetJ = [expectedDetJ, expectedSpacings1(i) * expectedSpacings2(j)];
                end
            end
            
            testCase.verifyEqual(detJ(:)', expectedDetJ, 'AbsTol', 1e-10);
            
            % Test face Jacobian determinants for face 1 (perpendicular to dim 1)
            detJFace1 = grid.computeFaceJacobianDeterminants([], 1);
            % Should be spacings in dim 2 for each element
            expectedDetJFace1 = [0.3, 0.3, 0.3, 0.7, 0.7, 0.7].';
            testCase.verifyEqual(detJFace1, expectedDetJFace1, 'AbsTol', 1e-10);
            
            % Test face Jacobian determinants for face 3 (perpendicular to dim 2)  
            detJFace3 = grid.computeFaceJacobianDeterminants([], 3);
            % Should be spacings in dim 1 for each element
            expectedDetJFace3 = [0.2, 0.5, 0.3, 0.2, 0.5, 0.3].';
            testCase.verifyEqual(detJFace3, expectedDetJFace3, 'AbsTol', 1e-10);
        end
    end
end