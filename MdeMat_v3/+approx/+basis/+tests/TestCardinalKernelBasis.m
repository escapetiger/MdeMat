classdef TestCardinalKernelBasis < matlab.unittest.TestCase
    % TESTCARDINALKERNALBASIS Unit tests for CardinalKernelBasis class.
    %
    %   This test suite verifies the functionality of the manual
    %   CardinalKernelBasis implementation including construction,
    %   barycentric weight computation, Lagrange interpolation properties,
    %   and partition of unity.

    properties
        uniformNodes1D     % Uniform nodes on [-1, 1]
        chebyshevNodes1D   % Chebyshev nodes on [-1, 1]
        testPoints1D       % Evaluation points for 1D tests
        tolerance          % Numerical tolerance for comparisons
    end

    methods (TestMethodSetup)
        function setupTestCase(testCase)
            testCase.tolerance = 1e-12;
            
            testCase.uniformNodes1D = {
                0, ...
                [-1, 1], ...
                [-1, 0, 1], ...
                [-1, -1/3, 1/3, 1], ...
                [-1, -1/2, 0, 1/2, 1]
            };
            
            testCase.chebyshevNodes1D = cell(1, 5);
            for n = 1:5
                k = 0:(n-1);
                testCase.chebyshevNodes1D{n} = -cos(pi * (2*k + 1) / (2*n));
            end
            
            testCase.testPoints1D = linspace(-1.5, 1.5, 21);
        end
    end

    methods (Test)
        function testConstructorBasic(testCase)
            nodes = [-1, 0, 1];
            basis = approx.basis.CardinalKernelBasis(nodes);
            
            testCase.verifyEqual(basis.nDims, 1);
            testCase.verifyEqual(basis.nCodims, 3);
            testCase.verifyEqual(size(basis.nodes), [1, 3]);
            testCase.verifyEqual(basis.nodes, nodes);
            testCase.verifyTrue(basis.isWellDefined);
            
            testCase.verifyEqual(size(basis.weights), [3, 1]);
            testCase.verifyTrue(all(isfinite(basis.weights)));
            testCase.verifyTrue(all(basis.weights ~= 0));
        end

        function testConstructorDifferentSizes(testCase)
            for nNodes = 1:5
                nodes = testCase.uniformNodes1D{nNodes};
                basis = approx.basis.CardinalKernelBasis(nodes);
                
                testCase.verifyEqual(basis.nDims, 1);
                testCase.verifyEqual(basis.nCodims, nNodes);
                testCase.verifyEqual(size(basis.nodes), [1, nNodes]);
                testCase.verifyEqual(size(basis.weights), [nNodes, 1]);
                testCase.verifyTrue(basis.isWellDefined);
            end
        end

        function testConstructorMultiDimensional(testCase)
            nodes2D = [0, 1, 0; 0, 0, 1];
            basis2D = approx.basis.CardinalKernelBasis(nodes2D);
            
            testCase.verifyEqual(basis2D.nDims, 2);
            testCase.verifyEqual(basis2D.nCodims, 3);
            testCase.verifyEqual(size(basis2D.nodes), [2, 3]);
            testCase.verifyEqual(size(basis2D.weights), [3, 1]);
            
            nodes3D = [0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
            basis3D = approx.basis.CardinalKernelBasis(nodes3D);
            
            testCase.verifyEqual(basis3D.nDims, 3);
            testCase.verifyEqual(basis3D.nCodims, 4);
            testCase.verifyEqual(size(basis3D.nodes), [3, 4]);
            testCase.verifyEqual(size(basis3D.weights), [4, 1]);
        end

        function testBaryCentricWeightsSymmetric(testCase)
            nodes = [-1, 0, 1];
            basis = approx.basis.CardinalKernelBasis(nodes);
            
            testCase.verifyEqual(basis.weights(1), basis.weights(3), ...
                'AbsTol', testCase.tolerance, ...
                'Symmetric nodes should have equal outer weights');
            
            nodes5 = [-1, -0.5, 0, 0.5, 1];
            basis5 = approx.basis.CardinalKernelBasis(nodes5);
            
            testCase.verifyEqual(basis5.weights(1), basis5.weights(5), ...
                'AbsTol', testCase.tolerance);
            testCase.verifyEqual(basis5.weights(2), basis5.weights(4), ...
                'AbsTol', testCase.tolerance);
        end

        function testLagrangeInterpolationProperty(testCase)
            for nNodes = 2:5
                nodes = testCase.uniformNodes1D{nNodes};
                basis = approx.basis.CardinalKernelBasis(nodes);
                
                Y = basis.evaluate(nodes);
                
                for i = 1:nNodes
                    for j = 1:nNodes
                        if i == j
                            testCase.verifyEqual(Y(i, j), 1, ...
                                'AbsTol', testCase.tolerance, ...
                                sprintf('L_%d(x_%d) should be 1', i, j));
                        else
                            testCase.verifyEqual(Y(i, j), 0, ...
                                'AbsTol', testCase.tolerance, ...
                                sprintf('L_%d(x_%d) should be 0', i, j));
                        end
                    end
                end
            end
        end

        function testLagrangeInterpolationPropertyChebyshev(testCase)
            for nNodes = 2:5
                nodes = testCase.chebyshevNodes1D{nNodes};
                basis = approx.basis.CardinalKernelBasis(nodes);
                
                Y = basis.evaluate(nodes);
                
                for i = 1:nNodes
                    for j = 1:nNodes
                        if i == j
                            testCase.verifyEqual(Y(i, j), 1, ...
                                'AbsTol', testCase.tolerance, ...
                                sprintf('Chebyshev L_%d(x_%d) should be 1', i, j));
                        else
                            testCase.verifyEqual(Y(i, j), 0, ...
                                'AbsTol', testCase.tolerance, ...
                                sprintf('Chebyshev L_%d(x_%d) should be 0', i, j));
                        end
                    end
                end
            end
        end

        function testPartitionOfUnity(testCase)
            for nNodes = 2:5
                nodes = testCase.uniformNodes1D{nNodes};
                basis = approx.basis.CardinalKernelBasis(nodes);
                
                x = testCase.testPoints1D;
                Y = basis.evaluate(x);
                sumY = sum(Y, 1);
                
                testCase.verifyEqual(sumY, ones(size(x)), ...
                    'AbsTol', testCase.tolerance, ...
                    sprintf('Partition of unity failed for %d nodes', nNodes));
            end
        end

        function testPartitionOfUnityChebyshev(testCase)
            for nNodes = 2:5
                nodes = testCase.chebyshevNodes1D{nNodes};
                basis = approx.basis.CardinalKernelBasis(nodes);
                
                x = testCase.testPoints1D;
                Y = basis.evaluate(x);
                sumY = sum(Y, 1);
                
                testCase.verifyEqual(sumY, ones(size(x)), ...
                    'AbsTol', testCase.tolerance, ...
                    sprintf('Chebyshev partition of unity failed for %d nodes', nNodes));
            end
        end

        function testEvaluationAtBoundaryPoints(testCase)
            nodes = [-1, 0, 1];
            basis = approx.basis.CardinalKernelBasis(nodes);
            
            Y_left = basis.evaluate(-1);
            testCase.verifyEqual(Y_left(1), 1, 'AbsTol', testCase.tolerance);
            testCase.verifyEqual(Y_left(2), 0, 'AbsTol', testCase.tolerance);
            testCase.verifyEqual(Y_left(3), 0, 'AbsTol', testCase.tolerance);
            
            Y_center = basis.evaluate(0);
            testCase.verifyEqual(Y_center(1), 0, 'AbsTol', testCase.tolerance);
            testCase.verifyEqual(Y_center(2), 1, 'AbsTol', testCase.tolerance);
            testCase.verifyEqual(Y_center(3), 0, 'AbsTol', testCase.tolerance);
            
            Y_right = basis.evaluate(1);
            testCase.verifyEqual(Y_right(1), 0, 'AbsTol', testCase.tolerance);
            testCase.verifyEqual(Y_right(2), 0, 'AbsTol', testCase.tolerance);
            testCase.verifyEqual(Y_right(3), 1, 'AbsTol', testCase.tolerance);
        end

        function testEvaluationMultiplePoints(testCase)
            nodes = [-1, 0, 1];
            basis = approx.basis.CardinalKernelBasis(nodes);
            
            x = [-0.5, 0.25, 0.75];
            Y = basis.evaluate(x);
            
            testCase.verifyEqual(size(Y), [3, 3]);
            
            testCase.verifyEqual(sum(Y, 1), ones(1, 3), ...
                'AbsTol', testCase.tolerance);
            
            for i = 1:length(x)
                Y_single = basis.evaluate(x(i));
                testCase.verifyEqual(Y(:, i), Y_single, ...
                    'AbsTol', testCase.tolerance);
            end
        end

        function testNumericalStabilityNearNodes(testCase)
            nodes = [-1, 0, 1];
            basis = approx.basis.CardinalKernelBasis(nodes);
            
            eps_vals = [1e-10, 1e-12, 1e-14];
            
            for eps_val = eps_vals
                x_near = nodes + eps_val;
                Y = basis.evaluate(x_near);
                
                sumY = sum(Y, 1);
                testCase.verifyEqual(sumY, ones(size(x_near)), ...
                    'AbsTol', 1e-10, ...
                    sprintf('Stability failed for eps = %.0e', eps_val));
                
                for i = 1:length(nodes)
                    testCase.verifyGreaterThan(Y(i, i), 0.99, ...
                        sprintf('Near-node stability failed for node %d', i));
                end
            end
        end

        function testExactNodeEvaluation(testCase)
            nodes = testCase.chebyshevNodes1D{4};
            basis = approx.basis.CardinalKernelBasis(nodes);
            
            x_mixed = [nodes(1), 0.5, nodes(3), -0.3, nodes(4)];
            Y = basis.evaluate(x_mixed);
            
            testCase.verifyEqual(Y(1, 1), 1, 'AbsTol', testCase.tolerance);
            testCase.verifyEqual(Y(3, 3), 1, 'AbsTol', testCase.tolerance);
            testCase.verifyEqual(Y(4, 5), 1, 'AbsTol', testCase.tolerance);
            
            testCase.verifyEqual(Y(2, 1), 0, 'AbsTol', testCase.tolerance);
            testCase.verifyEqual(Y(3, 1), 0, 'AbsTol', testCase.tolerance);
            testCase.verifyEqual(Y(4, 1), 0, 'AbsTol', testCase.tolerance);
        end

        function testMultiDimensionalEvaluation(testCase)
            nodes2D = [0, 1, 0.5; 0, 0, sqrt(3)/2];
            basis2D = approx.basis.CardinalKernelBasis(nodes2D);
            
            Y_vertices = basis2D.evaluate(nodes2D);
            expected = eye(3);
            testCase.verifyEqual(Y_vertices, expected, ...
                'AbsTol', testCase.tolerance, ...
                '2D Lagrange property failed');
            
            centroid = mean(nodes2D, 2);
            Y_centroid = basis2D.evaluate(centroid);
            
            testCase.verifyEqual(Y_centroid(1), Y_centroid(2), ...
                'AbsTol', 1e-10);
            testCase.verifyEqual(Y_centroid(2), Y_centroid(3), ...
                'AbsTol', 1e-10);
            
            testCase.verifyEqual(sum(Y_centroid), 1, ...
                'AbsTol', testCase.tolerance);
        end

        function testInputValidation(testCase)
            nodes = [-1, 0, 1];
            basis = approx.basis.CardinalKernelBasis(nodes);
            
            testCase.verifyError(...
                @() basis.evaluate([0.5; 0.5]), ...
                'core:function:Function:InvalidInput');
            
            testCase.verifyError(...
                @() basis.evaluate([]), ...
                'core:function:Function:InvalidInput');
            
            testCase.verifyError(...
                @() basis.gradient(0.5), ...
                'core:function:Function:NotScalarValued');
        end

        function testDifferentiationNotImplemented(testCase)
            nodes = [-1, 0, 1];
            basis = approx.basis.CardinalKernelBasis(nodes);
            
            testCase.verifyError(...
                @() basis.differentiate(0.5, [1]), ...
                'approx:basis:CardinalKernelBasis:NotImplemented');
            
            try
                basis.differentiate(0.5, [1]);
                testCase.verifyFail('Expected NotImplemented error');
            catch ME
                testCase.verifyTrue(contains(ME.message, 'not implemented'), ...
                    'Error message should indicate not implemented');
            end
        end

        function testJacobianFails(testCase)
            nodes = [-1, 0, 1];
            basis = approx.basis.CardinalKernelBasis(nodes);
            
            testCase.verifyError(...
                @() basis.jacobian(0.5), ...
                'approx:basis:CardinalKernelBasis:NotImplemented');
        end

        function testHighDegreePolynomial(testCase)
            n = 10;
            nodes = linspace(-1, 1, n);
            basis = approx.basis.CardinalKernelBasis(nodes);
            
            testCase.verifyEqual(basis.nCodims, n);
            
            x = linspace(-0.9, 0.9, 11);
            Y = basis.evaluate(x);
            sumY = sum(Y, 1);
            
            testCase.verifyEqual(sumY, ones(size(x)), ...
                'AbsTol', 1e-10, ...
                'High-degree partition of unity failed');
            
            test_indices = [1, n/2, n];
            for i = test_indices
                Y_node = basis.evaluate(nodes(i));
                testCase.verifyEqual(Y_node(i), 1, ...
                    'AbsTol', testCase.tolerance, ...
                    sprintf('High-degree interpolation failed at node %d', i));
            end
        end

        function testRobustnessRandomNodes(testCase)
            rng(42);
            
            for nNodes = [3, 5, 7]
                nodes = sort(-2 + 4*rand(1, nNodes));
                basis = approx.basis.CardinalKernelBasis(nodes);
                
                testCase.verifyTrue(basis.isWellDefined);
                testCase.verifyTrue(all(isfinite(basis.weights)));
                
                Y_nodes = basis.evaluate(nodes);
                expected_identity = eye(nNodes);
                testCase.verifyEqual(Y_nodes, expected_identity, ...
                    'AbsTol', testCase.tolerance, ...
                    sprintf('Random nodes interpolation failed for n=%d', nNodes));
                
                x_test = linspace(min(nodes)-0.5, max(nodes)+0.5, 15);
                Y_test = basis.evaluate(x_test);
                sumY_test = sum(Y_test, 1);
                testCase.verifyEqual(sumY_test, ones(size(x_test)), ...
                    'AbsTol', testCase.tolerance, ...
                    sprintf('Random nodes partition failed for n=%d', nNodes));
            end
        end

        function testMemoryEfficiency(testCase)
            nodes = linspace(-1, 1, 20);
            basis = approx.basis.CardinalKernelBasis(nodes);
            
            nPoints = 1000;
            x_large = linspace(-1.2, 1.2, nPoints);
            
            try
                Y_large = basis.evaluate(x_large);
                
                testCase.verifyEqual(size(Y_large), [20, nPoints]);
                
                indices = [1, nPoints/2, nPoints];
                for idx = indices
                    testCase.verifyEqual(sum(Y_large(:, idx)), 1, ...
                        'AbsTol', testCase.tolerance, ...
                        sprintf('Large evaluation partition failed at point %d', idx));
                end
                
            catch ME
                if contains(ME.identifier, 'OutOfMemory')
                    testCase.assumeFail('Insufficient memory for large evaluation test');
                else
                    rethrow(ME);
                end
            end
        end
    end
end