classdef TestRbfKernelBasis < matlab.unittest.TestCase
    % TESTRBFKERNELBASIS Unit tests for RbfKernelBasis class.
    %
    %   This test suite verifies the functionality of the RbfKernelBasis
    %   implementation including construction, kernel evaluation,
    %   derivative computation, and vectorized operations for different
    %   RBF kernel types.

    properties (TestParameter)
        kernelType = {'gaussian', 'multiquadric', 'inverse_multiquadric'}
        nDims = {1, 2, 3}
        epsilonType = {'scalar', 'vector'}
    end

    properties
        testNodes1D       % 1D test nodes
        testNodes2D       % 2D test nodes
        testNodes3D       % 3D test nodes
        testPoints1D      % 1D evaluation points
        testPoints2D      % 2D evaluation points
        testPoints3D      % 3D evaluation points
        tolerance         % Numerical tolerance
    end

    methods (TestMethodSetup)
        function setupTestCase(testCase)
            testCase.tolerance = 1e-12;
            
            testCase.testNodes1D = [-1, -0.5, 0, 0.5, 1];
            testCase.testNodes2D = [0, 1, 0.5, -0.5; 0, 0, 0.8, 0.3];
            testCase.testNodes3D = [0, 1, 0, -1; 0, 0, 1, 0; 0, 0, 0, 1];
            
            testCase.testPoints1D = linspace(-1.5, 1.5, 11);
            testCase.testPoints2D = [linspace(-1, 1, 5); linspace(-1, 1, 5)];
            testCase.testPoints3D = [linspace(-1, 1, 4); linspace(-1, 1, 4); linspace(-1, 1, 4)];
        end
    end

    methods (Test)
        function testConstructorBasic(testCase)
            nodes = testCase.testNodes2D;
            basis = approx.basis.RbfKernelBasis(nodes);
            
            testCase.verifyEqual(basis.nDims, 2);
            testCase.verifyEqual(basis.nCodims, 4);
            testCase.verifyEqual(size(basis.nodes), [2, 4]);
            testCase.verifyEqual(basis.nodes, nodes);
            testCase.verifyTrue(basis.isWellDefined);
            testCase.verifyEqual(basis.kernel, 'gaussian');
            testCase.verifyEqual(basis.epsilon, 1.0);
        end

        function testConstructorWithKernelType(testCase, kernelType)
            nodes = testCase.testNodes2D;
            basis = approx.basis.RbfKernelBasis(nodes, kernelType);
            
            testCase.verifyEqual(basis.kernel, kernelType);
            testCase.verifyEqual(basis.epsilon, 1.0);
            testCase.verifyEqual(basis.nCodims, size(nodes, 2));
        end

        function testConstructorWithEpsilon(testCase)
            nodes = testCase.testNodes2D;
            epsilon = 2.5;
            basis = approx.basis.RbfKernelBasis(nodes, 'gaussian', epsilon);
            
            testCase.verifyEqual(basis.epsilon, epsilon);
            testCase.verifyEqual(basis.kernel, 'gaussian');
        end

        function testConstructorWithVectorEpsilon(testCase)
            nodes = testCase.testNodes2D;
            epsilon = [1.0, 1.5, 2.0, 0.5];
            basis = approx.basis.RbfKernelBasis(nodes, 'multiquadric', epsilon);
            
            testCase.verifyEqual(basis.epsilon, epsilon(:));
            testCase.verifyEqual(length(basis.epsilon), basis.nCodims);
        end

        function testConstructorDifferentDimensions(testCase, nDims)
            switch nDims
                case 1
                    nodes = testCase.testNodes1D;
                case 2
                    nodes = testCase.testNodes2D;
                case 3
                    nodes = testCase.testNodes3D;
            end
            
            basis = approx.basis.RbfKernelBasis(nodes, 'gaussian', 1.0);
            
            testCase.verifyEqual(basis.nDims, nDims);
            testCase.verifyEqual(basis.nCodims, size(nodes, 2));
            testCase.verifyEqual(size(basis.nodes), size(nodes));
        end

        function testGaussianKernelEvaluation(testCase)
            nodes = [0, 1];
            epsilon = 1.0;
            basis = approx.basis.RbfKernelBasis(nodes, 'gaussian', epsilon);
            
            x = 0;
            Y = basis.evaluate(x);
            
            expected = [exp(0), exp(-1)];
            testCase.verifyEqual(Y, expected', ...
                'AbsTol', testCase.tolerance);
            
            x = 1;
            Y = basis.evaluate(x);
            expected = [exp(-1), exp(0)];
            testCase.verifyEqual(Y, expected', ...
                'AbsTol', testCase.tolerance);
        end

        function testMultiquadricKernelEvaluation(testCase)
            nodes = [0, 1];
            epsilon = 1.0;
            basis = approx.basis.RbfKernelBasis(nodes, 'multiquadric', epsilon);
            
            x = 0;
            Y = basis.evaluate(x);
            
            expected = [sqrt(1), sqrt(2)];
            testCase.verifyEqual(Y, expected', ...
                'AbsTol', testCase.tolerance);
        end

        function testInverseMultiquadricKernelEvaluation(testCase)
            nodes = [0, 1];
            epsilon = 1.0;
            basis = approx.basis.RbfKernelBasis(nodes, 'inverse_multiquadric', epsilon);
            
            x = 0;
            Y = basis.evaluate(x);
            
            expected = [1/sqrt(1), 1/sqrt(2)];
            testCase.verifyEqual(Y, expected', ...
                'AbsTol', testCase.tolerance);
        end

        function testEvaluationAtNodes(testCase, kernelType)
            nodes = testCase.testNodes1D;
            basis = approx.basis.RbfKernelBasis(nodes, kernelType, 1.0);
            
            for i = 1:length(nodes)
                Y = basis.evaluate(nodes(i));
                
                switch kernelType
                    case 'gaussian'
                        expected_val = 1;
                    case 'multiquadric'
                        expected_val = 1;
                    case 'inverse_multiquadric'
                        expected_val = 1;
                end
                
                testCase.verifyEqual(Y(i), expected_val, ...
                    'AbsTol', testCase.tolerance, ...
                    sprintf('%s kernel should equal %g at node %d', ...
                    kernelType, expected_val, i));
            end
        end

        function testEvaluationMultiplePoints(testCase, kernelType)
            nodes = testCase.testNodes2D;
            basis = approx.basis.RbfKernelBasis(nodes, kernelType, 1.5);
            
            X = testCase.testPoints2D;
            Y = basis.evaluate(X);
            
            testCase.verifyEqual(size(Y, 1), basis.nCodims);
            testCase.verifyEqual(size(Y, 2), size(X, 2));
            testCase.verifyTrue(all(isfinite(Y(:))));
            testCase.verifyTrue(all(Y(:) >= 0) || strcmp(kernelType, 'inverse_multiquadric'));
        end

        function testEvaluationWithVectorEpsilon(testCase, kernelType)
            nodes = testCase.testNodes2D;
            epsilon = [0.5, 1.0, 1.5, 2.0];
            basis = approx.basis.RbfKernelBasis(nodes, kernelType, epsilon);
            
            X = testCase.testPoints2D;
            Y = basis.evaluate(X);
            
            testCase.verifyEqual(size(Y), [basis.nCodims, size(X, 2)]);
            testCase.verifyTrue(all(isfinite(Y(:))));
        end

        function testKernelPositivity(testCase)
            nodes = testCase.testNodes2D;
            basis_gauss = approx.basis.RbfKernelBasis(nodes, 'gaussian', 1.0);
            basis_mq = approx.basis.RbfKernelBasis(nodes, 'multiquadric', 1.0);
            basis_imq = approx.basis.RbfKernelBasis(nodes, 'inverse_multiquadric', 1.0);
            
            X = testCase.testPoints2D;
            
            Y_gauss = basis_gauss.evaluate(X);
            testCase.verifyTrue(all(Y_gauss(:) > 0), ...
                'Gaussian kernel should be positive');
            
            Y_mq = basis_mq.evaluate(X);
            testCase.verifyTrue(all(Y_mq(:) > 0), ...
                'Multiquadric kernel should be positive');
            
            Y_imq = basis_imq.evaluate(X);
            testCase.verifyTrue(all(Y_imq(:) > 0), ...
                'Inverse multiquadric kernel should be positive');
        end

        function testGaussianDerivative(testCase)
            nodes = [0, 1];
            epsilon = 1.0;
            basis = approx.basis.RbfKernelBasis(nodes, 'gaussian', epsilon);
            
            x = 0.5;
            dY = basis.differentiate(x, [1]);
            
            r1 = abs(x - nodes(1));
            r2 = abs(x - nodes(2));
            
            expected1 = -2 * epsilon^2 * r1 * exp(-(epsilon * r1)^2) * abs(x - nodes(1)) / r1;
            expected2 = -2 * epsilon^2 * r2 * exp(-(epsilon * r2)^2) * abs(x - nodes(2)) / r2;
            
            testCase.verifyEqual(dY(1), expected1, 'AbsTol', testCase.tolerance);
            testCase.verifyEqual(dY(2), expected2, 'AbsTol', testCase.tolerance);
        end

        function testMultiquadricDerivative(testCase)
            nodes = [0, 1];
            epsilon = 1.0;
            basis = approx.basis.RbfKernelBasis(nodes, 'multiquadric', epsilon);
            
            x = 0.5;
            dY = basis.differentiate(x, [1]);
            
            r1 = abs(x - nodes(1));
            r2 = abs(x - nodes(2));
            
            expected1 = epsilon^2 * r1 / sqrt(1 + (epsilon * r1)^2) * abs(x - nodes(1)) / r1;
            expected2 = epsilon^2 * r2 / sqrt(1 + (epsilon * r2)^2) * abs(x - nodes(2)) / r2;
            
            testCase.verifyEqual(dY(1), expected1, 'AbsTol', testCase.tolerance);
            testCase.verifyEqual(dY(2), expected2, 'AbsTol', testCase.tolerance);
        end

        function testInverseMultiquadricDerivative(testCase)
            nodes = [0, 1];
            epsilon = 1.0;
            basis = approx.basis.RbfKernelBasis(nodes, 'inverse_multiquadric', epsilon);
            
            x = 0.5;
            dY = basis.differentiate(x, [1]);
            
            r1 = abs(x - nodes(1));
            r2 = abs(x - nodes(2));
            
            denom1 = (1 + (epsilon * r1)^2)^(3/2);
            denom2 = (1 + (epsilon * r2)^2)^(3/2);
            
            expected1 = -epsilon^2 * r1 / denom1 * abs(x - nodes(1)) / r1;
            expected2 = -epsilon^2 * r2 / denom2 * abs(x - nodes(2)) / r2;
            
            testCase.verifyEqual(dY(1), expected1, 'AbsTol', testCase.tolerance);
            testCase.verifyEqual(dY(2), expected2, 'AbsTol', testCase.tolerance);
        end

        function testDerivativeAtNodes(testCase, kernelType)
            nodes = testCase.testNodes1D;
            basis = approx.basis.RbfKernelBasis(nodes, kernelType, 1.0);
            
            for i = 1:length(nodes)
                dY = basis.differentiate(nodes(i), [1]);
                
                testCase.verifyEqual(dY(i), 0, 'AbsTol', testCase.tolerance, ...
                    sprintf('Derivative should be 0 at node %d for %s kernel', ...
                    i, kernelType));
            end
        end

        function testMultiDimensionalDerivative(testCase, kernelType)
            nodes = testCase.testNodes2D;
            basis = approx.basis.RbfKernelBasis(nodes, kernelType, 1.0);
            
            X = [0.3; 0.7];
            
            dY_dx1 = basis.differentiate(X, [1, 0]);
            dY_dx2 = basis.differentiate(X, [0, 1]);
            
            testCase.verifyEqual(size(dY_dx1), [basis.nCodims, 1]);
            testCase.verifyEqual(size(dY_dx2), [basis.nCodims, 1]);
            testCase.verifyTrue(all(isfinite(dY_dx1)));
            testCase.verifyTrue(all(isfinite(dY_dx2)));
        end

        function testJacobianComputation(testCase, kernelType)
            nodes = testCase.testNodes2D;
            basis = approx.basis.RbfKernelBasis(nodes, kernelType, 1.5);
            
            X = [0.2; 0.8];
            J = basis.jacobian(X);
            
            testCase.verifyEqual(size(J), [basis.nCodims, basis.nDims]);
            
            dY_dx1 = basis.differentiate(X, [1, 0]);
            dY_dx2 = basis.differentiate(X, [0, 1]);
            
            for i = 1:basis.nCodims
                testCase.verifyEqual(J(i, 1, 1), dY_dx1(i), ...
                    'AbsTol', testCase.tolerance);
                testCase.verifyEqual(J(i, 2, 1), dY_dx2(i), ...
                    'AbsTol', testCase.tolerance);
            end
        end

        function testVectorizedEvaluation(testCase, kernelType)
            nodes = testCase.testNodes2D;
            basis = approx.basis.RbfKernelBasis(nodes, kernelType, 1.0);
            
            X = testCase.testPoints2D;
            Y_vectorized = basis.evaluate(X);
            
            Y_individual = zeros(basis.nCodims, size(X, 2));
            for i = 1:size(X, 2)
                Y_individual(:, i) = basis.evaluate(X(:, i));
            end
            
            testCase.verifyEqual(Y_vectorized, Y_individual, ...
                'AbsTol', testCase.tolerance, ...
                'Vectorized evaluation should match individual evaluations');
        end

        function testVectorizedDerivatives(testCase, kernelType)
            nodes = testCase.testNodes2D;
            basis = approx.basis.RbfKernelBasis(nodes, kernelType, 1.0);
            
            X = testCase.testPoints2D;
            dY_vectorized = basis.differentiate(X, [1, 0]);
            
            dY_individual = zeros(basis.nCodims, size(X, 2));
            for i = 1:size(X, 2)
                dY_individual(:, i) = basis.differentiate(X(:, i), [1, 0]);
            end
            
            testCase.verifyEqual(dY_vectorized, dY_individual, ...
                'AbsTol', testCase.tolerance, ...
                'Vectorized derivatives should match individual computations');
        end

        function testScalarVsVectorEpsilon(testCase, kernelType)
            nodes = testCase.testNodes1D;
            epsilon_scalar = 1.5;
            epsilon_vector = repmat(epsilon_scalar, 1, length(nodes));
            
            basis_scalar = approx.basis.RbfKernelBasis(nodes, kernelType, epsilon_scalar);
            basis_vector = approx.basis.RbfKernelBasis(nodes, kernelType, epsilon_vector);
            
            X = testCase.testPoints1D;
            Y_scalar = basis_scalar.evaluate(X);
            Y_vector = basis_vector.evaluate(X);
            
            testCase.verifyEqual(Y_scalar, Y_vector, ...
                'AbsTol', testCase.tolerance, ...
                'Scalar and vector epsilon should give same results when equal');
        end

        function testKernelSymmetry(testCase, kernelType)
            nodes = [0, 1];
            basis = approx.basis.RbfKernelBasis(nodes, kernelType, 1.0);
            
            x1 = 0.3;
            x2 = 0.7;
            
            Y1 = basis.evaluate(x1);
            Y2 = basis.evaluate(x2);
            
            d1_to_0 = abs(x1 - 0);
            d1_to_1 = abs(x1 - 1);
            d2_to_0 = abs(x2 - 0);
            d2_to_1 = abs(x2 - 1);
            
            testCase.verifyEqual(d1_to_1, d2_to_0, ...
                'AbsTol', testCase.tolerance);
            testCase.verifyEqual(Y1(2), Y2(1), ...
                'AbsTol', testCase.tolerance, ...
                'Kernel should be symmetric in distance');
        end

        function testLargeScaleEvaluation(testCase)
            nodes = rand(2, 50);
            basis = approx.basis.RbfKernelBasis(nodes, 'gaussian', 1.0);
            
            X = rand(2, 1000);
            Y = basis.evaluate(X);
            
            testCase.verifyEqual(size(Y), [50, 1000]);
            testCase.verifyTrue(all(isfinite(Y(:))));
            testCase.verifyTrue(all(Y(:) > 0));
        end

        function testHighDimensionalNodes(testCase)
            nodes = rand(5, 10);
            basis = approx.basis.RbfKernelBasis(nodes, 'multiquadric', 2.0);
            
            testCase.verifyEqual(basis.nDims, 5);
            testCase.verifyEqual(basis.nCodims, 10);
            
            X = rand(5, 20);
            Y = basis.evaluate(X);
            
            testCase.verifyEqual(size(Y), [10, 20]);
            testCase.verifyTrue(all(isfinite(Y(:))));
        end

        function testNumericalStabilityNearZero(testCase, kernelType)
            nodes = [0];
            basis = approx.basis.RbfKernelBasis(nodes, kernelType, 1.0);
            
            x_near_zero = [1e-15, 1e-12, 1e-10, 1e-8];
            
            for x = x_near_zero
                Y = basis.evaluate(x);
                dY = basis.differentiate(x, [1]);
                
                testCase.verifyTrue(isfinite(Y), ...
                    sprintf('Evaluation should be finite near zero for %s', kernelType));
                testCase.verifyTrue(isfinite(dY), ...
                    sprintf('Derivative should be finite near zero for %s', kernelType));
            end
        end

        function testInputValidationDimensions(testCase)
            nodes = testCase.testNodes2D;
            basis = approx.basis.RbfKernelBasis(nodes, 'gaussian', 1.0);
            
            testCase.verifyError(...
                @() basis.evaluate([0.5; 0.5; 0.5]), ...
                'core:function:Function:InvalidInput');
            
            testCase.verifyError(...
                @() basis.evaluate([]), ...
                'core:function:Function:InvalidInput');
        end

        function testSingleNodeBasis(testCase, kernelType)
            nodes = [0.5];
            basis = approx.basis.RbfKernelBasis(nodes, kernelType, 1.0);
            
            testCase.verifyEqual(basis.nCodims, 1);
            
            Y = basis.evaluate(0.5);
            testCase.verifyEqual(length(Y), 1);
            
            switch kernelType
                case 'gaussian'
                    expected = 1;
                case 'multiquadric'
                    expected = 1;
                case 'inverse_multiquadric'
                    expected = 1;
            end
            
            testCase.verifyEqual(Y, expected, 'AbsTol', testCase.tolerance);
        end

        function testDifferentEpsilonValues(testCase)
            nodes = [0, 1];
            
            epsilons = [0.1, 0.5, 1.0, 2.0, 10.0];
            
            for eps_val = epsilons
                basis = approx.basis.RbfKernelBasis(nodes, 'gaussian', eps_val);
                
                Y = basis.evaluate(0.5);
                testCase.verifyTrue(all(isfinite(Y)), ...
                    sprintf('Should be finite for epsilon = %g', eps_val));
                testCase.verifyTrue(all(Y > 0), ...
                    sprintf('Gaussian should be positive for epsilon = %g', eps_val));
            end
        end

        function testMemoryEfficiency(testCase)
            nodes = rand(3, 100);
            basis = approx.basis.RbfKernelBasis(nodes, 'inverse_multiquadric', 1.0);
            
            X = rand(3, 10000);
            
            try
                Y = basis.evaluate(X);
                testCase.verifyEqual(size(Y), [100, 10000]);
                
                dY = basis.differentiate(X, [1, 0, 0]);
                testCase.verifyEqual(size(dY), [100, 10000]);
                
            catch ME
                if contains(ME.identifier, 'OutOfMemory')
                    testCase.assumeFail('Insufficient memory for large-scale test');
                else
                    rethrow(ME);
                end
            end
        end
    end
end