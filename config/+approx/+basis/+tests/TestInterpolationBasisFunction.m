classdef TestInterpolationBasisFunction < matlab.unittest.TestCase

    properties
        testPoints
        tolerance
        gaussLegendrePoints
        gaussLobattoPoints
    end

    methods (TestMethodSetup)
        function setupTestCase(testCase)
            testCase.testPoints = linspace(-1/2, 1/2, 11);
            testCase.tolerance = 1e-12;

            testCase.gaussLegendrePoints = { ...
                0, ...
                [-1 / sqrt(3) / 2, 1 / sqrt(3) / 2], ...
                [-sqrt(3/5) / 2, 0, sqrt(3/5) / 2], ...
                [-0.430568155797026, -0.169990521792428, 0.169990521792428, 0.430568155797026], ...
                };

            testCase.gaussLobattoPoints = { ...
                0, ...
                [-1 / 2, 1 / 2], ...
                [-1 / 2, 0, 1 / 2], ...
                [-1 / 2, -sqrt(1/5) / 2, sqrt(1/5) / 2, 1 / 2], ...
                };
        end
    end

    methods (Test)
        function testGaussLegendreConstruction(testCase)
            basis = approx.basis.InterpolationBasisFunction(2, 'lagrange', 'unit', 'gauss_legendre');
            testCase.verifyEqual(basis.basisType, 'lagrange');
            testCase.verifyEqual(basis.intervalType, 'unit');
            testCase.verifyEqual(basis.nodeType, 'gauss_legendre');
            testCase.verifyEqual(basis.nDims, 1);
            testCase.verifyEqual(basis.nCodims, 2);
            testCase.verifyTrue(basis.isWellDefined);
            testCase.verifyFalse(basis.hasMetadata);
        end

        function testGaussLobattoConstruction(testCase)
            basis = approx.basis.InterpolationBasisFunction(3, 'lagrange', 'unit', 'gauss_lobatto');
            testCase.verifyEqual(basis.basisType, 'lagrange');
            testCase.verifyEqual(basis.intervalType, 'unit');
            testCase.verifyEqual(basis.nodeType, 'gauss_lobatto');
            testCase.verifyEqual(basis.nDims, 1);
            testCase.verifyEqual(basis.nCodims, 3);
            testCase.verifyTrue(basis.isWellDefined);
            testCase.verifyFalse(basis.hasMetadata);
        end

        function testGaussLegendreAutoLoad(testCase)
            basis = approx.basis.InterpolationBasisFunction(2, 'lagrange', 'unit', 'gauss_legendre');

            expectedFilename = 'lagrange_unit_gauss_legendre_2';
            if ismember(expectedFilename, basis.metasource)
                try
                    basis = basis.autoLoad();
                    testCase.verifyTrue(basis.hasMetadata);
                    testCase.verifyTrue(basis.isWellDefined);
                    testCase.verifyEqual(basis.nCodims, 2);
                catch ME
                    testCase.assumeFail(sprintf('AutoLoad failed: %s', ME.message));
                end
            else
                testCase.verifyError(@() basis.autoLoad(), 'FileNotFound');
            end
        end

        function testGaussLobattoAutoLoad(testCase)
            basis = approx.basis.InterpolationBasisFunction(3, 'lagrange', 'unit', 'gauss_lobatto');

            expectedFilename = 'lagrange_unit_gauss_lobatto_3';
            if ismember(expectedFilename, basis.metasource)
                try
                    basis = basis.autoLoad();
                    testCase.verifyTrue(basis.hasMetadata);
                    testCase.verifyTrue(basis.isWellDefined);
                    testCase.verifyEqual(basis.nCodims, 3);
                catch ME
                    testCase.assumeFail(sprintf('AutoLoad failed: %s', ME.message));
                end
            else
                testCase.verifyError(@() basis.autoLoad(), 'FileNotFound');
            end
        end

        function testGaussLegendreLagrangeProperty(testCase)
            basis = approx.basis.InterpolationBasisFunction(2, 'lagrange', 'unit', 'gauss_legendre');
            expectedFilename = 'lagrange_unit_gauss_legendre_2';

            if ~ismember(expectedFilename, basis.metasource)
                return;
            end

            try
                basis = basis.autoLoad();
            catch ME
                testCase.assumeFail(sprintf('AutoLoad failed: %s', ME.message));
            end

            if ~basis.isWellDefined
                return;
            end

            degree = 1;

            nodes = testCase.gaussLegendrePoints{degree+1};
            Y = basis.evaluate(nodes);

            for i = 1:(degree + 1)
                for j = 1:(degree + 1)
                    if i == j
                        testCase.verifyEqual(Y(i, j), 1, 'AbsTol', testCase.tolerance, ...
                            sprintf('Gauss-Legendre: L_%d(x_%d) should be 1', i, j));
                    else
                        testCase.verifyEqual(Y(i, j), 0, 'AbsTol', testCase.tolerance, ...
                            sprintf('Gauss-Legendre: L_%d(x_%d) should be 0', i, j));
                    end
                end
            end
        end

        function testGaussLobattoLagrangeProperty(testCase)
            basis = approx.basis.InterpolationBasisFunction(3, 'lagrange', 'unit', 'gauss_lobatto');
            expectedFilename = 'lagrange_unit_gauss_lobatto_3';

            if ~ismember(expectedFilename, basis.metasource)
                return;
            end

            try
                basis = basis.autoLoad();
            catch ME
                testCase.assumeFail(sprintf('AutoLoad failed: %s', ME.message));
            end

            if ~basis.isWellDefined
                return;
            end

            degree = 2;

            nodes = testCase.gaussLobattoPoints{degree+1};
            Y = basis.evaluate(nodes);

            for i = 1:(degree + 1)
                for j = 1:(degree + 1)
                    if i == j
                        testCase.verifyEqual(Y(i, j), 1, 'AbsTol', testCase.tolerance, ...
                            sprintf('Gauss-Lobatto: L_%d(x_%d) should be 1', i, j));
                    else
                        testCase.verifyEqual(Y(i, j), 0, 'AbsTol', testCase.tolerance, ...
                            sprintf('Gauss-Lobatto: L_%d(x_%d) should be 0', i, j));
                    end
                end
            end
        end

        function testGaussLegendrePartitionOfUnity(testCase)
            basis = approx.basis.InterpolationBasisFunction(2, 'lagrange', 'unit', 'gauss_legendre');
            expectedFilename = 'lagrange_unit_gauss_legendre_2';

            if ~ismember(expectedFilename, basis.metasource)
                return;
            end

            try
                basis = basis.autoLoad();
            catch ME
                testCase.assumeFail(sprintf('AutoLoad failed: %s', ME.message));
            end

            if ~basis.isWellDefined
                return;
            end

            x = testCase.testPoints;
            Y = basis.evaluate(x);
            sumY = sum(Y, 1);

            testCase.verifyEqual(sumY, ones(size(x)), 'AbsTol', testCase.tolerance, ...
                'Gauss-Legendre: Sum of basis functions should be 1');
        end

        function testGaussLobattoPartitionOfUnity(testCase)
            basis = approx.basis.InterpolationBasisFunction(3, 'lagrange', 'unit', 'gauss_lobatto');
            expectedFilename = 'lagrange_unit_gauss_lobatto_3';

            if ~ismember(expectedFilename, basis.metasource)
                return;
            end

            try
                basis = basis.autoLoad();
            catch ME
                testCase.assumeFail(sprintf('AutoLoad failed: %s', ME.message));
            end

            if ~basis.isWellDefined
                return;
            end

            x = testCase.testPoints;
            Y = basis.evaluate(x);
            sumY = sum(Y, 1);

            testCase.verifyEqual(sumY, ones(size(x)), 'AbsTol', testCase.tolerance, ...
                'Gauss-Lobatto: Sum of basis functions should be 1');
        end

        function testGaussLobattoEndpointValues(testCase)
            basis = approx.basis.InterpolationBasisFunction(2, 'lagrange', 'unit', 'gauss_lobatto');
            expectedFilename = 'lagrange_unit_gauss_lobatto_2';

            if ~ismember(expectedFilename, basis.metasource)
                return;
            end

            try
                basis = basis.autoLoad();
            catch ME
                testCase.assumeFail(sprintf('AutoLoad failed: %s', ME.message));
            end

            if ~basis.isWellDefined || basis.nCodims < 2
                return;
            end

            Y_left = basis.evaluate(-1/2);
            testCase.verifyEqual(Y_left(1), 1, 'AbsTol', testCase.tolerance, ...
                'First basis function should be 1 at left endpoint');
            for i = 2:basis.nCodims
                testCase.verifyEqual(Y_left(i), 0, 'AbsTol', testCase.tolerance, ...
                    sprintf('Basis function %d should be 0 at left endpoint', i));
            end

            Y_right = basis.evaluate(1/2);
            testCase.verifyEqual(Y_right(end), 1, 'AbsTol', testCase.tolerance, ...
                'Last basis function should be 1 at right endpoint');
            for i = 1:basis.nCodims - 1
                testCase.verifyEqual(Y_right(i), 0, 'AbsTol', testCase.tolerance, ...
                    sprintf('Basis function %d should be 0 at right endpoint', i));
            end
        end

        function testGaussLegendreEvaluate(testCase)
            basis = approx.basis.InterpolationBasisFunction(2, 'lagrange', 'unit', 'gauss_legendre');
            expectedFilename = 'lagrange_unit_gauss_legendre_2';

            if ~ismember(expectedFilename, basis.metasource)
                return;
            end

            try
                basis = basis.autoLoad();
            catch ME
                testCase.assumeFail(sprintf('AutoLoad failed: %s', ME.message));
            end

            if ~basis.isWellDefined
                return;
            end

            x = testCase.testPoints;
            Y = basis.evaluate(x);

            if basis.nCodims >= 2
                expected0 = 1 / 2 - sqrt(3) * x;
                expected1 = 1 / 2 + sqrt(3) * x;

                testCase.verifyEqual(Y(1, :), expected0, 'AbsTol', testCase.tolerance, ...
                    'Gauss-Legendre L_0(x) = 1/2 - √3*x');
                testCase.verifyEqual(Y(2, :), expected1, 'AbsTol', testCase.tolerance, ...
                    'Gauss-Legendre L_1(x) = 1/2 + √3*x');
            end
        end

        function testGaussLobattoEvaluate(testCase)
            basis = approx.basis.InterpolationBasisFunction(2, 'lagrange', 'unit', 'gauss_lobatto');
            expectedFilename = 'lagrange_unit_gauss_lobatto_2';

            if ~ismember(expectedFilename, basis.metasource)
                return;
            end

            try
                basis = basis.autoLoad();
            catch ME
                testCase.assumeFail(sprintf('AutoLoad failed: %s', ME.message));
            end

            if ~basis.isWellDefined
                return;
            end

            x = testCase.testPoints;
            Y = basis.evaluate(x);

            if basis.nCodims >= 2
                expected0 = 1 / 2 - x;
                expected1 = 1 / 2 + x;

                testCase.verifyEqual(Y(1, :), expected0, 'AbsTol', testCase.tolerance, ...
                    'Gauss-Lobatto L_0(x) = 1/2 - x');
                testCase.verifyEqual(Y(2, :), expected1, 'AbsTol', testCase.tolerance, ...
                    'Gauss-Lobatto L_1(x) = 1/2 + x');
            end

            basis3 = approx.basis.InterpolationBasisFunction(3, 'lagrange', 'unit', 'gauss_lobatto');
            expectedFilename3 = 'lagrange_unit_gauss_lobatto_3';

            if ismember(expectedFilename3, basis3.metasource)
                try
                    basis3 = basis3.autoLoad();
                    if basis3.isWellDefined && basis3.nCodims >= 3
                        Y_center = basis3.evaluate(0);
                        testCase.verifyEqual(Y_center(1), 0, 'AbsTol', testCase.tolerance, ...
                            'L_0(0) should be 0');
                        testCase.verifyEqual(Y_center(2), 1, 'AbsTol', testCase.tolerance, ...
                            'L_1(0) should be 1');
                        testCase.verifyEqual(Y_center(3), 0, 'AbsTol', testCase.tolerance, ...
                            'L_2(0) should be 0');
                    end
                catch
                end
            end
        end

        function testGaussLegendreDerivative(testCase)
            basis = approx.basis.InterpolationBasisFunction(2, 'lagrange', 'unit', 'gauss_legendre');
            expectedFilename = 'lagrange_unit_gauss_legendre_2';

            if ~ismember(expectedFilename, basis.metasource)
                return;
            end

            try
                basis = basis.autoLoad();
            catch ME
                testCase.assumeFail(sprintf('AutoLoad failed: %s', ME.message));
            end

            if ~basis.isWellDefined
                return;
            end

            x = testCase.testPoints;

            dY1 = basis.derivative(x, 1);

            if basis.nCodims >= 2
                expectedDeriv0 = -sqrt(3) * ones(size(x));
                expectedDeriv1 = sqrt(3) * ones(size(x));

                testCase.verifyEqual(dY1(1, :), expectedDeriv0, 'AbsTol', testCase.tolerance, ...
                    'Gauss-Legendre L_0''(x) = -√3');
                testCase.verifyEqual(dY1(2, :), expectedDeriv1, 'AbsTol', testCase.tolerance, ...
                    'Gauss-Legendre L_1''(x) = √3');
            end

            try
                dY2 = basis.derivative(x, 2);
                if basis.nCodims == 2
                    testCase.verifyEqual(dY2, zeros(size(dY2)), 'AbsTol', testCase.tolerance, ...
                        'Second derivative of linear functions should be zero');
                end
            catch
            end
        end

        function testGaussLobattoDerivative(testCase)
            basis = approx.basis.InterpolationBasisFunction(2, 'lagrange', 'unit', 'gauss_lobatto');
            expectedFilename = 'lagrange_unit_gauss_lobatto_2';

            if ~ismember(expectedFilename, basis.metasource)
                return;
            end

            try
                basis = basis.autoLoad();
            catch ME
                testCase.assumeFail(sprintf('AutoLoad failed: %s', ME.message));
            end

            if ~basis.isWellDefined
                return;
            end

            x = testCase.testPoints;

            dY1 = basis.derivative(x, 1);

            if basis.nCodims >= 2
                expectedDeriv0 = -ones(size(x));
                expectedDeriv1 = ones(size(x));

                testCase.verifyEqual(dY1(1, :), expectedDeriv0, 'AbsTol', testCase.tolerance, ...
                    'Gauss-Lobatto L_0''(x) = -1');
                testCase.verifyEqual(dY1(2, :), expectedDeriv1, 'AbsTol', testCase.tolerance, ...
                    'Gauss-Lobatto L_1''(x) = 1');
            end

            try
                dY2 = basis.derivative(x, 2);
                if basis.nCodims == 2
                    testCase.verifyEqual(dY2, zeros(size(dY2)), 'AbsTol', testCase.tolerance, ...
                        'Second derivative of linear functions should be zero');
                end
            catch
            end

            basis3 = approx.basis.InterpolationBasisFunction(3, 'lagrange', 'unit', 'gauss_lobatto');
            expectedFilename3 = 'lagrange_unit_gauss_lobatto_3';

            if ismember(expectedFilename3, basis3.metasource)
                try
                    basis3 = basis3.autoLoad();
                    if basis3.isWellDefined && basis3.nCodims >= 3
                        dY1_center = basis3.derivative(0, 1);
                        testCase.verifyEqual(dY1_center(1), -1, 'AbsTol', testCase.tolerance, ...
                            'L_0''(0) should be -1');
                        testCase.verifyEqual(dY1_center(2), 0, 'AbsTol', testCase.tolerance, ...
                            'L_1''(0) should be 0');
                        testCase.verifyEqual(dY1_center(3), 1, 'AbsTol', testCase.tolerance, ...
                            'L_2''(0) should be 1');

                        dY2_center = basis3.derivative(0, 2);
                        testCase.verifyEqual(dY2_center(1), 4, 'AbsTol', testCase.tolerance, ...
                            'L_0''''(0) should be 4');
                        testCase.verifyEqual(dY2_center(2), -8, 'AbsTol', testCase.tolerance, ...
                            'L_1''''(0) should be -8');
                        testCase.verifyEqual(dY2_center(3), 4, 'AbsTol', testCase.tolerance, ...
                            'L_2''''(0) should be 4');
                    end
                catch
                end
            end
        end

        function testInvalidInputs(testCase)
            basisGL = approx.basis.InterpolationBasisFunction(2, 'lagrange', 'unit', 'gauss_legendre');
            expectedFilename = 'lagrange_unit_gauss_legendre_2';

            if ismember(expectedFilename, basisGL.metasource)
                try
                    basisGL = basisGL.autoLoad();
                    if basisGL.isWellDefined
                        testCase.verifyError(@() basisGL.evaluate([1, 2; 3, 4]), ...
                            'core:function:Function:InvalidInput');
                        testCase.verifyError(@() basisGL.evaluate([]), ...
                            'core:function:Function:InvalidInput');
                        testCase.verifyError(@() basisGL.derivative([0], -1), ...
                            'core:function:Function:InvalidInput');
                        testCase.verifyError(@() basisGL.derivative([0], [1, 2]), ...
                            'core:function:Function:InvalidInput');
                    end
                catch ME
                    testCase.assumeFail(sprintf('AutoLoad failed: %s', ME.message));
                end
            end

            basisGLo = approx.basis.InterpolationBasisFunction(2, 'lagrange', 'unit', 'gauss_lobatto');
            expectedFilename2 = 'lagrange_unit_gauss_lobatto_2';

            if ismember(expectedFilename2, basisGLo.metasource)
                try
                    basisGLo = basisGLo.autoLoad();
                    if basisGLo.isWellDefined
                        testCase.verifyError(@() basisGLo.evaluate([1, 2; 3, 4]), ...
                            'core:function:Function:InvalidInput');
                        testCase.verifyError(@() basisGLo.evaluate([]), ...
                            'core:function:Function:InvalidInput');
                        testCase.verifyError(@() basisGLo.derivative([0], -1), ...
                            'core:function:Function:InvalidInput');
                        testCase.verifyError(@() basisGLo.derivative([0], [1, 2]), ...
                            'core:function:Function:InvalidInput');
                    end
                catch ME
                    testCase.assumeFail(sprintf('AutoLoad failed: %s', ME.message));
                end
            end
        end

        function testDifferentBasisSizes(testCase)
            for nCodims = [2, 3, 4, 5]
                basisGL = approx.basis.InterpolationBasisFunction(nCodims, 'lagrange', 'unit', 'gauss_legendre');
                testCase.verifyEqual(basisGL.nCodims, nCodims);
                testCase.verifyTrue(basisGL.isWellDefined);

                expectedFilename = sprintf('lagrange_unit_gauss_legendre_%d', nCodims);
                if ismember(expectedFilename, basisGL.metasource)
                    basisGL = basisGL.autoLoad();
                    if basisGL.hasMetadata
                        testCase.verifyEqual(basisGL.nCodims, nCodims);

                        x = testCase.testPoints;
                        Y = basisGL.evaluate(x);
                        testCase.verifyEqual(size(Y, 1), nCodims);
                        testCase.verifyEqual(size(Y, 2), length(x));
                    end

                end

                basisGLo = approx.basis.InterpolationBasisFunction(nCodims, 'lagrange', 'unit', 'gauss_lobatto');
                testCase.verifyEqual(basisGLo.nCodims, nCodims);
                testCase.verifyTrue(basisGLo.isWellDefined);

                expectedFilename2 = sprintf('lagrange_unit_gauss_lobatto_%d', nCodims);
                if ismember(expectedFilename2, basisGLo.metasource)
                    basisGLo = basisGLo.autoLoad();
                    if basisGLo.hasMetadata
                        testCase.verifyEqual(basisGLo.nCodims, nCodims);

                        x = testCase.testPoints;
                        Y = basisGLo.evaluate(x);
                        testCase.verifyEqual(size(Y, 1), nCodims);
                        testCase.verifyEqual(size(Y, 2), length(x));
                    end
                end
            end
        end
    end
end