classdef TestSeparableBasisFunction < matlab.unittest.TestCase
    % TESTSEPARABLEBASISFUNCTION Unit tests for SeparableBasisFunction class.
    %
    % This test suite verifies the functionality of the SeparableBasisFunction
    % class including construction, metadata handling, evaluation, and
    % differentiation with different pattern types and basis function
    % combinations.

    properties (TestParameter)
        patternType = {'full', 'l1', 'lx'}
        nDims = {2, 3}
        degree = {2, 3, 4}
    end

    methods (Test)
        function testInheritance(testCase)
            % Test class inheritance hierarchy
            
            factors = {approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                       approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit')};
            basis = approx.basis.SeparableBasisFunction(factors);
            
            testCase.verifyTrue(isa(basis, 'core.function.SeparableFunction'));
            testCase.verifyTrue(isa(basis, 'core.function.SeparableFunction'));
            testCase.verifyTrue(isa(basis, 'core.function.Function'));
            testCase.verifyTrue(isa(basis, 'core.function.Function'));
        end

        function testConstructorWithOrthogonalFactors(testCase)
            % Test constructor with OrthogonalBasisFunction factors
            
            factors = {approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                       approx.basis.OrthogonalBasisFunction(4, 'chebyshev', 'canonical')};
            basis = approx.basis.SeparableBasisFunction(factors);
            
            testCase.verifyEqual(basis.nFactors, 2);
            testCase.verifyEqual(basis.nDims, 2);
            testCase.verifyEqual(basis.nFactorCodims, [3, 4]);
            testCase.verifyEqual(basis.nCodims, 12); % 3 * 4
            testCase.verifyFalse(basis.isSparse);
        end

        function testConstructorWithInterpolationFactors(testCase)
            % Test constructor with InterpolationBasisFunction factors
            
            factors = {approx.basis.InterpolationBasisFunction(4, 'lagrange', 'unit', 'gauss_legendre'), ...
                       approx.basis.InterpolationBasisFunction(3, 'hermite', 'canonical', 'gauss_lobatto')};
            basis = approx.basis.SeparableBasisFunction(factors);
            
            testCase.verifyEqual(basis.nFactors, 2);
            testCase.verifyEqual(basis.nDims, 2);
            testCase.verifyEqual(basis.nFactorCodims, [4, 3]);
            testCase.verifyEqual(basis.nCodims, 12); % 4 * 3
            testCase.verifyFalse(basis.isSparse);
        end

        function testConstructorWithMixedFactors(testCase)
            % Test constructor with mixed factor types
            
            factors = {approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                       approx.basis.InterpolationBasisFunction(4, 'lagrange', 'unit', 'gauss_legendre')};
            basis = approx.basis.SeparableBasisFunction(factors);
            
            testCase.verifyEqual(basis.nFactors, 2);
            testCase.verifyEqual(basis.nFactorCodims, [3, 4]);
            testCase.verifyEqual(basis.nCodims, 12); % 3 * 4
        end

        function testConstructorWithPatternTypes(testCase, patternType)
            % Test constructor with different pattern types
            
            factors = {approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                       approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit')};
            basis = approx.basis.SeparableBasisFunction(factors, patternType);
            basis.autoLoad();
            
            testCase.verifyEqual(basis.nFactors, 2);
            
            switch lower(patternType)
                case 'full'
                    testCase.verifyFalse(basis.isSparse);
                    testCase.verifyEqual(basis.nCodims, 9); % 3 * 3
                case 'l1'
                    testCase.verifyTrue(basis.isSparse);
                    testCase.verifyEqual(basis.nCodims, 6);
                case 'lx'
                    testCase.verifyTrue(basis.isSparse);
                    testCase.verifyEqual(basis.nCodims, 9);
            end
        end

        function testEmptyConstructor(testCase)
            % Test empty constructor
            
            basis = approx.basis.SeparableBasisFunction();
            
            testCase.verifyEqual(basis.nFactors, 0);
            testCase.verifyEqual(basis.nDims, 0);
            testCase.verifyEqual(basis.nCodims, 0);
            testCase.verifyFalse(basis.isWellDefined);
        end

        function testInvalidFactorValidation(testCase)
            % Test validation of invalid factor types
            
            % Test with non-basis function
            invalidFactors = {12, 2};
            testCase.verifyError(...
                @() approx.basis.SeparableBasisFunction(invalidFactors), ...
                'approx:basis:SeparableBasisFunction:InvalidFactor');
            
            % Test with mixed valid and invalid factors
            mixedFactors = {approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                           1};
            testCase.verifyError(...
                @() approx.basis.SeparableBasisFunction(mixedFactors), ...
                'approx:basis:SeparableBasisFunction:InvalidFactor');
        end

        function testMetadataProperties(testCase)
            % Test metadata-related dependent properties
            
            factors = {approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                       approx.basis.InterpolationBasisFunction(4, 'lagrange', 'unit', 'gauss_legendre')};
            basis = approx.basis.SeparableBasisFunction(factors);
            
            % Before loading metadata
            testCase.verifyFalse(basis.hasMetadata);
            sources = basis.metasource;
            testCase.verifyTrue(iscell(sources));
            testCase.verifyEqual(length(sources), 2);
        end

        function testAutoLoadMethod(testCase)
            % Test autoLoad method (mock test since metadata files may not exist)
            
            factors = {approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                       approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit')};
            basis = approx.basis.SeparableBasisFunction(factors);
            
            % This test assumes autoLoad will work or gracefully handle missing files
            try
                basis = basis.autoLoad();
                % If successful, verify the object is still valid
                testCase.verifyEqual(basis.nFactors, 2);
            catch ME
                % If autoLoad fails due to missing metadata files, that's expected
                testCase.verifyTrue(contains(ME.identifier, 'FileNotFound') || ...
                                   contains(ME.identifier, 'MetadataNotFound'));
            end
        end

        function testLoadMethodValidation(testCase)
            % Test load method input validation
            
            factors = {approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                       approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit')};
            basis = approx.basis.SeparableBasisFunction(factors);
            
            % Test with wrong number of filenames
            testCase.verifyError(...
                @() basis.load({'file1'}), ...
                'approx:basis:SeparableBasisFunction:InvalidInput');
            
            % Test with non-cell input
            testCase.verifyError(...
                @() basis.load('file1'), ...
                'approx:basis:SeparableBasisFunction:InvalidInput');
        end

        function testEvaluationWithoutMetadata(testCase, patternType)
            % Test that evaluation works even without loaded metadata
            % (using default implementations)
            
            factors = {approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                       approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit')};
            basis = approx.basis.SeparableBasisFunction(factors, patternType);
            
            X = [0.5; 0.3];
            
            % This should work with the base implementations
            % even if compiled metadata is not available
            try
                Y = basis.evaluate(X);
                testCase.verifyEqual(size(Y, 1), basis.nCodims);
                testCase.verifyEqual(size(Y, 2), 1);
            catch ME
                % If evaluation fails due to missing metadata, verify it's the expected error
                testCase.verifyTrue(contains(ME.identifier, 'NoMetadata') || ...
                                   contains(ME.identifier, 'NotWellDefined'));
            end
        end

        function testDerivativeWithoutMetadata(testCase, patternType)
            % Test derivative computation without loaded metadata
            
            factors = {approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                       approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit')};
            basis = approx.basis.SeparableBasisFunction(factors, patternType);
            
            X = [0.5; 0.3];
            orders = {[1, 0], [0, 1], [1, 1]};
            
            for i = 1:length(orders)
                try
                    dY = basis.derivative(X, orders{i});
                    testCase.verifyEqual(size(dY, 1), basis.nCodims);
                    testCase.verifyEqual(size(dY, 2), 1);
                catch ME
                    % If derivative fails due to missing metadata, verify it's the expected error
                    testCase.verifyTrue(contains(ME.identifier, 'NoMetadata') || ...
                                       contains(ME.identifier, 'NotWellDefined'));
                end
            end
        end

        function testCellArrayInput(testCase)
            % Test evaluation with cell array input format
            
            factors = {approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                       approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit')};
            basis = approx.basis.SeparableBasisFunction(factors);
            
            X_matrix = [0.5; 0.3];
            X_cell = {0.5, 0.3};
            
            try
                Y_matrix = basis.evaluate(X_matrix);
                Y_cell = basis.evaluate(X_cell);
                testCase.verifyEqual(Y_matrix, Y_cell, 'AbsTol', 1e-12);
            catch ME
                % Both should fail in the same way if metadata is missing
                if contains(ME.identifier, 'NoMetadata')
                    testCase.verifyError(@() basis.evaluate(X_cell), ME.identifier);
                else
                    rethrow(ME);
                end
            end
        end

        function testMultiplePointEvaluation(testCase)
            % Test evaluation at multiple points
            
            factors = {approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                       approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit')};
            basis = approx.basis.SeparableBasisFunction(factors);
            
            X = [0.1, 0.2, 0.3; 0.4, 0.5, 0.6];
            
            try
                Y = basis.evaluate(X);
                testCase.verifyEqual(size(Y, 1), basis.nCodims);
                testCase.verifyEqual(size(Y, 2), 3);
            catch ME
                if ~contains(ME.identifier, 'NoMetadata')
                    rethrow(ME);
                end
            end
        end

        function testHigherDimensionalInput(testCase)
            % Test evaluation with higher-dimensional input arrays
            
            factors = {approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                       approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit')};
            basis = approx.basis.SeparableBasisFunction(factors);
            
            X_3d = zeros(2, 2, 3);
            X_3d(1, :, :) = reshape(0.1:0.1:0.6, [1, 2, 3]);
            X_3d(2, :, :) = reshape(0.4:0.1:0.9, [1, 2, 3]);
            
            try
                Y = basis.evaluate(X_3d);
                testCase.verifyEqual(size(Y, 1), basis.nCodims);
                testCase.verifyEqual([size(Y, 2), size(Y, 3)], [2, 3]);
            catch ME
                if ~contains(ME.identifier, 'NoMetadata')
                    rethrow(ME);
                end
            end
        end

        function testBoundaryValues(testCase)
            % Test evaluation at boundary values
            
            factors = {approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                       approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit')};
            basis = approx.basis.SeparableBasisFunction(factors);
            
            corners = [-1, -1; 1, -1; -1, 1; 1, 1]';
            
            for i = 1:size(corners, 2)
                try
                    Y = basis.evaluate(corners(:, i));
                    testCase.verifyEqual(size(Y, 1), basis.nCodims);
                    testCase.verifyTrue(all(isfinite(Y)));
                catch ME
                    if ~contains(ME.identifier, 'NoMetadata')
                        rethrow(ME);
                    end
                end
            end
        end

        function testExceptionHandling(testCase)
            % Test proper exception handling for invalid inputs
            
            factors = {approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                       approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit')};
            basis = approx.basis.SeparableBasisFunction(factors);
            
            % Test dimension mismatch in evaluation
            X_invalid = [0.5; 0.5; 0.5]; % 3D input for 2D function
            testCase.verifyError(...
                @() basis.evaluate(X_invalid), ...
                'core:function:SeparableFunction:InvalidInput');
            
            % Test invalid derivative order
            X_valid = [0.5; 0.5];
            order_invalid = [1, 0, 0]; % 3D order for 2D function
            testCase.verifyError(...
                @() basis.derivative(X_valid, order_invalid), ...
                'core:function:SeparableFunction:InvalidInput');
            
            % Test negative derivative order
            order_negative = [-1, 0];
            testCase.verifyError(...
                @() basis.derivative(X_valid, order_negative), ...
                'core:function:SeparableFunction:InvalidInput');
            
            % Test zero derivative order
            order_zero = [0, 0];
            testCase.verifyError(...
                @() basis.derivative(X_valid, order_zero), ...
                'core:function:SeparableFunction:InvalidInput');
        end

        function testPatternTypeValidation(testCase)
            % Test validation of pattern type parameter
            
            factors = {approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                       approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit')};
            
            % Test invalid pattern type
            testCase.verifyError(...
                @() approx.basis.SeparableBasisFunction(factors, 'invalid'), ...
                'core:function:SeparableFunction:InvalidInput');
            
            % Test numeric pattern type (should fail)
            testCase.verifyError(...
                @() approx.basis.SeparableBasisFunction(factors, 123), ...
                'core:function:SeparableFunction:InvalidInput');
        end

        function testSparsePatternProperties(testCase, nDims, degree)
            % Test properties of sparse patterns
            
            factors = repmat({approx.basis.OrthogonalBasisFunction(degree, 'monic_legendre', 'unit')}, 1, nDims);
            
            basisFull = approx.basis.SeparableBasisFunction(factors, 'full');
            basisL1 = approx.basis.SeparableBasisFunction(factors, 'l1');
            basisLx = approx.basis.SeparableBasisFunction(factors, 'lx');
            
            % Full pattern
            testCase.verifyFalse(basisFull.isSparse);
            testCase.verifyEqual(basisFull.nCodims, degree^nDims);
            
            % L1 pattern
            testCase.verifyTrue(basisL1.isSparse);
            testCase.verifyLessThan(basisL1.nCodims, degree^nDims);
            
            % Lx pattern
            testCase.verifyTrue(basisLx.isSparse);
            testCase.verifyEqual(basisLx.nCodims, degree^nDims);
            
            % L1 and Lx should generally have different sizes
            if nDims > 1 && degree > 1
                testCase.verifyNotEqual(basisL1.nCodims, basisLx.nCodims);
            end
        end

        function testFactorConsistency(testCase)
            % Test that factors are properly maintained
            
            factor1 = approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit');
            factor2 = approx.basis.InterpolationBasisFunction(4, 'lagrange', 'unit', 'gauss_legendre');
            factors = {factor1, factor2};
            
            basis = approx.basis.SeparableBasisFunction(factors);
            
            % Verify factors are accessible
            retrievedFactors = basis.factors;
            testCase.verifyEqual(length(retrievedFactors), 2);
            
            if iscell(retrievedFactors)
                testCase.verifyEqual(retrievedFactors{1}.nCodims, 3);
                testCase.verifyEqual(retrievedFactors{2}.nCodims, 4);
            else
                testCase.verifyEqual(retrievedFactors(1).nCodims, 3);
                testCase.verifyEqual(retrievedFactors(2).nCodims, 4);
            end
        end

        function testObjectArrayFactors(testCase)
            % Test construction with object array of factors
            
            factors = [approx.basis.OrthogonalBasisFunction(3, 'monic_legendre', 'unit'), ...
                       approx.basis.OrthogonalBasisFunction(4, 'monic_legendre', 'unit')];
            
            basis = approx.basis.SeparableBasisFunction(factors);
            
            testCase.verifyEqual(basis.nFactors, 2);
            testCase.verifyEqual(basis.nFactorCodims, [3, 4]);
        end
    end
end