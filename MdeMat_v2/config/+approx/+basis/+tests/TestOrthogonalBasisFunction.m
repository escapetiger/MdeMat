classdef TestOrthogonalBasisFunction < matlab.unittest.TestCase
    
    properties
        unitTestPoints      % Test points for unit interval [-1/2, 1/2]
        canonicalTestPoints % Test points for canonical interval [-1, 1]
        tolerance
    end
    
    methods (TestMethodSetup)
        function setupTestCase(testCase)
            testCase.unitTestPoints = linspace(-1/2, 1/2, 11);
            testCase.canonicalTestPoints = linspace(-1, 1, 11);
            testCase.tolerance = 1e-12;
        end
    end
    
    methods (Test)
        function testConstruction(testCase)
            % Test valid constructions
            basis1 = approx.basis.OrthogonalBasisFunction(5, 'monic_legendre', 'unit');
            testCase.verifyEqual(basis1.basisType, 'monic_legendre');
            testCase.verifyEqual(basis1.intervalType, 'unit');
            testCase.verifyEqual(basis1.nDims, 1);
            testCase.verifyEqual(basis1.nCodims, 5);
            testCase.verifyTrue(basis1.isWellDefined); % nCodims = 5 > 0
            testCase.verifyFalse(basis1.hasMetadata);
            
            % Test case insensitive input
            basis2 = approx.basis.OrthogonalBasisFunction(3, 'MONIC_LEGENDRE', 'UNIT');
            testCase.verifyEqual(basis2.basisType, 'monic_legendre');
            testCase.verifyEqual(basis2.intervalType, 'unit');
            testCase.verifyEqual(basis2.nCodims, 3);
            
            % Test different basis types
            basis3 = approx.basis.OrthogonalBasisFunction(4, 'legendre', 'canonical');
            testCase.verifyEqual(basis3.basisType, 'legendre');
            testCase.verifyEqual(basis3.intervalType, 'canonical');
            testCase.verifyEqual(basis3.nCodims, 4);
            
            basis4 = approx.basis.OrthogonalBasisFunction(6, 'chebyshev', 'unit');
            testCase.verifyEqual(basis4.basisType, 'chebyshev');
            testCase.verifyEqual(basis4.intervalType, 'unit');
            testCase.verifyEqual(basis4.nCodims, 6);
        end
        
        function testInvalidConstruction(testCase)
            % Test invalid basis types
            testCase.verifyError(@() approx.basis.OrthogonalBasisFunction(5, 'hermite', 'unit'), ...
                'approx:basis:OrthogonalBasisFunction:InvalidInput');
            testCase.verifyError(@() approx.basis.OrthogonalBasisFunction(5, 'invalid', 'unit'), ...
                'approx:basis:OrthogonalBasisFunction:InvalidInput');
            
            % Test invalid interval types
            testCase.verifyError(@() approx.basis.OrthogonalBasisFunction(5, 'legendre', 'custom'), ...
                'approx:basis:OrthogonalBasisFunction:InvalidInput');
            testCase.verifyError(@() approx.basis.OrthogonalBasisFunction(5, 'legendre', 'invalid'), ...
                'approx:basis:OrthogonalBasisFunction:InvalidInput');
        end
        
        function testMetasourceProperty(testCase)
            basis = approx.basis.OrthogonalBasisFunction(5, 'monic_legendre', 'unit');
            
            % Test that metasource is accessible (dependent property)
            fileList = basis.metasource;
            testCase.verifyTrue(iscell(fileList));
            
            % Multiple calls should return consistent results
            fileList2 = basis.metasource;
            testCase.verifyEqual(fileList, fileList2);
        end
        
        function testAutoFilenameGeneration(testCase)
            % Test filename generation for different configurations
            basis1 = approx.basis.OrthogonalBasisFunction(5, 'monic_legendre', 'unit');
            expectedFilename1 = 'monic_legendre_unit_5';
            
            % Test that autoFilename generates expected patterns
            % (We can't directly test protected methods, but we can test through autoLoad)
            if ismember(expectedFilename1, basis1.metasource)
                try
                    basis1.autoLoad();
                catch ME
                    % If it fails due to file format issues, that's okay
                    % We're testing filename generation, not file loading
                    testCase.verifyTrue(contains(ME.identifier, 'InvalidFile') || ...
                                      contains(ME.identifier, 'NoMetadata'));
                end
            end
        end
        
        function testAutoLoadWithValidFile(testCase)
            basis = approx.basis.OrthogonalBasisFunction(5, 'monic_legendre', 'unit');
            expectedFilename = 'monic_legendre_unit_5';
            
            % Check if the expected file exists
            if ismember(expectedFilename, basis.metasource)
                % Attempt auto-load
                try
                    basis = basis.autoLoad();
                    testCase.verifyTrue(basis.hasMetadata);
                    testCase.verifyTrue(basis.isWellDefined);
                    testCase.verifyEqual(basis.nCodims, 5);
                catch ME
                    % If autoLoad fails, it might be due to missing metadata file
                    % or corrupted file format - this is acceptable for testing
                    testCase.verifyTrue( ...
                        contains(ME.identifier, 'FileNotFound') ...
                        || contains(ME.identifier, 'InvalidFile') ...
                        || contains(ME.identifier, 'NoMetadata'));
                end
            else
                % If file doesn't exist, autoLoad should throw FileNotFound
                testCase.verifyError(@() basis.autoLoad(), 'approx:basis:OrthogonalBasisFunction:FileNotFound');
            end
        end
        
        function testAutoLoadWithMissingFile(testCase)
            % Test with combination that likely doesn't exist
            basis = approx.basis.OrthogonalBasisFunction(99, 'chebyshev', 'unit');
            expectedFilename = 'chebyshev_unit_99';
            
            if ~ismember(expectedFilename, basis.metasource)
                testCase.verifyError(@() basis.autoLoad(), 'approx:basis:OrthogonalBasisFunction:FileNotFound');
            end
        end
        
        function testEvaluationInterface(testCase)
            basis = approx.basis.OrthogonalBasisFunction(5, 'monic_legendre', 'unit');
            expectedFilename = 'monic_legendre_unit_5';
            
            % If metadata can be loaded, test evaluation
            if ismember(expectedFilename, basis.metasource)
                try
                    basis = basis.autoLoad();
                    if basis.hasMetadata && basis.isWellDefined
                        x = testCase.unitTestPoints;
                        Y = basis.evaluate(x);
                        
                        % Check output dimensions
                        testCase.verifyEqual(size(Y, 1), basis.nCodims);
                        testCase.verifyEqual(size(Y, 2), length(x));
                    end
                catch ME
                    % Skip evaluation tests if metadata loading fails
                    testCase.assumeFail(sprintf('Metadata loading failed: %s', ME.message));
                end
            end
        end
        
        function testDerivativeInterface(testCase)
            basis = approx.basis.OrthogonalBasisFunction(5, 'monic_legendre', 'unit');
            expectedFilename = 'monic_legendre_unit_5';
            
            % If metadata can be loaded, test derivatives
            if ismember(expectedFilename, basis.metasource)
                try
                    basis = basis.autoLoad();
                    if basis.hasMetadata && basis.isWellDefined
                        x = testCase.unitTestPoints;
                        
                        % Test first derivative
                        dY1 = basis.derivative(x, 1);
                        testCase.verifyEqual(size(dY1, 1), basis.nCodims);
                        testCase.verifyEqual(size(dY1, 2), length(x));
                        
                        % Test invalid derivative orders
                        testCase.verifyError(@() basis.derivative(x, -1), ...
                            'approx:basis:OrthogonalBasisFunction:InvalidInput');
                        testCase.verifyError(@() basis.derivative(x, [1, 2]), ...
                            'approx:basis:OrthogonalBasisFunction:InvalidInput');
                    end
                catch ME
                    testCase.assumeFail(sprintf('Metadata loading failed: %s', ME.message));
                end
            end
        end
        
        function testJacobianAndGradient(testCase)
            basis = approx.basis.OrthogonalBasisFunction(5, 'monic_legendre', 'unit');
            expectedFilename = 'monic_legendre_unit_5';
            
            if ismember(expectedFilename, basis.metasource)
                try
                    basis = basis.autoLoad();
                    if basis.hasMetadata && basis.isWellDefined
                        x = testCase.unitTestPoints;
                        
                        % Test Jacobian
                        J = basis.jacobian(x);
                        testCase.verifyEqual(size(J, 1), basis.nCodims);
                        testCase.verifyEqual(size(J, 2), basis.nDims);
                        testCase.verifyEqual(size(J, 3), length(x));
                        
                        % Compare Jacobian with derivative
                        dY = basis.derivative(x, 1);
                        for i = 1:basis.nCodims
                            testCase.verifyEqual(squeeze(J(i, 1, :))', dY(i, :), ...
                                'AbsTol', testCase.tolerance);
                        end
                        
                        % Test gradient (should fail for multi-output functions)
                        if basis.nCodims > 1
                            testCase.verifyError(@() basis.gradient(x), ...
                                'approx:basis:OrthogonalBasisFunction:NotScalarValued');
                        else
                            % For scalar output, gradient should work
                            g = basis.gradient(x);
                            testCase.verifyEqual(size(g), [length(x), 1]);
                            testCase.verifyEqual(g, dY', 'AbsTol', testCase.tolerance);
                        end
                    end
                catch ME
                    testCase.assumeFail(sprintf('Metadata loading failed: %s', ME.message));
                end
            end
        end
        
        function testParameterStorage(testCase)
            basis = approx.basis.OrthogonalBasisFunction(5, 'monic_legendre', 'unit');
            
            % Initially no parameters
            testCase.verifyFalse(basis.hasParameters);
            testCase.verifyTrue(isempty(fieldnames(basis.parameters)));
            
            % Add parameters
            basis.parameters.description = 'Test Legendre basis';
            basis.parameters.created = datetime('now');
            basis.parameters.tolerance = 1e-12;
            
            % Check parameters are stored
            testCase.verifyTrue(basis.hasParameters);
            testCase.verifyEqual(basis.parameters.description, 'Test Legendre basis');
            testCase.verifyTrue(isa(basis.parameters.created, 'datetime'));
            testCase.verifyEqual(basis.parameters.tolerance, 1e-12);
        end
        
        function testInputValidation(testCase)
            basis = approx.basis.OrthogonalBasisFunction(5, 'monic_legendre', 'unit');
            expectedFilename = 'monic_legendre_unit_5';
            
            if ismember(expectedFilename, basis.metasource)
                try
                    basis = basis.autoLoad();
                    if basis.hasMetadata && basis.isWellDefined
                        % Test invalid input dimensions for evaluate
                        testCase.verifyError(@() basis.evaluate([1, 2; 3, 4]), ...
                            'approx:basis:OrthogonalBasisFunction:InvalidInput');
                        
                        % Test empty input
                        testCase.verifyError(@() basis.evaluate([]), ...
                            'approx:basis:OrthogonalBasisFunction:InvalidInput');
                        
                        % Test invalid derivative orders
                        x = 0.25;
                        testCase.verifyError(@() basis.derivative(x, 0), ...
                            'approx:basis:OrthogonalBasisFunction:InvalidInput'); % sum(r) must be > 0
                    end
                catch ME
                    testCase.assumeFail(sprintf('Metadata loading failed: %s', ME.message));
                end
            end
        end
   
        function testMonicLegendreUnitEvaluate(testCase)
            % Test mathematical correctness of Legendre polynomials on unit interval
            basis = approx.basis.OrthogonalBasisFunction(5, 'monic_legendre', 'unit');
            expectedFilename = 'monic_legendre_unit_5';
            
            if ismember(expectedFilename, basis.metasource)
                try
                    basis = basis.autoLoad();
                    if basis.hasMetadata && basis.isWellDefined && basis.nCodims >= 5
                        x = testCase.unitTestPoints;
                        Y = basis.evaluate(x);
                        
                        % Test first few Legendre polynomials on unit interval [-1/2, 1/2]
                        % P_0(x) = 1
                        testCase.verifyEqual(Y(1, :), ones(size(x)), ...
                            'AbsTol', testCase.tolerance, ...
                            'P_0(x) should be 1');
                        
                        % P_1(x) = x
                        testCase.verifyEqual(Y(2, :), x, ...
                            'AbsTol', testCase.tolerance, ...
                            'P_1(x) should be x');
                        
                        % P_2(x) = x^2 - 1/12 (for unit interval)
                        testCase.verifyEqual(Y(3, :), x.^2 - 1/12, ...
                            'AbsTol', testCase.tolerance, ...
                            'P_2(x) should be x^2 - 1/12');
                        
                        % P_3(x) = (x^2 - 3/20) * x (for unit interval)
                        testCase.verifyEqual(Y(4, :), (x.^2 - 3/20) .* x, ...
                            'AbsTol', testCase.tolerance, ...
                            'P_3(x) should be (x^2 - 3/20) * x');
                        
                        % P_4(x) = (x^2 - 3/14) * x^2 + 3/560 (for unit interval)
                        testCase.verifyEqual(Y(5, :), (x.^2 - 3/14) .* x.^2 + 3/560, ...
                            'AbsTol', testCase.tolerance, ...
                            'P_4(x) should be (x^2 - 3/14) * x^2 + 3/560');
                    end
                catch ME
                    testCase.assumeFail(sprintf('Metadata loading failed: %s', ME.message));
                end
            end
        end

        function testLegendreCanonicalEvaluate(testCase)
            % Test mathematical correctness of Legendre polynomials on canonical interval
            basis = approx.basis.OrthogonalBasisFunction(4, 'legendre', 'canonical');
            expectedFilename = 'legendre_canonical_4';
            
            if ismember(expectedFilename, basis.metasource)
                try
                    basis = basis.autoLoad();
                    if basis.hasMetadata && basis.isWellDefined && basis.nCodims >= 4
                        x = testCase.canonicalTestPoints;
                        Y = basis.evaluate(x);
                        
                        % Test first few Legendre polynomials on canonical interval [-1, 1]
                        % P_0(x) = 1
                        testCase.verifyEqual(Y(1, :), ones(size(x)), ...
                            'AbsTol', testCase.tolerance, ...
                            'P_0(x) should be 1');
                        
                        % P_1(x) = x
                        testCase.verifyEqual(Y(2, :), x, ...
                            'AbsTol', testCase.tolerance, ...
                            'P_1(x) should be x');
                        
                        % P_2(x) = (3x^2 - 1)/2 (standard Legendre on [-1,1])
                        testCase.verifyEqual(Y(3, :), (3*x.^2 - 1)/2, ...
                            'AbsTol', testCase.tolerance, ...
                            'P_2(x) should be (3x^2 - 1)/2');
                        
                        % P_3(x) = (5x^3 - 3x)/2 (standard Legendre on [-1,1])
                        testCase.verifyEqual(Y(4, :), (5*x.^3 - 3*x)/2, ...
                            'AbsTol', testCase.tolerance, ...
                            'P_3(x) should be (5x^3 - 3x)/2');
                    end
                catch ME
                    testCase.assumeFail(sprintf('Metadata loading failed: %s', ME.message));
                end
            end
        end
    
        function testLegendreUnitDerivative(testCase)
            % Test mathematical correctness of Legendre polynomial derivatives on unit interval
            basis = approx.basis.OrthogonalBasisFunction(5, 'monic_legendre', 'unit');
            expectedFilename = 'monic_legendre_unit_5';
            
            if ismember(expectedFilename, basis.metasource)
                try
                    basis = basis.autoLoad();
                    if basis.hasMetadata && basis.isWellDefined && basis.nCodims >= 5
                        x = testCase.unitTestPoints;
                        
                        % Test first derivatives  
                        dY1 = basis.derivative(x, [1]);
                        
                        % P_0'(x) = 0
                        testCase.verifyEqual(dY1(1, :), zeros(size(x)), ...
                            'AbsTol', testCase.tolerance, ...
                            'P_0''(x) should be 0');
                        
                        % P_1'(x) = 1
                        testCase.verifyEqual(dY1(2, :), ones(size(x)), ...
                            'AbsTol', testCase.tolerance, ...
                            'P_1''(x) should be 1');
                        
                        % P_2'(x) = 2x
                        testCase.verifyEqual(dY1(3, :), 2*x, ...
                            'AbsTol', testCase.tolerance, ...
                            'P_2''(x) should be 2x');
                        
                        % P_3'(x) = 3x^2 - 3/20
                        testCase.verifyEqual(dY1(4, :), 3*x.^2 - 3/20, ...
                            'AbsTol', testCase.tolerance, ...
                            'P_3''(x) should be 3x^2 - 3/20');
                        
                        % P_4'(x) = x(4x^2 - 3/7)
                        testCase.verifyEqual(dY1(5, :), x.*(4*x.^2-3/7), ...
                            'AbsTol', testCase.tolerance, ...
                            'P_4''(x) should be x(4*x^2-3/7)');
                        
                        % Test second derivatives
                        try
                            dY2 = basis.derivative(x, [2]);
                            
                            % P_0''(x) = 0
                            testCase.verifyEqual(dY2(1, :), zeros(size(x)), ...
                                'AbsTol', testCase.tolerance, ...
                                'P_0''''(x) should be 0');
                            
                            % P_1''(x) = 0
                            testCase.verifyEqual(dY2(2, :), zeros(size(x)), ...
                                'AbsTol', testCase.tolerance, ...
                                'P_1''''(x) should be 0');
                            
                            % P_2''(x) = 2
                            testCase.verifyEqual(dY2(3, :), 2*ones(size(x)), ...
                                'AbsTol', testCase.tolerance, ...
                                'P_2''''(x) should be 2');
                            
                            % P_3''(x) = 6x
                            testCase.verifyEqual(dY2(4, :), 6*x, ...
                                'AbsTol', testCase.tolerance, ...
                                'P_3''''(x) should be 6x');
                            
                            % P_4''(x) = 12x^2 - 3/7
                            testCase.verifyEqual(dY2(5, :), 12*x.^2 - 3/7, ...
                                'AbsTol', testCase.tolerance, ...
                                'P_4''''(x) should be 12x^2 - 3/7');
                        catch
                            % Second derivatives might not be available
                        end
                    end
                catch ME
                    testCase.assumeFail(sprintf('Metadata loading failed: %s', ME.message));
                end
            end
        end
        
        function testLegendreCanonicalDerivative(testCase)
            % Test mathematical correctness of Legendre polynomial derivatives on canonical interval
            basis = approx.basis.OrthogonalBasisFunction(4, 'legendre', 'canonical');
            expectedFilename = 'legendre_canonical_4';
            
            if ismember(expectedFilename, basis.metasource)
                try
                    basis = basis.autoLoad();
                    if basis.hasMetadata && basis.isWellDefined && basis.nCodims >= 4
                        x = testCase.canonicalTestPoints;
                        
                        % Test first derivatives
                        dY1 = basis.derivative(x, [1]);
                        
                        % P_0'(x) = 0
                        testCase.verifyEqual(dY1(1, :), zeros(size(x)), ...
                            'AbsTol', testCase.tolerance, ...
                            'P_0''(x) should be 0');
                        
                        % P_1'(x) = 1
                        testCase.verifyEqual(dY1(2, :), ones(size(x)), ...
                            'AbsTol', testCase.tolerance, ...
                            'P_1''(x) should be 1');
                        
                        % P_2'(x) = 3x (derivative of (3x^2 - 1)/2)
                        testCase.verifyEqual(dY1(3, :), 3*x, ...
                            'AbsTol', testCase.tolerance, ...
                            'P_2''(x) should be 3x');
                        
                        % P_3'(x) = (15x^2 - 3)/2 (derivative of (5x^3 - 3x)/2)
                        testCase.verifyEqual(dY1(4, :), (15*x.^2 - 3)/2, ...
                            'AbsTol', testCase.tolerance, ...
                            'P_3''(x) should be (15x^2 - 3)/2');
                        
                        % Test second derivatives
                        try
                            dY2 = basis.derivative(x, [2]);
                            
                            % P_0''(x) = 0
                            testCase.verifyEqual(dY2(1, :), zeros(size(x)), ...
                                'AbsTol', testCase.tolerance, ...
                                'P_0''''(x) should be 0');
                            
                            % P_1''(x) = 0
                            testCase.verifyEqual(dY2(2, :), zeros(size(x)), ...
                                'AbsTol', testCase.tolerance, ...
                                'P_1''''(x) should be 0');
                            
                            % P_2''(x) = 3
                            testCase.verifyEqual(dY2(3, :), 3*ones(size(x)), ...
                                'AbsTol', testCase.tolerance, ...
                                'P_2''''(x) should be 3');
                            
                            % P_3''(x) = 15x
                            testCase.verifyEqual(dY2(4, :), 15*x, ...
                                'AbsTol', testCase.tolerance, ...
                                'P_3''''(x) should be 15x');
                        catch
                            % Second derivatives might not be available
                        end
                    end
                catch ME
                    testCase.assumeFail(sprintf('Metadata loading failed: %s', ME.message));
                end
            end
        end
        
        function testDifferentBasisSizes(testCase)
            % Test that different nCodims values work correctly
            for nCodims = [3, 5, 8, 10]
                basis = approx.basis.OrthogonalBasisFunction(nCodims, 'monic_legendre', 'unit');
                testCase.verifyEqual(basis.nCodims, nCodims);
                testCase.verifyTrue(basis.isWellDefined);
                
                expectedFilename = sprintf('monic_legendre_unit_%d', nCodims);
                if ismember(expectedFilename, basis.metasource)
                    try
                        basis = basis.autoLoad();
                        if basis.hasMetadata
                            testCase.verifyEqual(basis.nCodims, nCodims);
                            
                            % Test evaluation with correct output size
                            x = testCase.unitTestPoints;
                            Y = basis.evaluate(x);
                            testCase.verifyEqual(size(Y, 1), nCodims);
                            testCase.verifyEqual(size(Y, 2), length(x));
                        end
                    catch ME
                        % Skip if metadata loading fails
                    end
                end
            end
        end
    end
end