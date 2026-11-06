classdef TestNumericBasisResources < matlab.unittest.TestCase
    % Runtime resource schema, numeric identities and expression boundaries.

    methods (Test)
        function interpolationNodesAreActualGaussAndLobattoRoots(testCase)
            ordering = {[1], [1,2], [2,1,3], [1,4,2,3], [3,2,4,1,5], [1,6,2,5,3,4]};
            for count = 1:6
                legendre = legendrePolynomial(count);
                gaussian = sort(real(roots(legendre)))/2;
                gaussian = gaussian(ordering{count});
                for family = ["legendre", "lobatto"]
                    if family == "legendre"
                        nodes = gaussian;
                    elseif count == 1
                        nodes = 0;
                    else
                        nodes = [-0.5; sort(real(roots(polyder(legendrePolynomial(count-1)))))/2; 0.5];
                    end
                    name = sprintf('unit_gauss_%s_lagrange_%d', family, count);
                    basis = mdemat.common.approximation.CompiledFunction(1, count, name);
                    testCase.verifyEqual(basis.eval(nodes.'), eye(count), ...
                        sprintf('%s interpolation identity at independently computed roots.', name), ...
                        'AbsTol', 8e-12);
                end
                element = mdemat.common.fem.BH1OrthotopeElement.nodal(1,count);
                weights = 1./((1-4*gaussian.^2).*polyval(polyder(legendre),2*gaussian).^2);
                testCase.verifyEqual(element.Approximator.Mass, diag(weights), 'AbsTol', 1e-11);
            end
        end
        function bundledDataContainsOnlyPortableValues(testCase)
            directory = TestNumericBasisResources.resourceDirectory();
            files = dir(fullfile(directory, '*.mat'));
            testCase.verifyNumElements(files, 18);
            for index = 1:numel(files)
                loaded = load(fullfile(directory, files(index).name), 'metadata');
                metadata = loaded.metadata;
                testCase.verifyTrue(TestNumericBasisResources.isPortable(metadata), ...
                    files(index).name);
                basis = mdemat.common.approximation.CompiledFunction( ...
                    metadata.nDims, metadata.nCodims, metadata.basisName);
                testCase.verifyTrue(contains(basis.Metadata.resourceFile, ...
                    fullfile('resources', 'basis')));
                testCase.verifyEqual(basis.Metadata.derivativeOrders, ...
                    0:size(metadata.functionText, 2)-1);
            end
        end

        function monicLegendreValuesAndDerivatives(testCase)
            basis = mdemat.common.approximation.CompiledFunction(1, 4, 'monic_unit_legendre');
            x = linspace(-0.49, 0.49, 23);
            testCase.verifyEqual(basis.eval(x), ...
                [ones(size(x)); x; x.^2-1/12; x.^3-3*x/20], 'AbsTol', 2e-14);
            testCase.verifyEqual(basis.diff(x, 1), ...
                [zeros(size(x)); ones(size(x)); 2*x; 3*x.^2-3/20], 'AbsTol', 2e-14);
            testCase.verifyEqual(basis.diff(x, 3), ...
                [zeros(3, numel(x)); 6*ones(size(x))], 'AbsTol', 2e-14);
        end

        function interpolationPreservesConstants(testCase)
            x = linspace(-0.49, 0.49, 23);
            for family = {'legendre', 'lobatto'}
                for degreeCount = 1:6
                    name = sprintf('unit_gauss_%s_lagrange_%d', family{1}, degreeCount);
                    basis = mdemat.common.approximation.CompiledFunction(1, degreeCount, name);
                    testCase.verifyEqual(sum(basis.eval(x), 1), ones(size(x)), ...
                        'AbsTol', 2e-11);
                    testCase.verifyEqual(sum(basis.diff(x, 1), 1), zeros(size(x)), ...
                        'AbsTol', 2e-10);
                end
            end
        end

        function undeclaredResourceNamesAreRejected(testCase)
            for name = {'../normal_legendre', 'normal_legendre/../normal_legendre', ...
                    'unknown_basis', 'unit_gauss_legendre_lagrange_7'}
                testCase.verifyError(@() mdemat.common.approximation.CompiledFunction(1, 1, name{1}), ...
                    'MdeMat:BasisResource:InvalidName');
            end
            testCase.verifyError(@() mdemat.common.approximation.CompiledFunction( ...
                2, 37, 'normal_spherical_harmonic'), ...
                'MdeMat:BasisResource:CapacityExceeded');
        end

        function schemaAndArithmeticGrammarAreEnforced(testCase)
            loaded = load(fullfile(TestNumericBasisResources.resourceDirectory(), ...
                'normal_legendre.mat'), 'metadata');
            metadata = loaded.metadata;
            bad = metadata;
            bad.expressions = {'unapproved field'};
            testCase.verifyError(@() mdemat.common.approximation.CompiledFunction.compileResource(bad, 1, 1), ...
                'MdeMat:BasisResource:InvalidSchema');
            bad = metadata;
            bad.derivativeOrders = 1;
            testCase.verifyError(@() mdemat.common.approximation.CompiledFunction.compileResource(bad, 1, 1), ...
                'MdeMat:BasisResource:InvalidSchema');
            for text = {'@(x)system(x)', '@(x)sym(x)', '@(x)eval(x)', ...
                    '@(x)x(1)', '@(x)x;1', '@(x)sin(x)+y', '@(x)x.field', ...
                    '@(x)cos', '@(x)pi(1)', '@(x,y)x', '@(x)[x,x]'}
                bad = metadata;
                bad.functionText{1} = text{1};
                testCase.verifyError(@() mdemat.common.approximation.CompiledFunction.compileResource(bad, 1, 1), ...
                    'MdeMat:BasisResource:InvalidExpression', text{1});
            end
        end
    end

    methods (Static, Access = private)
        function directory = resourceDirectory()
            root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
            directory = fullfile(root, 'toolbox', 'resources', 'basis');
        end

        function result = isPortable(value)
            if isnumeric(value) || islogical(value) || ischar(value)
                result = true;
            elseif iscell(value)
                result = all(cellfun(@TestNumericBasisResources.isPortable, value), 'all');
            elseif isstruct(value)
                contents = struct2cell(value);
                result = all(cellfun(@TestNumericBasisResources.isPortable, contents), 'all');
            else
                result = false;
            end
        end
    end
end

function coefficients = legendrePolynomial(degree)
previous = 1;
if degree == 0, coefficients = previous; return; end
coefficients = [1,0];
for k = 2:degree
    next = ((2*k-1)*[coefficients,0]-(k-1)*[0,0,previous])/k;
    previous = coefficients;
    coefficients = next;
end
end
