classdef TestSymbolicSphericalHarmonics < matlab.unittest.TestCase
    methods (Test, TestTags = {'Symbolic'})
        function symbolicModesMatchNumericConvention(testCase)
            testCase.assertTrue(logical(license('test', 'Symbolic_Toolbox')), ...
                'Symbolic Math Toolbox is needed to regenerate symbolic expressions.');
            [points, weights] = TestSymbolicSphericalHarmonics.sphereRule();
            variables = sym('z', [2, 1], 'real');
            for normalized = [false, true]
                symbolic = mdemat.common.symbolic.SphericalHarmonicBasis(variables, 2, ...
                    isNormalized=normalized);
                symbolic.compile();
                values = zeros(9, size(points, 2));
                for i = 1:9
                    values(i, :) = symbolic.Handles{i, 1}(points(1, :), points(2, :));
                end
                direct = mdemat.common.spectral.SphericalHarmonicBasis(2, ...
                    isNormalized=normalized);
                testCase.verifyEqual(values, direct.eval(points), 'AbsTol', 2e-12);
                if normalized
                    testCase.verifyEqual((values .* weights) * values.', eye(9), ...
                        'AbsTol', 2e-12);
                end
            end
        end
    end

    methods (Static, Access = private)
        function [points, weights] = sphereRule()
            n = 12;
            k = (1:n-1).';
            offDiagonal = k ./ sqrt(4*k.^2-1);
            [vectors, values] = eig(diag(offDiagonal, 1) + diag(offDiagonal, -1));
            mu = diag(values);
            muWeights = 2 * vectors(1, :).^2;
            phi = (0:31) * (2*pi/32);
            [muGrid, phiGrid] = ndgrid(mu, phi);
            points = [muGrid(:).'; phiGrid(:).'];
            weights = kron(ones(1, numel(phi))*(2*pi/numel(phi)), muWeights);
        end
    end
end
