classdef TestSphericalHarmonics < matlab.unittest.TestCase
    % The sphere measure in (mu,phi) is dmu*dphi, mu=cos(theta).
    % These checks use an independent Golub--Welsch/trapezoidal product
    % rule rather than the production angular quadratures.

    methods (Test)
        function directAndCompiledRealModesAreOrthonormal(testCase)
            [points, weights] = TestSphericalHarmonics.sphereRule();
            direct = mdemat.common.spectral.SphericalHarmonicBasis(5, isNormalized=true);
            compiled = mdemat.common.spectral.CompiledOrthogonalBasis(36, ...
                'normal_spherical_harmonic');
            values = direct.eval(points);
            compiledValues = compiled.eval(points);
            testCase.verifyEqual((values .* weights) * values.', eye(36), ...
                'AbsTol', 2e-12);
            testCase.verifyEqual((compiledValues .* weights) * compiledValues.', ...
                eye(36), 'AbsTol', 2e-12);
            testCase.verifyEqual(compiledValues, values, 'AbsTol', 2e-12);

            mu = points(1, :);
            phi = points(2, :);
            cartesian = [mu; sqrt(1-mu.^2).*cos(phi); sqrt(1-mu.^2).*sin(phi)];
            expected = [ones(1, size(points, 2))/sqrt(4*pi); ...
                sqrt(3/(4*pi))*cartesian];
            testCase.verifyEqual(values(1:4, :), expected, 'AbsTol', 2e-13);
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
