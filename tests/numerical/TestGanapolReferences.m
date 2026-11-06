classdef TestGanapolReferences < matlab.unittest.TestCase
    methods (TestClassSetup)
        function enableIndependentReference(test)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            test.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'benchmarks', 'ganapol')));
            test.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'benchmarks', 'ganapol', 'utils')));
        end
    end
    methods (Test)
        function inverseTransformAtOriginIsIntegratedOnce(test)
            for acceleration = [0, 0.45, 1]
                for order = [0, 4, 12]
                    actual = cifft(@(k) pi./(k.^2+4), 0, order, acceleration);
                    test.verifyEqual(actual, pi/4, AbsTol=2e-12);
                end
            end
            for x = [-1, 1]
                test.verifyEqual(cifft(@(k) pi./(k.^2+4), x, 20), ...
                    pi*exp(-2*abs(x))/4, AbsTol=3e-5);
            end
        end
        function combinatoricsQuadratureAndExtrapolation(test)
            test.verifyEqual(binom(2, [0, 1, 2, 3]), [1, 2, 1, 0]);
            test.verifyEqual(binom(0:3, 0:3), [1,0,0,0;1,1,0,0;1,2,1,0;1,3,3,1]);
            test.verifyEqual(euler_knopp(@(k) 0.8.^k, 40, 0.5), 1/1.8, AbsTol=1e-11);
            [x, w] = gausslegd(1, 2, 5, false);
            test.verifyEqual(x, 3.5); test.verifyEqual(w, 3);
            test.verifyEqual(gaussint(@(x) x^5, -2, 3, 3), (3^6-2^6)/6, AbsTol=2e-12);
            [value, report] = romberg(@sin, 0, pi);
            test.verifyEqual(value, 2, AbsTol=1e-10); test.verifyTrue(report.Converged);
        end
        function originalFourierDriverMatchesAnalyticTransforms(test)
            originalPath = path;
            result = ex_cifft();
            test.verifySize(result.Computed, [1, 9]);
            test.verifyTrue(all(isfinite(result.Computed)));
            exact = [pi*exp(-0.1)/2, atan(2), pi/(4*cosh(pi)), ...
                pi*exp(-2)/4, pi*exp(-1)*(cos(1)+sin(1))];
            test.verifyEqual(result.Computed([1, 3, 4, 7, 8]), exact, AbsTol=4e-6);
            test.verifyEqual(result.UnconvergedComparisons, [5, 6]);
            test.verifyEqual(path, originalPath);
        end
        function collisionlessPlaneSphereCylinderHaveCorrectNormalization(test)
            x = [0.2, 0.7, 2];
            plane = ganapolReference(x);
            test.verifyEqual(plane, expint(x)/2, AbsTol=1e-14);
            for scattering = [0, 0.3, 0.9]
                test.verifyEqual(ganapolReference(0, ScatteringRatio=scattering), Inf);
            end
            sphere = ganapolReference(x, Geometry="sphere");
            test.verifyEqual(4*pi*x.^2.*sphere, exp(-x), AbsTol=1e-14);
            cylinder = ganapolReference(x, Geometry="cylinder", ProjectionPoints=64);
            expected = arrayfun(@(r) integral(@(theta) exp(-r./cos(theta)), ...
                0, pi/2, 'AbsTol', 1e-12)/(2*pi*r), x);
            test.verifyEqual(cylinder, expected, AbsTol=2e-11);
            test.verifyEqual(ganapolReference(0, Geometry="cylinder"), Inf);
        end
        function beamAndStepRecoverIndependentCollisionlessSolutions(test)
            x = [-1, -0.2, 0, 0.2, 1];
            beam = ganapolReference(x, Source="beam", BeamDirection=0.5);
            test.verifyEqual(beam, [0,0,1,2*exp(-0.4),2*exp(-2)], AbsTol=1e-14);
            step = ganapolReference(x, Source="step");
            test.verifyEqual(step+fliplr(step), ones(size(step)), AbsTol=2e-14);
            derivative = (ganapolReference(0.7+1e-5, Source="step") ...
                -ganapolReference(0.7-1e-5, Source="step"))/(2e-5);
            test.verifyEqual(derivative, expint(0.7)/2, RelTol=1e-9);
        end
        function scatteringReferenceConverges(test)
            a = ganapolReference(0.8, ScatteringRatio=0.3, AngularPoints=32, FourierOrder=18);
            b = ganapolReference(0.8, ScatteringRatio=0.3, AngularPoints=48, FourierOrder=22);
            test.verifyGreaterThan(a, expint(0.8)/2);
            test.verifyEqual(a, b, AbsTol=8e-5);
        end
        function allDriversMatchNumericHistoricalFixturesWithoutGraphics(test)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            folder = test.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture());
            figures = findall(groot, 'Type', 'figure');
            originalPath = path;
            drivers = { ...
                'plane_isotropic_source_vc', 'plane_beam_source_vc', ...
                'plane_beam_source_vf', 'plane_beam_source_vmu', ...
                'plane_step_source_vc', 'sphere_isotropic_source_vc', ...
                'cylinder_isotropic_source_vc'};
            for index = 1:numel(drivers)
                name = drivers{index};
                output = fullfile(folder.Folder, [name, '.json']);
                result = feval(['ex_og1d_', name], Plot=false, MaxPoints=2, OutputFile=output);
                test.verifySize(result.Flux, [numel(result.Flux(:, 1)), 2]);
                test.verifyTrue(all(isfinite(result.Flux(:, 2))));
                test.verifyEqual(mdemat.tools.io.readRecord(output), result);
                if contains(name, 'cylinder')
                    test.verifyEqual(result.Flux(:, 1), Inf);
                    continue; % The original repository contains no cylinder MAT fixture.
                end
                filename = ['data_og1d_', name, '.mat'];
                if strcmp(name, 'plane_beam_source_vf')
                    filename = 'data_og1d_plane_beam_souce_vf.mat';
                end
                % Historic fixtures also contain a saved function handle g.
                % Load only named numerical data, never deserialize that callback.
                variables = {'x', 'phi0', 'c', 'fs', 'n'};
                if strcmp(name, 'plane_beam_source_vmu')
                    variables{end+1} = 'mu0';
                end
                data = load(fullfile(root, 'benchmarks', 'ganapol', filename), variables{:});
                test.verifyEqual(result.Points, data.x(1:2));
                test.verifyEqual(result.ScatteringRatio, data.c);
                test.verifyEqual(result.Phase, data.fs);
                test.verifyEqual(result.AngularPoints, data.n);
                % At zero the old inverse transform repeated its nonoscillatory
                % integral, and the sphere fixture contains NaN. Positive
                % radius data retain the original full parameter sweep.
                test.verifyEqual(result.Flux(:, 2), data.phi0(:, 2), ...
                    RelTol=1e-10, AbsTol=2e-11);
                if isfield(data, 'mu0')
                    test.verifyEqual(result.BeamDirection, data.mu0);
                end
                if contains(name, 'sphere')
                    test.verifyTrue(all(isnan(data.phi0(:, 1))));
                    test.verifyEqual(result.Flux(:, 1), Inf(size(result.Flux(:, 1))));
                elseif strcmp(name, 'plane_isotropic_source_vc')
                    test.verifyEqual(result.Flux(:, 1), Inf(size(result.Flux(:, 1))));
                end
            end
            test.verifyEqual(findall(groot, 'Type', 'figure'), figures);
            test.verifyEqual(path, originalPath);
        end
        function scatteringCylinderMatchesIndependentContinuousSpectrum(test)
            radius = 0.8;
            for scattering = [0.3, 0.7]
                expected = isotropicCylinderSpectrum(radius, scattering);
                coarse = ganapolReference(radius, Geometry="cylinder", Phase=1, ...
                    ScatteringRatio=scattering, AngularPoints=24, FourierOrder=14, ProjectionPoints=20);
                fine = ganapolReference(radius, Geometry="cylinder", Phase=1, ...
                    ScatteringRatio=scattering, AngularPoints=64, FourierOrder=28, ProjectionPoints=48);
                test.verifyGreaterThan(fine, ganapolReference(radius, Geometry="cylinder"));
                test.verifyLessThan(abs(fine-expected), 0.05*abs(coarse-expected));
                test.verifyEqual(fine, expected, AbsTol=5e-7);
            end
        end
    end
end

function flux = isotropicCylinderSpectrum(radius, scattering)
% Independent continuous-angle inversion of T(k)/(1-c*T(k)), T=atan(k)/k.
% The pole k=i*s obeys c*atanh(s)/s=1. The two sides of the branch cut
% give D(mu)=(1-c*mu*atanh(mu))^2+(pi*c*mu/2)^2. Projecting each spherical
% exp(-s*r)/r term along the line produces 2*K0(s*r), without Fourier
% segmentation, a discrete-ordinates system, or the reference's z quadrature.
s = fzero(@(s) 1-scattering*atanh(s)/s, [sqrt(eps), 1-eps]);
residue = s^2/(scattering^2*(s/(1-s^2)-atanh(s)));
continuous = @(mu) besselk(0, radius./mu)./(mu.^2 .* ...
    ((1-scattering*mu.*atanh(mu)).^2 + (pi*scattering*mu/2).^2));
flux = residue*s*besselk(0, s*radius)/pi + ...
    integral(continuous, 0, 1, AbsTol=1e-11, RelTol=1e-10)/(2*pi);
end
