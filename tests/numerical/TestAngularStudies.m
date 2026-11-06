classdef TestAngularStudies < matlab.unittest.TestCase
    methods (TestClassSetup)
        function enableActualStudyEntrypoints(testCase)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'toolbox', 'examples', 'radiation')));
        end
    end

    methods (Test)
        function fourierProjectionAndNodalResidualUseTheSameAngleChart(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            before = findall(groot, 'Type', 'figure');
            result = ex00_func_approx(ReferencePoints=192, ModalOrders=[2, 3], ...
                NodalDofs=[4, 8], HybridConfigurations=[2, 4; 3, 8], ...
                MakePlots=false, OutputDirectory="", Verbose=false);
            theta = result.theta;
            equilibrium = 1/(2*pi);
            lowModes = equilibrium*(1+0.5*cos(theta)+0.2*sin(theta));
            exact = lowModes+0.3*equilibrium*cos(2*theta);
            testCase.verifyEqual(result.exactValues(1, :), exact, 'AbsTol', 2e-15);
            testCase.verifyEqual(reshape(result.modalValues(1, 1, :), 1, []), ...
                lowModes, 'AbsTol', 2e-13);
            testCase.verifyEqual(reshape(result.modalValues(1, 2, :), 1, []), ...
                exact, 'AbsTol', 2e-13);
            omittedModeError = 0.3/sqrt(2+0.5^2+0.3^2+0.2^2);
            testCase.verifyEqual(result.errModal(1, 1), omittedModeError, 'AbsTol', 2e-13);
            testCase.verifyLessThan(result.errModal(1, 2), 2e-13);
            testCase.verifyEqual(sum(result.exactValues, 2)*(2*pi/192), ...
                [1; 1; 1; 1.5], 'AbsTol', 2e-13);

            % Original angular nodes retain their pi/6 shift and ordering.
            % Compare evaluated fields, never assume that RBF DOFs are samples.
            for index = 1:2
                count = result.nvList(index);
                columns = 1+mod(192/12+(0:count-1)*(192/count), 192);
                samples = reshape(result.nodalValues(:, index, :), 4, []);
                testCase.verifyEqual(samples(:, columns), ...
                    result.exactValues(:, columns), 'AbsTol', 2e-12);
                hybrid = reshape(result.hybridValues(:, index, :), 4, []);
                testCase.verifyEqual(hybrid(:, columns), ...
                    result.exactValues(:, columns), 'AbsTol', 2e-12);
            end
            testCase.verifyEqual(findall(groot, 'Type', 'figure'), before);
            testCase.verifyFalse(isfile('results.mat'));
        end

        function continuousTransportReferenceConvergesAcrossAllRetainedStudies(testCase)
            result = ex05_generalization(ReferencePoints=192, Dofs=[9, 17], ...
                MakePlots=false, OutputDirectory="", Verbose=false);
            testCase.verifySize(result.values, [7, 5, 2, 192]);
            testCase.verifyTrue(all(isfinite(result.errL2), 'all'));
            testCase.verifyNumElements(result.caseNames, 7);
            testCase.verifyNumElements(result.configurationNames, 5);
            theta = result.theta;
            exact = 1+cos(pi*(0.25-cos(theta))).*cos(pi*(-0.15-sin(theta)));
            testCase.verifyEqual(result.exactValues(1, :), exact, 'AbsTol', 2e-15);

            % Product-to-sum plus the continuous integral of exp(i*r*cos(theta))
            % gives this exact mean, independently of quadrature/basis matrices.
            exactMean = 1+cos(pi*0.25)*cos(pi*(-0.15))*besselj(0, pi*sqrt(2));
            fine = reshape(result.values(1, 1, 2, :), 1, []);
            testCase.verifyEqual(mean(fine), exactMean, 'AbsTol', 2e-12);
            relativeError = norm(fine-exact)/norm(exact);
            testCase.verifyLessThan(relativeError, 2e-3);
            testCase.verifyLessThan(relativeError, result.errL2(1, 1, 1)/20);
            testCase.verifyEqual(relativeError, result.errL2(1, 1, 2), 'AbsTol', 2e-15);

            for index = 1:3
                width = [0.20, 0.10, 0.05];
                normalizer = width(index)*sqrt(2*pi)*erf(pi/(sqrt(2)*width(index)));
                % pi/6 is the seventeenth point on this independent reference grid.
                testCase.verifyEqual(result.exactValues(4+index, 17), ...
                    1/normalizer, 'RelTol', 2e-14);
                testCase.verifyEqual(sum(result.exactValues(4+index, :))*(2*pi/192), ...
                    1, 'AbsTol', 2e-12);
            end
        end

        function controlledPlotsPreserveExistingFiguresAndOutputIsExplicit(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            visible = get(groot, 'DefaultFigureVisible');
            testCase.addTeardown(@() set(groot, 'DefaultFigureVisible', visible));
            set(groot, 'DefaultFigureVisible', 'off');
            original = findall(groot, 'Type', 'figure');
            testCase.addTeardown(@() TestAngularStudies.deleteCreatedFigures(original));
            sentinel = figure('Visible', 'off', 'Name', 'Caller figure');
            ex00_func_approx(ReferencePoints=24, ModalOrders=2, NodalDofs=4, ...
                HybridConfigurations=[2, 4], MakePlots=true, OutputDirectory="", Verbose=false);
            result = ex05_generalization(ReferencePoints=24, Dofs=3, ...
                MakePlots=true, OutputDirectory="", Verbose=false);
            testCase.verifyTrue(isgraphics(sentinel));
            testCase.verifyEqual(get(sentinel, 'Name'), 'Caller figure');
            testCase.verifySize(result.meanL2, [5, 1]);
            testCase.verifyNumElements(findall(groot, 'Type', 'figure'), numel(original)+25);
            testCase.verifyFalse(isfile('results.mat'));

            directory = fullfile(pwd, 'explicit-output');
            saved = ex05_generalization(ReferencePoints=24, Dofs=3, MakePlots=false, ...
                OutputDirectory=directory, Verbose=false);
            contents = load(fullfile(directory, 'results.mat'), 'results');
            testCase.verifyEqual(contents.results, saved);
            testCase.verifyError(@() ex05_generalization(ReferencePoints=24, Dofs=3, ...
                MakePlots=false, OutputDirectory=directory, Verbose=false), ...
                'MdeMat:AngularStudy:FileExists');
            testCase.verifyError(@() ex05_generalization(Dofs=4, MakePlots=false), ...
                'MdeMat:AngularStudy:InvalidDofs');
        end
    end

    methods (Static, Access = private)
        function deleteCreatedFigures(original)
            figures = findall(groot, 'Type', 'figure');
            delete(setdiff(figures, original));
        end
    end
end
