classdef TestMergeResults < matlab.unittest.TestCase
    methods (TestClassSetup)
        function enableActualHelper(testCase)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'toolbox', 'examples', 'plasma')));
        end
    end

    methods (Test)
        function allGroupingModesRetainRowsParametersAndHistoricalNames(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            input = fullfile(pwd, 'input');
            mkdir(input);
            for epsilon = [1, 0.01]
                for tau0 = [0.1, 1]
                    marker = 100*epsilon+tau0;
                    data = table([2; 1], [marker; marker+0.5], ...
                        'VariableNames', {'row', 'value'});
                    name = sprintf('LDI_FDGHS2_mass_epsilon%.0e_tau0%.0e.csv', epsilon, tau0);
                    writetable(data, fullfile(input, name));
                end
            end
            % Decimal/scientific mix is a separate metric, not silently ignored.
            decimal = table([3; 4], [7; 8], 'VariableNames', {'row', 'value'});
            writetable(decimal, fullfile(input, 'LDI_FDGHS2_energy_epsilon0.25_tau0.5E+0.csv'));
            writetable(decimal, fullfile(input, 'prefix_LDI_FDGHS2_noise_epsilon1_tau01.csv'));
            writetable(decimal, fullfile(input, 'LDI_FDGHS2_noise_epsilon1_tau01_tail.csv'));

            tauDirectory = fullfile(pwd, 'by-tau');
            epsilonDirectory = fullfile(pwd, 'by-epsilon');
            metricDirectory = fullfile(pwd, 'by-metric');
            mergeResults(input, tauDirectory, 'tau0');
            mergeResults(input, epsilonDirectory, 'epsilon');
            mergeResults(input, metricDirectory, 'metric');

            for tau0 = [0.1, 1]
                file = sprintf('LDI_FDGHS2_mass_tau0%.0e_merged.csv', tau0);
                actual = readtable(fullfile(tauDirectory, file));
                testCase.verifyEqual(actual.Properties.VariableNames, {'epsilon', 'row', 'value'});
                testCase.verifyEqual(actual.epsilon, [1; 1; 0.01; 0.01], 'AbsTol', eps);
                testCase.verifyEqual(actual.row, [2; 1; 2; 1]);
                testCase.verifyEqual(actual.value, [100+tau0; 100.5+tau0; 1+tau0; 1.5+tau0], ...
                    'AbsTol', 1e-12);
            end
            for epsilon = [0.01, 1]
                file = sprintf('LDI_FDGHS2_mass_epsilon%.0e_merged.csv', epsilon);
                actual = readtable(fullfile(epsilonDirectory, file));
                testCase.verifyEqual(actual.Properties.VariableNames, {'tau0', 'row', 'value'});
                testCase.verifyEqual(actual.tau0, [0.1; 0.1; 1; 1], 'AbsTol', eps);
                testCase.verifyEqual(actual.row, [2; 1; 2; 1]);
                testCase.verifyEqual(actual.value, ...
                    100*epsilon+[0.1; 0.6; 1; 1.5], 'AbsTol', 1e-12);
            end
            actual = readtable(fullfile(metricDirectory, 'LDI_FDGHS2_mass_all_merged.csv'));
            testCase.verifyEqual(actual.Properties.VariableNames, {'epsilon', 'tau0', 'row', 'value'});
            testCase.verifyEqual(actual.epsilon, [0.01*ones(4, 1); ones(4, 1)], 'AbsTol', eps);
            testCase.verifyEqual(actual.tau0, repmat([0.1; 0.1; 1; 1], 2, 1), 'AbsTol', eps);
            testCase.verifyEqual(actual.row, repmat([2; 1], 4, 1));
            testCase.verifyEqual(actual.value, [1.1; 1.6; 2; 2.5; 100.1; 100.6; 101; 101.5], ...
                'AbsTol', 1e-12);
            actual = readtable(fullfile(metricDirectory, 'LDI_FDGHS2_energy_all_merged.csv'));
            testCase.verifyEqual(actual.epsilon, [0.25; 0.25]);
            testCase.verifyEqual(actual.tau0, [0.5; 0.5]);
            testCase.verifyEqual(actual.row, [3; 4]);
            testCase.verifyNumElements(dir(fullfile(metricDirectory, '*.csv')), 2);

            for group = ["tau0", "epsilon", "metric"]
                directory = fullfile(pwd, "by-"+group+"-overwrite");
                mergeResults(input, directory, group);
                testCase.verifyError(@() mergeResults(input, directory, group), ...
                    'mdemat:examples:FileExists');
                mergeResults(input, directory, group, Overwrite=true);
            end
        end
    end
end
