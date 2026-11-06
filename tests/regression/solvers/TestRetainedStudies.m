classdef TestRetainedStudies < matlab.unittest.TestCase
    %TESTRETAINEDSTUDIES Execute the original research builders at small size.
    properties (TestParameter)
        entry = {'elliptic/ex01_periodic', 'elliptic/ex02_dirichlet', ...
            'hyperbolic/ex01_periodic', 'hyperbolic/ex02_dirichlet', ...
            'hyperbolic/ex03_discontinuous', 'parabolic/ex01_periodic', ...
            'parabolic/ex02_dirichlet', 'parabolic/ex03_inhomogenous', ...
            'radiation/ex01_convergence', 'radiation/ex02_isotropic_dirichlet', ...
            'radiation/ex03_gaussian', 'radiation/ex04_beam', ...
            'plasma/ex01_convergence', 'plasma/ex01_time_convergence', ...
            'plasma/ex01_reference_convergence', 'plasma/ex01_order_diagnostics', ...
            'plasma/ex02_asymptotic_preserving', 'plasma/ex03_collision'}
    end
    methods (Test)
        function originalEntryRunsItsRealModelAndReturnsAcceptedResults(testCase, entry)
            toolboxRoot = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
            source = fullfile(toolboxRoot, 'toolbox', 'examples', entry);
            folder = fileparts(source);
            [~, name] = fileparts(source);
            originalPath = path;
            cleanup = onCleanup(@() path(originalPath)); %#ok<NASGU>
            addpath(folder, '-begin');
            clear(name);
            factory = str2func(name);
            testCase.assertEqual(which(name), [source, '.m']);
            figures = findall(groot, 'Type', 'figure');
            [state, result] = factory(Smoke=true);
            testCase.verifyNotEmpty(state);
            testCase.verifyNotEmpty(result);
            testCase.verifyTrue(isa(state, 'mdemat.data.variable.State'));
            fields = fieldnames(state.Dofs);
            for index = 1:numel(fields)
                coefficients = state.Dofs.(fields{index});
                testCase.verifyTrue(all(isfinite(coefficients(:))));
            end
            testCase.verifyEqual(findall(groot, 'Type', 'figure'), figures);
            if startsWith(entry, 'radiation/') || startsWith(entry, 'plasma/')
                directory = tempname;
                mkdir(directory);
                remove = onCleanup(@() rmdir(directory, 's')); %#ok<NASGU>
                file = fullfile(directory, 'accepted.mat');
                state.save(file, omitCallbacks=true);
                restored = mdemat.internal.execution.loadStudyState(file, state);
                testCase.verifyEqual(restored.Dofs, state.Dofs);
                testCase.verifyEqual(restored.History, state.History);
                if ismethod(state, 'density')
                    nodes = state.XDisc.Element.Volume.Nodes;
                    testCase.verifyEqual(restored.density(nodes), state.density(nodes), 'AbsTol', 2e-13);
                end
            end
            if isstruct(result) && isfield(result, 'execution')
                testCase.verifyEqual(result.execution.steps, 2);
                testCase.verifyEqual(result.times(end), 0.002, 'AbsTol', eps);
            end
            if strcmp(entry, 'radiation/ex01_convergence')
                % A velocity sequence is an aggregate over all selected
                % components, so its identity must not claim component 1.
                uName = ['U',sprintf('_%d',1:state.NVMacroDofs),'_sequence'];
                gName = ['G',sprintf('_%d',1:state.NVMicroDofs),'_sequence'];
                testCase.assertTrue(isfield(result,uName) && isfield(result,gName));
                testCase.verifyNotEmpty(result.(uName));
                testCase.verifyNotEmpty(result.(gName));
                testCase.verifyTrue(all(isfinite(result.(uName).L2_Error)));
                testCase.verifyTrue(all(isfinite(result.(gName).L2_Error)));
                testCase.verifyFalse(isfield(result, 'F1'));
            end
        end

        function discontinuousStudyPersistsOnlyNumericSamples(testCase)
            root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
            originalPath = path;
            restore = onCleanup(@() path(originalPath)); %#ok<NASGU>
            addpath(fullfile(root, 'toolbox', 'examples', 'hyperbolic'), '-begin');
            directory = tempname;
            mkdir(directory);
            remove = onCleanup(@() rmdir(directory, 's')); %#ok<NASGU>
            filename = fullfile(directory, 'sample.mat');
            [~, result] = ex03_discontinuous(Smoke=true, OutputFile=filename);
            sample = mdemat.internal.execution.readStudySamples(filename);
            testCase.verifyEqual(sample, result.Samples);
            testCase.verifyEqual(size(sample.Points, 2), size(sample.Values.U, 1));
            testCase.verifyTrue(all(isfinite(sample.Values.U(:))));
            listing = whos('-file', filename);
            testCase.verifyTrue(all(ismember({listing.class}, ...
                {'double', 'single', 'logical', 'uint8'})));
        end
    end
end
