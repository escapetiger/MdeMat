classdef TestNewtonDiagnostics < matlab.unittest.TestCase
    methods (Test)
        function restartedGmresReportsTheReturnedIterateAndPreservesRawCounters(testCase)
            matrix = diag([1, 2, 4, 8]);
            rhs = ones(4, 1);
            restart = 2;
            [reference, flag, ~, raw, residuals] = gmres(matrix, rhs, restart, 1e-10, 100);
            testCase.assertEqual(flag, 0);
            testCase.assertGreaterThan(raw(1), 1);
            optimizer = mdemat.common.nonlinear.JfnkOptimizer(maxIter=1, ...
                useMatrix=true, jacFun=@(~) matrix, fTol=1e-7, ...
                krylovOpts=struct('restart', restart, 'tol', 1e-10, 'maxIt', 100));
            [value, report] = optimizer.solve(@(x) matrix*x-rhs, zeros(4, 1));
            testCase.verifyEqual(value, reference, 'AbsTol', 1e-13);
            testCase.verifyEqual(report.krylovStats.rawIterations{1}, raw);
            % A successful GMRES residual trace has the initial residual
            % followed by one value for each completed inner iteration.
            testCase.verifyEqual(report.krylovStats.linearSteps(1), numel(residuals)-1);
            testCase.verifyGreaterThan(report.krylovStats.linearSteps(1), max(raw));
        end
    end
end
