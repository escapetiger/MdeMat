classdef TestLinearAlgebra < matlab.unittest.TestCase
    methods (Test)
        function directResidualAndMultipleRhs(test)
            A = [4, 1; 1, 3];
            b = [1, 0; 2, 0];
            solver = mdemat.common.linalg.LinearSolver();
            [x, report] = solver.solve(A=A, b=b);
            test.verifyEqual(A * x, b, AbsTol=1e-13);
            test.verifyEqual(report.relativeResidual, norm(A*x-b, 'fro') / norm(b, 'fro'));
            test.verifyTrue(report.converged);
            solver.setLhs([5, 1; 1, 3]);
            test.verifyFalse(solver.IsPrecomputed);
        end

        function iterativeComplexMultipleRhs(test)
            A = sparse([5, 1i; -1i, 4]);
            b = [1+2i, 3; 4, 5-1i];
            solver = mdemat.common.linalg.LinearSolver(nnzTh=1, tol=1e-12);
            [x, report] = solver.solve(A=A, b=b);
            test.verifyLessThan(norm(A*x-b, 'fro') / norm(b, 'fro'), 1e-12);
            test.verifyTrue(report.converged);
            test.verifySize(x, [2, 2]);
        end

        function conditionEstimatorPreservesRng(test)
            old = rng;
            restore = onCleanup(@() rng(old)); %#ok<NASGU>
            A = spdiags([-ones(8, 1), 4*ones(8, 1), -ones(8, 1)], -1:1, 8, 8);
            before = rng;
            solver = mdemat.common.linalg.LinearSolver(nnzTh=1, condTh=1, tol=1e-12);
            x = solver.solve(A=A, b=ones(8, 1));
            test.verifyEqual(rng, before);
            test.verifyLessThan(norm(A*x-1), 1e-11);
        end

        function hashDistinguishesRepresentations(test)
            A = reshape(1:6, 2, 3);
            h = mdemat.tools.base.hash(A);
            test.verifyEqual(mdemat.tools.base.hash(A), h);
            test.verifyNotEqual(mdemat.tools.base.hash(A.'), h);
            test.verifyNotEqual(mdemat.tools.base.hash(single(A)), h);
            test.verifyNotEqual(mdemat.tools.base.hash(sparse(A)), h);
            test.verifyNotEqual(mdemat.tools.base.hash(A+1i), h);
            test.verifyEqual(mdemat.tools.base.hash(sparse(A)), mdemat.tools.base.hash(sparse(A)));
        end

        function rejectsInvalidSystem(test)
            solver = mdemat.common.linalg.LinearSolver();
            test.verifyError(@() solver.solve(A=eye(2), b=1), 'mdemat:linalg:InvalidRhs');
            test.verifyError(@() solver.solve(A=eye(2), b=[NaN; 1]), ...
                'mdemat:linalg:NonfiniteSystem');
        end
    end
end
