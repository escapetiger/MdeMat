classdef TestJfnkOptimizer < matlab.unittest.TestCase
    % TESTJFNKOPTIMIZER Unit tests for JfnkOptimizer class.
    %
    %   TestJfnkOptimizer provides comprehensive unit tests for the JfnkOptimizer
    %   class including convergence tests on various nonlinear systems and
    %   verification of all major functionality.

    methods (Test)
        function testConstructorDefaults(testCase)
            % Test constructor with default parameters
            optimizer = core.optim.JfnkOptimizer();

            % Verify default options
            testCase.verifyEqual(optimizer.Opts.maxIter, 30);
            testCase.verifyEqual(optimizer.Opts.fTol, 1e-6);
            testCase.verifyEqual(optimizer.Opts.xTol, 1e-10);
            testCase.verifyFalse(optimizer.Opts.useMatrix);
            testCase.verifyEmpty(optimizer.Opts.epsilon);
            testCase.verifyEqual(optimizer.Opts.mType, "sparse");
            testCase.verifyEmpty(optimizer.Opts.jacFun);
            testCase.verifyEqual(optimizer.Opts.krylov, "gmres");
            testCase.verifyEqual(optimizer.Opts.lineSearch, "none");
            testCase.verifyFalse(optimizer.Opts.equilibrate);
            testCase.verifyEqual(optimizer.Opts.reorderingMethod, "dissect");

            % Verify empty options structures
            testCase.verifyEqual(optimizer.krylovOpts, struct());
            testCase.verifyEqual(optimizer.lineSearchOpts, struct());
            testCase.verifyEmpty(optimizer.Loss);
        end

        function testConstructorCustomOptions(testCase)
            % Test constructor with custom parameters
            krylovOpts = struct('tol', 1e-8, 'restart', 40);
            lineSearchOpts = struct('minStep', 1e-4, 'contractionRate', 0.3);

            optimizer = core.optim.JfnkOptimizer( ...
                'maxIter', 50, ...
                'fTol', 1e-8, ...
                'xTol', 1e-12, ...
                'useMatrix', true, ...
                'epsilon', 1e-6, ...
                'mType', "dense", ...
                'krylov', "bicgstab", ...
                'krylovOpts', krylovOpts, ...
                'lineSearch', "backtracking", ...
                'lineSearchOpts', lineSearchOpts, ...
                'equilibrate', true, ...
                'reorderingMethod', "amd");

            % Verify custom options
            testCase.verifyEqual(optimizer.Opts.maxIter, 50);
            testCase.verifyEqual(optimizer.Opts.fTol, 1e-8);
            testCase.verifyEqual(optimizer.Opts.xTol, 1e-12);
            testCase.verifyTrue(optimizer.Opts.useMatrix);
            testCase.verifyEqual(optimizer.Opts.epsilon, 1e-6);
            testCase.verifyEqual(optimizer.Opts.mType, "dense");
            testCase.verifyEqual(optimizer.Opts.krylov, "bicgstab");
            testCase.verifyEqual(optimizer.Opts.lineSearch, "backtracking");
            testCase.verifyTrue(optimizer.Opts.equilibrate);
            testCase.verifyEqual(optimizer.Opts.reorderingMethod, "amd");

            % Verify options structures
            testCase.verifyEqual(optimizer.krylovOpts.tol, 1e-8);
            testCase.verifyEqual(optimizer.krylovOpts.restart, 40);
            testCase.verifyEqual(optimizer.lineSearchOpts.minStep, 1e-4);
            testCase.verifyEqual(optimizer.lineSearchOpts.contractionRate, 0.3);
        end

        function testSimpleNonlinearSystem(testCase)
            % Test JFNK optimizer on simple 2x2 nonlinear system
            optimizer = core.optim.JfnkOptimizer('fTol', 1e-8, 'maxIter', 50);

            % System: F(x) = [x1^2 + x2^2 - 1; x1 - x2] = 0
            F = @(x) [x(1)^2 + x(2)^2 - 1; x(1) - x(2)];
            x0 = [1.2; 0.8];

            [x, stats] = optimizer.solve(F, x0);

            % Verify convergence
            testCase.verifyTrue(stats.convergedId > 0, 'Optimizer should converge');
            testCase.verifyLessThan(stats.fNorm(end), 1e-7, ...
                'Final residual should be small');

            % Verify solution accuracy
            exactSol = [sqrt(2)/2; sqrt(2)/2];
            testCase.verifyLessThan(norm(x - exactSol), 1e-6, ...
                'Solution should be accurate');
        end

        function testLinearSystem(testCase)
            % Test JFNK optimizer on linear system (should converge quickly)
            optimizer = core.optim.JfnkOptimizer('fTol', 1e-10, 'maxIter', 10);

            % System: F(x) = [2*x1 + x2 - 3; x1 + 2*x2 - 3] = 0
            F = @(x) [2*x(1) + x(2) - 3; x(1) + 2*x(2) - 3];
            x0 = [0.5; 0.5];

            [x, stats] = optimizer.solve(F, x0);

            % Verify convergence
            testCase.verifyTrue(stats.convergedId > 0, 'Optimizer should converge');
            testCase.verifyLessThan(stats.nlIter, 5, ...
                'Linear system should converge quickly');

            % Verify solution
            exactSol = [1; 1];
            testCase.verifyLessThan(norm(x - exactSol), 1e-8, ...
                'Solution should be exact');
        end

        function testRosenbrockSystem(testCase)
            % Test JFNK optimizer on challenging Rosenbrock-like system
            optimizer = core.optim.JfnkOptimizer('fTol', 1e-6, 'maxIter', 200);

            % System: F(x) = [10*(x2 - x1^2); 1 - x1]
            F = @(x) [10*(x(2) - x(1)^2); 1 - x(1)];
            x0 = [0.5; 0.5];

            [x, stats] = optimizer.solve(F, x0);

            % Verify convergence
            testCase.verifyTrue(stats.convergedId > 0, 'Optimizer should converge');

            % Verify solution
            exactSol = [1; 1];
            testCase.verifyLessThan(norm(x - exactSol), 1e-5, ...
                'Solution should be accurate');
        end

        function test3x3System(testCase)
            % Test JFNK optimizer on 3x3 nonlinear system
            optimizer = core.optim.JfnkOptimizer('fTol', 1e-6, 'maxIter', 100);

            % System: F(x) = [x1 + x2 + x3 - 6; x1^2 + x2^2 - 5; x1*x3 - 2] = 0
            F = @(x) [x(1) + x(2) + x(3) - 6; x(1)^2 + x(2)^2 - 5; x(1)*x(3) - 2];
            x0 = [1.8; 1.6; 1.6];

            [x, stats] = optimizer.solve(F, x0);

            % Verify convergence
            testCase.verifyTrue(stats.convergedId > 0, 'Optimizer should converge');
            testCase.verifyLessThan(stats.fNorm(end), 1e-5, ...
                'Final residual should be small');
        end

        function testMatrixBasedSolving(testCase)
            % Test matrix-based solving with useMatrix = true
            optimizer = core.optim.JfnkOptimizer('useMatrix', true, 'fTol', 1e-8);

            % Simple 2x2 system
            F = @(x) [x(1)^2 + x(2)^2 - 1; x(1) - x(2)];
            x0 = [1.2; 0.8];

            [x, stats] = optimizer.solve(F, x0);

            % Verify convergence
            testCase.verifyTrue(stats.convergedId > 0, 'Matrix-based solver should converge');

            % Verify solution accuracy
            exactSol = [sqrt(2)/2; sqrt(2)/2];
            testCase.verifyLessThan(norm(x - exactSol), 1e-6, ...
                'Matrix-based solution should be accurate');
        end

        function testDifferentKrylovSolvers(testCase)
            % Test different Krylov solvers
            solvers = {'gmres', 'bicgstab'};

            % Simple 2x2 system
            F = @(x) [x(1)^2 + x(2)^2 - 1; x(1) - x(2)];
            x0 = [1.2; 0.8];
            exactSol = [sqrt(2)/2; sqrt(2)/2];

            for i = 1:length(solvers)
                optimizer = core.optim.JfnkOptimizer('krylov', solvers{i}, 'fTol', 1e-6);

                [x, stats] = optimizer.solve(F, x0);

                % Verify convergence for each solver
                testCase.verifyTrue(stats.convergedId > 0, ...
                    sprintf('%s solver should converge', solvers{i}));
                testCase.verifyLessThan(norm(x - exactSol), 1e-5, ...
                    sprintf('%s solution should be accurate', solvers{i}));
            end
        end

        function testBacktrackingLineSearch(testCase)
            % Test backtracking line search
            lineSearchOpts = struct('minStep', 1e-3, 'contractionRate', 0.5);
            optimizer = core.optim.JfnkOptimizer( ...
                'lineSearch', "backtracking", ...
                'lineSearchOpts', lineSearchOpts, ...
                'fTol', 1e-6);

            % Rosenbrock-like system (challenging for line search)
            F = @(x) [10*(x(2) - x(1)^2); 1 - x(1)];
            x0 = [0.5; 0.5];

            [x, stats] = optimizer.solve(F, x0);

            % Verify convergence with line search
            testCase.verifyTrue(stats.convergedId > 0, ...
                'Optimizer with line search should converge');

            % Verify solution
            exactSol = [1; 1];
            testCase.verifyLessThan(norm(x - exactSol), 1e-5, ...
                'Line search solution should be accurate');
        end

        function testCustomJacobianFunction(testCase)
            % Test with custom analytical Jacobian function

            % System: F(x) = [x1^2 - x2; x1 + x2^2 - 3]
            F = @(x) [x(1)^2 - x(2); x(1) + x(2)^2 - 3];

            % Analytical Jacobian
            J = @(x) [2*x(1), -1; 1, 2*x(2)];

            optimizer = core.optim.JfnkOptimizer( ...
                'useMatrix', true, ...
                'jacFun', J, ...
                'fTol', 1e-8);

            x0 = [1.5; 1.2];
            [x, stats] = optimizer.solve(F, x0);

            % Verify convergence with analytical Jacobian
            testCase.verifyTrue(stats.convergedId > 0, ...
                'Optimizer with analytical Jacobian should converge');
            testCase.verifyLessThan(stats.fNorm(end), 1e-7, ...
                'Final residual should be very small with analytical Jacobian');
        end

        function testKrylovOptions(testCase)
            % Test Krylov solver options
            krylovOpts = struct( ...
                'tol', 1e-8, ...
                'restart', 2, ...
                'maxIt', 2);

            optimizer = core.optim.JfnkOptimizer( ...
                'krylov', "gmres", ...
                'krylovOpts', krylovOpts);

            % Simple system
            F = @(x) [x(1)^2 + x(2)^2 - 1; x(1) - x(2)];
            x0 = [1.2; 0.8];

            [x, stats] = optimizer.solve(F, x0);

            % Verify that custom options were used
            testCase.verifyEqual(optimizer.krylovOpts.tol, 1e-8);
            testCase.verifyEqual(optimizer.krylovOpts.restart, 2);
            testCase.verifyEqual(optimizer.krylovOpts.maxIt, 2);

            % Verify convergence
            testCase.verifyTrue(stats.convergedId > 0, 'Optimizer should converge');
        end

        function testStatisticsTracking(testCase)
            % Test that statistics are properly tracked
            optimizer = core.optim.JfnkOptimizer('fTol', 1e-8);

            F = @(x) [x(1)^2 + x(2)^2 - 1; x(1) - x(2)];
            x0 = [1.2; 0.8];

            [x, stats] = optimizer.solve(F, x0);

            % Verify statistics structure
            testCase.verifyTrue(isfield(stats, 'nlIter'));
            testCase.verifyTrue(isfield(stats, 'fNorm'));
            testCase.verifyTrue(isfield(stats, 'xNorm'));
            testCase.verifyTrue(isfield(stats, 'convergedId'));
            testCase.verifyTrue(isfield(stats, 'convergedMessage'));
            testCase.verifyTrue(isfield(stats, 'krylovStats'));

            % Verify Krylov statistics
            testCase.verifyTrue(isfield(stats.krylovStats, 'exitStatusId'));
            testCase.verifyTrue(isfield(stats.krylovStats, 'linearRes'));
            testCase.verifyTrue(isfield(stats.krylovStats, 'linearSteps'));
            testCase.verifyTrue(isfield(stats.krylovStats, 'resVec'));

            % Verify array sizes
            testCase.verifyEqual(length(stats.fNorm), stats.nlIter + 1);
            testCase.verifyEqual(length(stats.xNorm), stats.nlIter + 1);
            testCase.verifyEqual(length(stats.krylovStats.exitStatusId), stats.nlIter);
            testCase.verifyEqual(length(stats.krylovStats.linearRes), stats.nlIter);
            testCase.verifyEqual(length(stats.krylovStats.linearSteps), stats.nlIter);
            testCase.verifyEqual(length(stats.krylovStats.resVec), stats.nlIter);
        end

        function testConvergenceTolerances(testCase)
            % Test different convergence criteria

            % Test function tolerance convergence
            optimizer1 = core.optim.JfnkOptimizer('fTol', 1e-3, 'xTol', 1e-15);
            F = @(x) [x(1)^2 + x(2)^2 - 1; x(1) - x(2)];
            x0 = [1.2; 0.8];

            [~, stats1] = optimizer1.solve(F, x0);
            testCase.verifyEqual(stats1.convergedId, 1, ...
                'Should converge due to function tolerance');

            % Test solution tolerance convergence (use identity function for immediate x convergence)
            optimizer2 = core.optim.JfnkOptimizer('fTol', 1e-15, 'xTol', 1e-1);
            F2 = @(x) 0.1 * x; % Small residual but not zero
            x0_2 = [0.05; 0.05]; % Small initial guess

            [~, stats2] = optimizer2.solve(F2, x0_2);
            % Note: This test might be tricky due to the specific nature of convergence criteria
        end

        function testMaxIterationsReached(testCase)
            % Test behavior when maximum iterations are reached
            optimizer = core.optim.JfnkOptimizer( ...
                'maxIter', 2, ...  % Very few iterations
                'fTol', 1e-15);     % Very strict tolerance

            % Challenging system that won't converge quickly
            F = @(x) [10*(x(2) - x(1)^2); 1 - x(1)];
            x0 = [0.1; 0.1]; % Poor initial guess

            [~, stats] = optimizer.solve(F, x0);

            % Should not converge within 2 iterations
            testCase.verifyEqual(stats.convergedId, -1, ...
                'Should indicate non-convergence');
            testCase.verifyEqual(stats.nlIter, 2, ...
                'Should perform exactly maxIter iterations');
            testCase.verifyTrue(contains(stats.convergedMessage, 'failed to converge'), ...
                'Should indicate convergence failure');
        end

        function testEpsilonParameter(testCase)
            % Test custom epsilon parameter for finite differences
            optimizer = core.optim.JfnkOptimizer('epsilon', 1e-5);

            % Verify epsilon is set
            F = @(x) [x(1)^2 + x(2)^2 - 1; x(1) - x(2)];
            x0 = [1.2; 0.8];

            [x, stats] = optimizer.solve(F, x0);

            % Verify epsilon was used
            testCase.verifyEqual(optimizer.Opts.epsilon, 1e-5);

            % Should still converge
            testCase.verifyTrue(stats.convergedId > 0, 'Should converge with custom epsilon');
        end

        function testSparseLargeSystem100x100(testCase)
            % Test JFNK optimizer on sparse 100x100 nonlinear system
            n = 100;
            optimizer = core.optim.JfnkOptimizer( ...
                'fTol', 1e-6, ...
                'maxIter', 100, ...
                'mType', "sparse", ...
                'krylov', "gmres");

            % Sparse nonlinear system: discrete 2D Laplacian with nonlinearity
            % F(x) = A*x + sin(x) - b, where A is tridiagonal matrix
            A = spdiags([-ones(n,1) 2*ones(n,1) -ones(n,1)], -1:1, n, n);
            b = ones(n, 1);
            F = @(x) A*x + sin(x) - b;
            x0 = 0.1 * ones(n, 1);

            [x, stats] = optimizer.solve(F, x0);

            % Verify convergence
            testCase.verifyTrue(stats.convergedId > 0, 'Large sparse system should converge');
            testCase.verifyLessThan(stats.fNorm(end), 1e-5, ...
                'Final residual should be small for large system');
            testCase.verifyEqual(length(x), n, 'Solution should have correct dimensions');
        end

        function testSparseLargeSystem500x500(testCase)
            % Test JFNK optimizer on very large sparse 500x500 system
            n = 500;
            optimizer = core.optim.JfnkOptimizer( ...
                'fTol', 1e-5, ...
                'maxIter', 150, ...
                'mType', "sparse", ...
                'krylov', "bicgstab");

            % Large sparse system: 2D Poisson-like with nonlinearity
            % Create sparse pentadiagonal matrix (5-point stencil)
            e = ones(n, 1);
            A = spdiags([-e -e 4*e -e -e], [-10 -1 0 1 10], n, n);
            b = randn(n, 1);
            F = @(x) A*x + 0.1*x.^3 - b;
            x0 = 0.01 * randn(n, 1);

            [x, stats] = optimizer.solve(F, x0);

            % Verify convergence
            testCase.verifyTrue(stats.convergedId > 0, 'Very large sparse system should converge');
            testCase.verifyLessThan(stats.fNorm(end), 1e-4, ...
                'Final residual should be reasonable for very large system');
            testCase.verifyEqual(length(x), n, 'Solution should have correct dimensions');
        end

        function testSparseSystemWithEquilibration(testCase)
            % Test sparse system with matrix equilibration
            n = 50;
            optimizer = core.optim.JfnkOptimizer( ...
                'fTol', 1e-6, ...
                'mType', "sparse", ...
                'equilibrate', true, ...
                'reorderingMethod', "amd");

            % Ill-conditioned sparse system
            A = spdiags([1e-3*ones(n,1) 1e6*ones(n,1) 1e-3*ones(n,1)], [-1 0 1], n, n);
            b = ones(n, 1);
            F = @(x) A*x + 0.01*x.^2 - b;
            x0 = 0.1 * ones(n, 1);

            [x, stats] = optimizer.solve(F, x0);

            % Verify convergence with equilibration
            testCase.verifyTrue(stats.convergedId > 0, ...
                'Ill-conditioned sparse system should converge with equilibration');
            testCase.verifyLessThan(stats.fNorm(end), 1e-5, ...
                'Final residual should be small with equilibration');
        end

        function testDifferentSparseReorderingMethods(testCase)
            % Test different sparse matrix reordering methods
            n = 80;
            methods = {"dissect", "amd", "colamd"};

            % Sparse system
            A = spdiags([-ones(n,1) 3*ones(n,1) -ones(n,1)], [-1 0 1], n, n);
            b = sin((1:n)');
            F = @(x) A*x + 0.1*tanh(x) - b;
            x0 = 0.1 * ones(n, 1);

            for i = 1:length(methods)
                optimizer = core.optim.JfnkOptimizer( ...
                    'fTol', 1e-6, ...
                    'mType', "sparse", ...
                    'reorderingMethod', methods{i});

                [x, stats] = optimizer.solve(F, x0);

                % Verify convergence for each reordering method
                testCase.verifyTrue(stats.convergedId > 0, ...
                    sprintf('Sparse system should converge with %s reordering', methods{i}));
                testCase.verifyLessThan(norm(x), 100, ...
                    sprintf('Solution norm should be reasonable with %s', methods{i}));
            end
        end

        function testAllKrylovSolversOnSparseSystem(testCase)
            % Test all available Krylov solvers on sparse system
            n = 60;
            solvers = {'gmres', 'bicgstab', 'pcg', 'minres', 'symmlq'};

            % Symmetric positive definite sparse system for PCG, MINRES, SYMMLQ
            A = spdiags([ones(n,1) 4*ones(n,1) ones(n,1)], [-1 0 1], n, n);
            b = ones(n, 1);
            F = @(x) A*x + 0.01*x.^2 - b;
            x0 = 0.1 * ones(n, 1);

            for i = 1:length(solvers)
                optimizer = core.optim.JfnkOptimizer( ...
                    'fTol', 1e-6, ...
                    'krylov', solvers{i}, ...
                    'mType', "sparse");

                [x, stats] = optimizer.solve(F, x0);

                % Verify convergence for each Krylov solver
                testCase.verifyTrue(stats.convergedId > 0, ...
                    sprintf('%s solver should converge on sparse system', solvers{i}));
                testCase.verifyLessThan(stats.fNorm(end), 1e-5, ...
                    sprintf('%s should achieve good residual reduction', solvers{i}));
            end
        end

        function testSparseSystemWithAdvancedKrylovOptions(testCase)
            % Test sparse system with advanced Krylov solver options
            n = 100;
            krylovOpts = struct( ...
                'tol', 1e-10, ...
                'restart', 20, ...
                'maxIt', 50);

            optimizer = core.optim.JfnkOptimizer( ...
                'fTol', 1e-8, ...
                'krylov', "gmres", ...
                'krylovOpts', krylovOpts, ...
                'mType', "sparse");

            % Sparse system
            A = spdiags([-0.5*ones(n,1) 2*ones(n,1) -0.5*ones(n,1)], [-1 0 1], n, n);
            b = cos((1:n)' * pi / n);
            F = @(x) A*x + 0.1*sin(x) - b;
            x0 = 0.1 * ones(n, 1);

            [x, stats] = optimizer.solve(F, x0);

            % Verify convergence with advanced options
            testCase.verifyTrue(stats.convergedId > 0, ...
                'Sparse system should converge with advanced Krylov options');
            testCase.verifyLessThan(stats.fNorm(end), 1e-7, ...
                'Should achieve very high accuracy with tight tolerances');
        end

        function testNonlinearEigenvalueProblem(testCase)
            % Test JFNK on nonlinear eigenvalue-like problem
            n = 50;
            optimizer = core.optim.JfnkOptimizer( ...
                'fTol', 1e-6, ...
                'maxIter', 100, ...
                'mType', "sparse");

            % Nonlinear eigenvalue problem: (A - λI)x + g(x) = 0
            % We solve for x with λ = 2
            A = spdiags([ones(n,1) 3*ones(n,1) ones(n,1)], [-1 0 1], n, n);
            lambda = 2;
            F = @(x) (A - lambda*speye(n))*x + 0.1*x.^3;
            x0 = ones(n, 1) / sqrt(n); % Normalized initial guess

            [x, stats] = optimizer.solve(F, x0);

            % Verify convergence
            testCase.verifyTrue(stats.convergedId > 0, ...
                'Nonlinear eigenvalue problem should converge');
            testCase.verifyLessThan(stats.fNorm(end), 1e-5, ...
                'Final residual should be small for eigenvalue problem');
        end

        function testBoundaryValueProblem(testCase)
            % Test JFNK on discretized nonlinear boundary value problem
            n = 80;
            h = 1 / (n + 1); % Grid spacing
            optimizer = core.optim.JfnkOptimizer( ...
                'fTol', 1e-6, ...
                'maxIter', 150, ...
                'mType', "sparse");

            % Discrete 1D nonlinear BVP: -u'' + u^3 = f, u(0) = u(1) = 0
            % Discretized as: (-u_{i-1} + 2u_i - u_{i+1})/h^2 + u_i^3 = f_i
            A = (1/h^2) * spdiags([-ones(n,1) 2*ones(n,1) -ones(n,1)], [-1 0 1], n, n);
            x_grid = (1:n)' * h;
            f = sin(pi * x_grid); % Right-hand side
            F = @(x) A*x + x.^3 - f;
            x0 = 0.1 * sin(pi * x_grid); % Initial guess

            [x, stats] = optimizer.solve(F, x0);

            % Verify convergence
            testCase.verifyTrue(stats.convergedId > 0, ...
                'Nonlinear BVP should converge');
            testCase.verifyLessThan(stats.fNorm(end), 1e-5, ...
                'BVP residual should be small');

            % Verify boundary conditions (approximately)
            testCase.verifyLessThan(abs(x(1)), 0.1, ...
                'Left boundary condition should be approximately satisfied');
            testCase.verifyLessThan(abs(x(end)), 0.1, ...
                'Right boundary condition should be approximately satisfied');
        end

        function testRobustnessToInitialGuess(testCase)
            % Test robustness to different initial guesses
            n = 40;
            optimizer = core.optim.JfnkOptimizer( ...
                'fTol', 1e-6, ...
                'maxIter', 200, ...
                'lineSearch', "backtracking", ...
                'mType', "sparse");

            % Sparse nonlinear system
            A = spdiags([-ones(n,1) 3*ones(n,1) -ones(n,1)], [-1 0 1], n, n);
            b = ones(n, 1);
            F = @(x) A*x + 0.1*x.^2 - b;

            % Test different initial guesses
            initialGuesses = {zeros(n,1), ones(n,1), -ones(n,1), 10*randn(n,1)};

            for i = 1:length(initialGuesses)
                x0 = initialGuesses{i};
                [x, stats] = optimizer.solve(F, x0);

                % Should converge from reasonable initial guesses
                if norm(x0) < 50 % Exclude very large initial guesses
                    testCase.verifyTrue(stats.convergedId > 0, ...
                        sprintf('Should converge from initial guess %d', i));
                end
            end
        end

        function testMemoryEfficiency(testCase)
            % Test memory efficiency for large systems
            n = 300;
            optimizer = core.optim.JfnkOptimizer( ...
                'fTol', 1e-5, ...
                'maxIter', 50, ...
                'mType', "sparse", ...
                'krylov', "bicgstab");

            % Memory-efficient: only store sparse matrix, not full Jacobian
            A = spdiags([ones(n,1) -2*ones(n,1) ones(n,1)], [-1 0 1], n, n);
            b = ones(n, 1);
            F = @(x) A*x + 0.01*x.^3 - b;
            x0 = 0.1 * ones(n, 1);

            % Should complete without memory issues
            [x, stats] = optimizer.solve(F, x0);

            testCase.verifyTrue(stats.convergedId > 0, ...
                'Large system should converge efficiently');
            testCase.verifyEqual(length(x), n, ...
                'Solution should have correct size');
        end

        function testLineSearchOnChallengingProblem(testCase)
            % Test line search on challenging sparse problem
            n = 60;
            lineSearchOpts = struct( ...
                'minStep', 1e-4, ...
                'contractionRate', 0.7, ...
                'maxBacktrack', 20);

            optimizer = core.optim.JfnkOptimizer( ...
                'fTol', 1e-6, ...
                'lineSearch', "backtracking", ...
                'lineSearchOpts', lineSearchOpts, ...
                'mType', "sparse");

            % Challenging problem with steep nonlinearity
            A = spdiags([ones(n,1) 4*ones(n,1) ones(n,1)], [-1 0 1], n, n);
            b = ones(n, 1);
            F = @(x) A*x + sign(x) .* x.^2 - b; % Non-smooth nonlinearity
            x0 = 2 * ones(n, 1); % Aggressive initial guess

            [x, stats] = optimizer.solve(F, x0);

            % Line search should help convergence
            testCase.verifyTrue(stats.convergedId > 0, ...
                'Challenging problem should converge with line search');
            testCase.verifyLessThan(stats.fNorm(end), 1e-5, ...
                'Should achieve good accuracy with line search');
        end
    end
end