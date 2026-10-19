classdef LinearSolver < handle
    % LINEARSOLVER Adpative linear solver.
    %
    %   LinearSolver automatically detects matrix properties like sparsity
    %   and positive-definiteness to choose between direct methods and
    %   iterative methods for optimal performance.
    %
    % Examples:
    %   % Create solver and solve system
    %   solver = core.linalg.LinearSolver();
    %   A = sprand(1000, 1000, 0.01) + speye(1000);
    %   b = rand(1000, 1);
    %   x = solver.solve(A, b);

    properties (Constant)
        SIZE_THRESHOLD = 1e5 % Size threshold
    end

    properties
        tol % Tolerance for iterative solvers
        maxIter % Maximum iterations for iterative solvers
        spd % Flag indicating if matrix is SPD (empty = auto-detect)
    end

    methods
        function obj = LinearSolver()
            % LINEARSOLVER Constructor for LinearSolver.
            %
            %   obj = LinearSolver() creates empty LinearSolver.
            %
            % Outputs:
            %   obj - Constructed LinearSolver object

            obj.tol = 1e-6;
            obj.maxIter = 500;
            obj.spd = [];
        end

        function x = solve(obj, A, b)
            % SOLVE Solve linear system @a A*x = @a b.
            %
            %   x = solve(obj, A, b) automatically selects and applies
            %   appropriate solver for linear system @a A*x = @a b.
            %
            % Inputs:
            %   obj - The LinearSolver object
            %   A - Left-hand side square matrix
            %   b - Right-hand side vector
            %
            % Outputs:
            %   x - Solution vector to the system @a A*x = @a b

            nnzs = nnz(A);

            if nnzs > obj.SIZE_THRESHOLD
                x = obj.iterative(A, b);
            else
                x = obj.direct(A, b);
            end
        end
    end

    methods (Access = private)
        function f = isSpd(obj, A)
            % ISSPD Check if A is symmetric positive definite.
            f = obj.spd;
            if isempty(f)
                if ~issymmetric(A)
                    f = false;
                else
                    [~, p] = chol(A);
                    f = (p == 0);
                end
            end
        end

        function x = iterative(obj, A, b)
            % ITERATIVE Call itervative methods.
            if obj.isSpd(A)
                [x, f, rel, iter] = pcg(A, b, obj.tol, obj.maxIter);
                core.except.verify(f == 0, 'ConvergenceFailed', ...
                    ['PCG did not converge. ' ...
                    'Flag: %d, Relative: %g, Iterations: %d'], f, rel, iter);
            elseif issparse(A)
                [x, f, rel, iter] = gmres(A, b, [], obj.tol, obj.maxIter);
                core.except.verify(f == 0, 'ConvergenceFailed', ...
                    ['GMRES did not converge. ' ...
                    'Flag: %d, Relative: %g, Iterations: %d'], f, rel, iter);
            else
                [x, f, rel, iter] = bicgstab(A, b, obj.tol, obj.maxIter);
                core.except.verify(f == 0, 'ConvergenceFailed', ...
                    ['BICGSTAB did not converge. ' ...
                    'Flag: %d, Relative: %g, Iterations: %d'], f, rel, iter);
            end
        end

        function x = direct(obj, A, b)
            % DIRECT Call direct methods.
            if obj.isSpd(A)
                try
                    R = chol(A);
                    x = R \ (R' \ b);
                catch ME
                    core.except.verify(0, 'CholFailed', ...
                        ['Cholesky factorization failed: %s. ' ...
                        'Falling back to backslash.'], ME.message);
                    x = A \ b;
                end
            else
                x = A \ b;
            end
        end
    end
end