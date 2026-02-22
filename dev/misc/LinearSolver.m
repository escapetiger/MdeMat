classdef LinearSolver < handle
    % LINEARSOLVER Adaptive linear solver.
    %
    %   LinearSolver automatically detects matrix properties like sparsity
    %   and positive-definiteness to choose between direct methods and
    %   iterative methods for optimal performance.

    properties
        Alpha % Size threshold for choosing direct vs iterative methods
        Beta % Condition number threshold for preconditioning
        Tol % Tolerance for iterative solvers
        MaxIter % Maximum iterations for iterative solvers
        IsSpd % Flag indicating if matrix is SPD (empty = auto-detect)
    end

    methods
        function obj = LinearSolver(options)
            % LINEARSOLVER Construct an instance of LinearSolver.
            %
            %   obj = LinearSolver() creates LinearSolver with default
            %   settings.
            %
            %   obj = LinearSolver(Name=Value) creates LinearSolver with
            %   specified options.

            arguments
                options.alpha{mustBeNumeric, mustBePositive} = inf
                options.beta{mustBeNumeric, mustBePositive} = 1e4
                options.tol{mustBeNumeric, mustBePositive} = 1e-8
                options.maxIter{mustBeInteger, mustBePositive} = 500
                options.isSpd{mustBeNumericOrLogical} = []
            end

            obj.Alpha = options.alpha;
            obj.Beta = options.beta;
            obj.Tol = options.tol;
            obj.MaxIter = options.maxIter;
            obj.IsSpd = options.isSpd;
        end

        function x = solve(obj, A, b)
            % SOLVE Solve linear system @a A*x = @a b.
            %
            %   x = solve(obj, A, b) automatically selects and applies
            %   appropriate solver for linear system @a A*x = @a b.
            %
            % Notes:
            %   Automatically detects matrix properties and selects the
            %   most appropriate solver method. Matrix @a A must be square.

            arguments
                obj core.linalg.LinearSolver
                A{mustBeNumeric, mustBeSquare(A)}
                b{mustBeNumeric}
            end

            nnzs = nnz(A);

            if nnzs > obj.Alpha
                x = obj.iterative(A, b);
            else
                x = obj.direct(A, b);
            end
        end
    end

    methods (Access = private)
        function tf = checkSpd(obj, A)
            % ISSPD Check if A is symmetric positive definite.
            tf = obj.IsSpd;
            if isempty(tf)
                if ~issymmetric(A)
                    tf = false;
                else
                    [~, p] = chol(A);
                    tf = (p == 0);
                end
            end
        end

        function x = iterative(obj, A, b)
            % ITERATIVE Call iterative methods.
            if obj.checkSpd(A)
                [x, f, rel, iter] = pcg(A, b, obj.Tol, obj.MaxIter);
                core.except.verify(f == 0, 'ConvergenceFailed', ...
                    ['PCG did not converge. ', ...
                    'Flag: %d, Relative: %g, Iterations: %d'], f, rel, iter);
            elseif issparse(A)
                n = min(obj.MaxIter, size(A, 1));
                restart = 50;
                [x, f, rel, iter] = gmres(A, b, restart, obj.Tol, n);

                core.except.verify(f == 0, 'ConvergenceFailed', ...
                    ['GMRES did not converge. ', ...
                    'Flag: %d, Relative: %g, Iterations: %d'], f, rel, iter);
            else
                [x, f, rel, iter] = bicgstab(A, b, obj.Tol, obj.MaxIter);

                core.except.verify(f == 0, 'ConvergenceFailed', ...
                    ['BICGSTAB did not converge. ', ...
                    'Flag: %d, Relative: %g, Iterations: %d'], f, rel, iter);
            end
        end

        function x = direct(obj, A, b)
            % DIRECT Call direct methods.
            if obj.checkSpd(A)
                try
                    R = chol(A);
                    x = R \ (R' \ b);
                catch ME
                    core.except.verify(0, 'CholFailed', ...
                        ['Cholesky factorization failed: %s. ', ...
                        'Falling back to backslash.'], ME.message);
                    x = A \ b;
                end
            else
                x = A \ b;
            end
        end
    end
end

function mustBeSquare(A)
if size(A, 1) ~= size(A, 2)
    error('Matrix must be square');
end
end