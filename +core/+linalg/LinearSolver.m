classdef LinearSolver < handle
    % LINEARSOLVER Adaptive linear solver with stateful caching.
    %
    %   LinearSolver automatically detects matrix properties like sparsity
    %   and condition number to choose between  methods and iterative
    %   methods with optional preconditioning for optimal performance.
    %   Uses MATLAB's optimized decomposition function for  factorization.
    %
    %   The solver uses sparsity and condition number thresholds to
    %   determine the optimal solution strategy, with factorization caching
    %   for repeated solves with the same matrix.

    properties
        NnzTh % Sparsity threshold for direct vs iterative method selection
        CondTh % Condition number threshold for preconditioner activation
        KrylovOpts % Options structure for iterative solver
        Lhs % Cached left-hand side matrix
        Rhs % Cached right-hand side vector
        Decomp % Cached matrix decomposition with metadata
        IsPrecomputed % Flag indicating if precomputation is current
        LhsHash % Hash of current LHS matrix for change detection
    end

    methods
        function obj = LinearSolver(options)
            % LINEARSOLVER Construct an instance of LinearSolver.
            %
            %   obj = LinearSolver() creates LinearSolver with default
            %   settings.
            %
            %   obj = LinearSolver(Name=Value) creates LinearSolver with
            %   specified options using @a Name and @a Value pairs.

            arguments
                options.nnzTh{mustBeNumeric, mustBePositive} = 1e5
                options.condTh{mustBeNumeric, mustBePositive} = inf
                options.maxIter{mustBeInteger, mustBePositive} = 1000
                options.tol{mustBeNumeric, mustBePositive} = 1e-8
                options.restart{mustBeInteger, mustBePositive} = 50
            end

            obj.NnzTh = options.nnzTh;
            obj.CondTh = options.condTh;

            %< Initialize options structure
            obj.KrylovOpts = struct();
            obj.KrylovOpts.maxIter = options.maxIter;
            obj.KrylovOpts.tol = options.tol;
            obj.KrylovOpts.restart = options.restart;

            %< Initialize stateful properties
            obj.Lhs = [];
            obj.Rhs = [];
            obj.Decomp = struct('type', 'none', 'data', [], 'isSpd', false);
            obj.IsPrecomputed = false;
            obj.LhsHash = [];
        end

        function obj = setLhs(obj, A)
            % SETLhs Set the left-hand side matrix.
            %
            %   obj = setLhs(obj, A) stores the left-hand side matrix @a A
            %   for repeated solves. Uses hash-based change detection to
            %   avoid unnecessary refactorization when the same matrix is
            %   set multiple times.

            arguments
                obj core.linalg.LinearSolver
                A{mustBeNumeric, mustBeSquare(A)}
            end

            %< Compute hash of new matrix
            newHash = core.linalg.hash(A);
            % newHash = '';
            
            %< Check if matrix has actually changed
            if ~isempty(obj.LhsHash) && strcmp(obj.LhsHash, newHash)
                return; % Same matrix, no need to update
            end

            %< Update matrix and invalidate cached factorization
            obj.Lhs = A;
            obj.LhsHash = newHash;
            obj.IsPrecomputed = false;
            obj.Decomp = struct('type', 'none', 'data', [], 'isSpd', false);
        end

        function obj = setRhs(obj, b)
            % SETRhs Set the right-hand side vector.
            %
            %   obj = setRhs(obj, b) stores the right-hand side vector @a b
            %   for solving.

            arguments
                obj core.linalg.LinearSolver
                b{mustBeNumeric}
            end

            obj.Rhs = b;
        end

        function obj = precompute(obj)
            % PRECOMPUTE Factorize the left-hand side matrix.
            %
            %   obj = precompute(obj) computes and caches the matrix
            %   factorization for the stored left-hand side matrix.
            %   Subsequent solves will use the cached factorization.

            arguments
                obj core.linalg.LinearSolver
            end

            core.except.verify(~isempty(obj.Lhs), 'InvalidState', ...
                'Left-hand side matrix must be set before precomputing.');

            if obj.IsPrecomputed
                return;
            end

            %< Compute factorization based on matrix properties
            if nnz(obj.Lhs) <= obj.NnzTh
                %< Use builtin factorization for small/dense matrices
                obj.Decomp.data = decomposition(obj.Lhs);
                obj.Decomp.type = 'auto';
                obj.Decomp.isSpd = false;
            else
                %< Fast SPD check for large sparse matrices
                obj.Decomp.isSpd = obj.isSpd(obj.Lhs);

                %< Setup preconditioner for large sparse matrices
                if issparse(obj.Lhs) && obj.condest(obj.Lhs) > obj.CondTh
                    p = amd(obj.Lhs);
                    Ap = obj.Lhs(p, p);
                    opts = struct();
                    opts.type = 'ilutp';
                    opts.droptol = 1e-6;
                    opts.milu = 'row';
                    opts.udiag = 0;
                    opts.thresh = 0.1;
                    [L, U] = ilu(Ap, opts);
                    obj.Decomp.data.L = L;
                    obj.Decomp.data.U = U;
                    obj.Decomp.data.p = p;
                    obj.Decomp.type = 'ilu';
                else
                    obj.Decomp.type = 'none';
                end
            end

            obj.IsPrecomputed = true;
        end

        function [x, stats] = solve(obj, options)
            % SOLVE Solve linear system using adaptive method selection.
            %
            %   [x, stats] = solve(obj) solves the linear system using the cached
            %   left-hand side matrix and right-hand side vector. Uses
            %   cached factorization if available. Returns solution @a x and
            %   statistics structure @a stats.
            %
            %   [x, stats] = solve(obj, b=b) solves using cached left-hand side
            %   matrix and the provided right-hand side vector @a b.
            %
            %   [x, stats] = solve(obj, A=A, b=b) solves the linear system @a A * x
            %   = @a b using automatic method selection.

            arguments
                obj core.linalg.LinearSolver
                options.A{mustBeNumeric, mustBeSquare} = []
                options.b{mustBeNumeric} = []
            end

            if ~isempty(options.A)
                obj.setLhs(options.A);
            end

            if ~isempty(options.b)
                obj.setRhs(options.b);
            end

            if ~obj.IsPrecomputed
                obj.precompute();
            end

            % Initialize statistics
            stats = struct();
            stats.matrixSize = size(obj.Lhs);
            stats.matrixNnz = nnz(obj.Lhs);
            stats.converged = true;

            if nnz(obj.Lhs) <= obj.NnzTh
                stats.method = 'direct';
                [x, stats] = obj.direct(obj.Rhs, stats);
            else
                stats.method = 'iterative';
                [x, stats] = obj.iterative(obj.Rhs, stats);
            end
        end
    end

    methods (Access = private)
        function [x, stats] = direct(obj, b, stats)
            % DIRECT Solve using cached factorization.
            stats.decompositionType = obj.Decomp.type;
            stats.iterations = 0;
            stats.relativeResidual = 0;
            stats.flag = 0;

            switch obj.Decomp.type
                case 'auto'
                    x = obj.Decomp.data \ b;
                    stats.solverType = 'decomposition';
                otherwise
                    x = obj.Lhs \ b;
                    stats.solverType = 'backslash';
            end
        end

        function [x, stats] = iterative(obj, b, stats)
            % ITERATIVE Call iterative methods using cached matrix.
            A = obj.Lhs;
            stats.isSpd = obj.Decomp.isSpd;
            stats.preconditioner = obj.Decomp.type;

            if obj.Decomp.isSpd
                maxIter = obj.KrylovOpts.maxIter;
                tol = obj.KrylovOpts.tol;
                [x, f, rel, iter] = pcg(A, b, tol, maxIter);

                stats.solverType = 'pcg';
                stats.flag = f;
                stats.relativeResidual = rel;
                stats.iterations = iter;
                stats.converged = (f == 0);

                core.except.verify(f == 0, 'ConvergenceFailed', ...
                    ['PCG did not converge. ', ...
                    'Flag: %d, Relative: %g, Iterations: %d'], f, rel, iter);
                return;
            end

            maxIter = min(obj.KrylovOpts.maxIter, size(A, 1));
            restart = obj.KrylovOpts.restart;
            tol = obj.KrylovOpts.tol;

            switch obj.Decomp.type
                case 'ilu'
                    L = obj.Decomp.data.L;
                    U = obj.Decomp.data.U;
                    p = obj.Decomp.data.p;
                    Ap = A(p, p);
                    bp = b(p);
                    [xp, f, rel, iter] = gmres(Ap, bp, restart, tol, ...
                        maxIter, L, U);
                    x(p) = xp;
                    x = x.';
                    stats.solverType = 'gmres_ilu';
                    stats.reordering = true;
                otherwise
                    [x, f, rel, iter] = gmres(A, b, restart, tol, maxIter);
                    stats.solverType = 'gmres';
                    stats.reordering = false;
            end

            stats.flag = f;
            stats.relativeResidual = rel;
            stats.iterations = iter;
            stats.converged = (f == 0);
            stats.restart = restart;
            stats.tolerance = tol;
            stats.maxIterations = maxIter;

            core.except.verify(f == 0, 'ConvergenceFailed', ...
                ['GMRES did not converge. ', ...
                'Flag: %d, Relative: %g, Iterations: %d'], f, rel, iter);
        end

    end

    methods (Static, Access = private)
        function tf = isDiagonal(A)
            % ISDIAGONAL Check if matrix A is diagonal.
            if issparse(A)
                % Sparse case: only check positions of nonzeros
                [i, j, ~] = find(A);
                tf = all(i == j);
            else
                % Dense case: check off-diagonal elements directly
                tf = all(A(~eye(size(A, 1))) == 0);
            end
        end

        function tf = isSpd(A)
            % ISSPD Fast heuristic SPD check for large sparse matrices.
            if isempty(A) || ~issparse(A) || ~issymmetric(A)
                tf = false;
                return;
            end

            %< Check diagonal positivity (necessary condition)
            d = diag(A);
            if any(d <= 0)
                tf = false;
                return;
            end

            %< Check diagonal dominance (sufficient condition for SPD)
            offDiagSums = sum(abs(A), 2) - abs(d);
            if all(d > offDiagSums)
                tf = true;
                return;
            end

            %< Gershgorin circle test for positive eigenvalues
            minEig = min(d - offDiagSums);
            if minEig <= 0
                tf = false;
                return;
            end

            %< Heuristic: if matrix has reasonable condition and structure, assume SPD
            %< This trades some accuracy for speed on large matrices
            tf = true;
        end

        function kappa = condest(A)
            % CONDEST Fast condition number estimate for large sparse
            % matrices Uses power iteration and inverse power iteration for
            % largest/smallest eigenvalues

            n = size(A, 1);
            maxIter = 20; % Keep iterations low for speed
            tol = 1e-3; % Relaxed tolerance for estimate

            %< Estimate largest eigenvalue using power iteration
            v = randn(n, 1);
            v = v / norm(v);

            for iter = 1:maxIter
                vOld = v;
                v = A * v;
                maxLambda = norm(v);
                v = v / maxLambda;

                if norm(v-vOld) < tol
                    break;
                end
            end

            %< Estimate smallest eigenvalue using diagonal dominance
            %< approximation This is much faster than inverse power
            %< iteration
            d = abs(diag(A));
            offDiagSums = sum(abs(A), 2) - d;

            %< Gershgorin circle estimate for smallest eigenvalue
            minLambdaEst = min(d-offDiagSums);

            %< If matrix might be nearly singular, use a conservative
            %< estimate
            if minLambdaEst <= 0
                % Use Frobenius norm scaling approach
                minLambdaEst = norm(A, 'fro') / n^2;
            end

            %< Condition number estimate
            cond_est = abs(maxLambda/minLambdaEst);

            %< Cap at reasonable value to avoid infinity
            kappa = min(cond_est, 1e16);
        end
    end
end

function mustBeSquare(A)
if size(A, 1) ~= size(A, 2)
    error('Matrix must be square');
end
end
