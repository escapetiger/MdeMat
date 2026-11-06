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
        LhsHash % Optional diagnostic identity; equality controls reuse
    end

    properties (SetAccess = private)
        RecordHash (1, 1) logical % Explicitly enable diagnostic array hashing
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
                options.recordHash (1, 1) logical = false
            end

            obj.NnzTh = options.nnzTh;
            obj.CondTh = options.condTh;
            obj.RecordHash = options.recordHash;

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
            %   for repeated solves. Exact equality controls reuse. Hashing
            %   is an optional diagnostic and is disabled by default.

            arguments
                obj mdemat.common.linalg.LinearSolver
                A{mustBeNumeric, mustBeSquare(A)}
            end

            if ~isempty(obj.Lhs) && isequaln(obj.Lhs, A)
                return;
            end
            newHash = [];
            if obj.RecordHash, newHash = mdemat.tools.base.hash(A); end
            
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
                obj mdemat.common.linalg.LinearSolver
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
                obj mdemat.common.linalg.LinearSolver
            end

            assert(all(~isempty(obj.Lhs), 'all'), 'mdemat:common:linalg:LinearSolver:InvalidState', ...
                '%s', 'Left-hand side matrix must be set before precomputing.');

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
                if issparse(obj.Lhs) && isfinite(obj.CondTh) && ...
                        obj.estimateCondition(obj.Lhs) > obj.CondTh
                    if obj.Decomp.isSpd
                        obj.Decomp.data.L = ichol(obj.Lhs, struct('type', 'ict', 'droptol', 1e-6));
                        obj.Decomp.type = 'ichol';
                        obj.IsPrecomputed = true;
                        return;
                    end
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
                obj mdemat.common.linalg.LinearSolver
                options.A{mustBeNumeric, mustBeSquare} = []
                options.b{mustBeNumeric} = []
            end

            if ~isempty(options.A)
                obj.setLhs(options.A);
            end

            if ~isempty(options.b)
                obj.setRhs(options.b);
            end

            if size(obj.Rhs, 1) ~= size(obj.Lhs, 1) || isempty(obj.Rhs)
                error('mdemat:linalg:InvalidRhs', ...
                    'Right-hand side must have one row per matrix row.');
            end
            if any(~isfinite(obj.Rhs(:))) || any(~isfinite(nonzeros(obj.Lhs)))
                error('mdemat:linalg:NonfiniteSystem', 'System entries must be finite.');
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
            residual = norm(obj.Lhs * x - b, 'fro');
            denominator = norm(b, 'fro');
            if denominator > 0
                residual = residual / denominator;
            end
            stats.relativeResidual = residual;
            stats.converged = all(isfinite(x(:))) && residual <= max(obj.KrylovOpts.tol, 100 * eps);
            stats.flag = double(~stats.converged);
            if ~stats.converged
                error('mdemat:linalg:ConvergenceFailed', ...
                    'Direct solve failed its residual check (%g).', residual);
            end
        end

        function [x, stats] = iterative(obj, b, stats)
            % ITERATIVE Call iterative methods using cached matrix.
            if size(b, 2) > 1
                x = zeros(size(b), 'like', b);
                reports = cell(1, size(b, 2));
                for column = 1:size(b, 2)
                    [x(:, column), reports{column}] = obj.iterative(b(:, column), stats);
                end
                stats.rhsReports = [reports{:}];
                stats.converged = all([stats.rhsReports.converged]);
                stats.flag = max([stats.rhsReports.flag]);
                stats.relativeResidual = max([stats.rhsReports.relativeResidual]);
                return;
            end
            A = obj.Lhs;
            stats.isSpd = obj.Decomp.isSpd;
            stats.preconditioner = obj.Decomp.type;

            if obj.Decomp.isSpd
                maxIter = obj.KrylovOpts.maxIter;
                tol = obj.KrylovOpts.tol;
                if strcmp(obj.Decomp.type, 'ichol')
                    L = obj.Decomp.data.L;
                    [x, f, rel, iter] = pcg(A, b, tol, maxIter, L, L');
                else
                    [x, f, rel, iter] = pcg(A, b, tol, maxIter);
                end

                stats.solverType = 'pcg';
                stats.flag = f;
                stats.relativeResidual = rel;
                stats.iterations = iter;
                stats.converged = (f == 0);

                assert(all(f == 0, 'all'), 'mdemat:common:linalg:LinearSolver:ConvergenceFailed', ...
                    ['PCG did not converge. ', ...
                    'Flag: %d, Relative: %g, Iterations: %d'], f, rel, iter);
                return;
            end

            maxIter = obj.KrylovOpts.maxIter;
            restart = min(obj.KrylovOpts.restart, size(A, 1));
            tol = obj.KrylovOpts.tol;

            switch obj.Decomp.type
                case 'ilu'
                    L = obj.Decomp.data.L;
                    U = obj.Decomp.data.U;
                    p = obj.Decomp.data.p;
                    Ap = A(p, p);
                    bp = b(p, :);
                    [xp, f, rel, iter] = gmres(Ap, bp, restart, tol, ...
                        maxIter, L, U);
                    x = zeros(size(b), 'like', xp);
                    x(p, :) = xp;
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

            assert(all(f == 0, 'all'), 'mdemat:common:linalg:LinearSolver:ConvergenceFailed', ...
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
            if isempty(A) || ~issparse(A) || ~ishermitian(A)
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

        function kappa = estimateCondition(A)
            % Built-in 1-norm estimator with caller RNG restored on failure too.
            state = rng;
            restore = onCleanup(@() rng(state));
            kappa = condest(A);
        end

    end
end

function mustBeSquare(A)
if size(A, 1) ~= size(A, 2)
    error('Matrix must be square');
end
end
