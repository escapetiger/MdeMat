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
        AmgOpts % Options structure for AMG preconditioner
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
                options.amgTheta{mustBeInRange(options.amgTheta, 0, 1)} = 0.35
                options.amgMaxLevels{mustBeInteger, mustBePositive} = 5
                options.amgCoarseSize{mustBeInteger, mustBePositive} = 100
                options.amgMaxCoarseRatio{mustBeInRange(options.amgMaxCoarseRatio, 0, 1)} = 0.8
                options.amgMinCoarseRatio{mustBeInRange(options.amgMinCoarseRatio, 0, 1)} = 0.05
                options.amgMaxPreStep{mustBeNumeric, mustBeNonnegative} = 1
                options.amgMaxPostStep{mustBeNumeric, mustBeNonnegative} = 1
                options.amgJacDamp{mustBeNumeric, mustBePositive} = 0.7
            end

            obj.NnzTh = options.nnzTh;
            obj.CondTh = options.condTh;

            %< Initialize options structure
            obj.KrylovOpts = struct();
            obj.KrylovOpts.maxIter = options.maxIter;
            obj.KrylovOpts.tol = options.tol;
            obj.KrylovOpts.restart = options.restart;

            %< AMG parameters
            obj.AmgOpts = struct();
            obj.AmgOpts.theta = options.amgTheta;
            obj.AmgOpts.maxLevels = options.amgMaxLevels;
            obj.AmgOpts.coarseSize = options.amgCoarseSize;
            obj.AmgOpts.maxCoarseRatio = options.amgMaxCoarseRatio;
            obj.AmgOpts.minCoarseRatio = options.amgMinCoarseRatio;
            obj.AmgOpts.maxPreStep = options.amgMaxPreStep;
            obj.AmgOpts.maxPostStep = options.amgMaxPostStep;
            obj.AmgOpts.jacDamp = options.amgJacDamp;

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

        function amgData = amgSetup(obj, A)
            amgData = cell(1, obj.AmgOpts.maxLevels);

            A0 = A;
            maxRetries = 5; %< Limit retries to avoid infinite loops

            fprintf('Building AMG hierarchy:\n');

            for level = 1:obj.AmgOpts.maxLevels
                n = size(A0, 1);

                if n <= obj.AmgOpts.coarseSize
                    break;
                end

                fprintf('  Level %d (%d dofs): ', level, n);

                amgData{level} = struct();
                amgData{level}.A = A0;
                amgData{level}.size = n;

                retryCount = 0;
                successful = false;

                while ~successful && retryCount < maxRetries
                    [P, Ac] = obj.amgCoarsen(A0);
                    nc = size(Ac, 1);
                    coarseRatio = nc / n;

                    if coarseRatio < obj.AmgOpts.minCoarseRatio
                        obj.AmgOpts.theta = min(0.9, obj.AmgOpts.theta*1.5); %< less aggressive
                        retryCount = retryCount + 1;
                        continue;
                    end

                    if coarseRatio > obj.AmgOpts.maxCoarseRatio
                        obj.AmgOpts.theta = max(0.1, obj.AmgOpts.theta*0.8); %< more aggressive
                        retryCount = retryCount + 1;
                        continue;
                    end

                    successful = true;
                    fprintf('-> %d (ratio=%.2f)\n', nc, coarseRatio);
                    amgData{level}.P = P;
                    amgData{level}.R = P';
                    A0 = Ac;
                end

                if ~successful
                    fprintf('failed, stopping\n');
                    break;
                end
            end

            amgData{level}.A = A0;
            amgData{level}.size = size(A0, 1);
            amgData = amgData(1:level);

            %< Setup smoothers
            for level = 1:length(amgData)
                A0 = amgData{level}.A;
                n0 = size(A0, 1);

                if level == length(amgData) && n0 <= 1000
                    %<  solve for small coarsest level
                    amgData{level}.smoother = '';
                    if condest(A0) < 1e12
                        [L, U, PLU] = lu(double(A0));
                        amgData{level}.L = L;
                        amgData{level}.U = U;
                        amgData{level}.PLU = PLU;
                        amgData{level}.useLU = true;
                    else
                        amgData{level}.useLU = false;
                    end
                else
                    %< Jacobi smoother
                    amgData{level}.smoother = 'jacobi';
                    d = diag(A0);
                    dSafe = d;
                    smallDiag = abs(d) < 1e-14;
                    if any(smallDiag)
                        rowNorms = sum(abs(A0), 2);
                        dSafe(smallDiag) = max(rowNorms(smallDiag), 1e-12);
                    end
                    nd = length(d);
                    nds = length(dSafe);
                    amgData{level}.DInv = spdiags(1./dSafe, 0, nds, nds);
                    amgData{level}.N = A0 - spdiags(d, 0, nd, nd);
                end
            end
        end

        function [P, Ac] = amgCoarsen(obj, A)
            %< Build coarse level using simple aggregation
            n = size(A, 1);

            %< Build strong connections
            S = obj.amgStrongConnect(A);

            %< Simple aggregation-based coarsening
            [J, nc] = obj.amgAggregate(S);

            %< Build prolongation
            P = sparse(1:n, J, 1, n, nc);

            %< Galerkin coarse operator
            R = P';
            Ac = R * A * P;

            %< Ensure coarse matrix is well-conditioned
            Ac = Ac + 1e-12 * speye(nc);
        end

        function S = amgStrongConnect(obj, A)
            %< Build strong connection matrix
            [i, j, v] = find(A);
            n = size(A, 1);

            %< Remove diagonal
            offDiag = (i ~= j);
            iOff = i(offDiag);
            jOff = j(offDiag);
            vOff = abs(v(offDiag));

            %< Find strong connections
            maxRow = accumarray(iOff, vOff, [n, 1], @max, 0);
            maxRow(maxRow == 0) = 1;
            strongMask = vOff >= obj.AmgOpts.theta * maxRow(iOff);

            S = sparse(iOff(strongMask), jOff(strongMask), 1, n, n);
            S = S | S'; %< Make symmetric
        end

        function [J, nc] = amgAggregate(~, S)
            %< Optimized simple maximal independent set aggregation
            n = size(S, 1);
            J = zeros(n, 1);
            visited = false(n, 1);
            nc = 0;

            %< Pre-extract sparse matrix structure for efficiency
            [rows, cols] = find(S);

            %< Build adjacency list for faster neighbor access
            adjList = cell(n, 1);
            for i = 1:n
                adjList{i} = [];
            end

            for k = 1:length(rows)
                i = rows(k);
                j = cols(k);
                if i ~= j % Skip diagonal
                    adjList{i}(end +1) = j;
                end
            end

            %< Process nodes in order of decreasing connectivity (greedy strategy)
            degrees = zeros(n, 1);
            for i = 1:n
                degrees(i) = length(adjList{i});
            end
            [~, order] = sort(degrees, 'descend');

            for idx = 1:n
                i = order(idx);
                if ~visited(i)
                    nc = nc + 1;
                    J(i) = nc;
                    visited(i) = true;

                    %< Add unvisited neighbors efficiently
                    neighbors = adjList{i};
                    if ~isempty(neighbors)
                        unvisitedMask = ~visited(neighbors);
                        unvisitedNeighbors = neighbors(unvisitedMask);
                        if ~isempty(unvisitedNeighbors)
                            J(unvisitedNeighbors) = nc;
                            visited(unvisitedNeighbors) = true;
                        end
                    end
                end
            end
        end

        function x = amgVCycle(obj, amgData, level, b)
            %< V-cycle implementation
            if level == length(amgData)
                %< Coarsest level solve
                x = decomposition(amgData{level}.A) \ b;
                return;
            end

            %< Initialize
            x = zeros(size(b));

            %< Pre-smoothing
            x = obj.amgSmooth(amgData{level}, b, x, obj.AmgOpts.maxPreStep);

            %< Coarse grid correction
            A = amgData{level}.A;
            r = b - A * x;

            %< Restrict
            R = amgData{level}.R;
            rc = R * r;

            %< Recursive solve
            ec = obj.amgVCycle(amgData, level+1, rc);

            %< Prolongate and correct
            P = amgData{level}.P;
            x = x + P * ec;

            %< Post-smoothing
            x = obj.amgSmooth(amgData{level}, b, x, obj.AmgOpts.maxPostStep);
        end

        function x = amgSmooth(obj, data, b, x, maxIter)
            if strcmp(data.smoother, 'jacobi')
                DInv = data.DInv;
                N = data.N;
                omega = obj.AmgOpts.jacDamp;

                for i = 1:maxIter
                    x = (1 - omega) * x + omega * DInv * (b - N * x);
                end
            end
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
