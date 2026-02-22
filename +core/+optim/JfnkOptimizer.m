classdef JfnkOptimizer < handle
    % JFNKOPTIMIZER Jacobian-Free Newton-Krylov optimizer for nonlinear systems.
    %
    %   JfnkOptimizer implements the Jacobian-Free Newton-Krylov (JFNK) method
    %   for solving nonlinear systems F(x) = 0. The method uses finite
    %   differences to approximate Jacobian-vector products and GMRES to
    %   solve the linear systems arising in Newton's method.
    %
    % Notes:
    %   The JFNK method is particularly effective for large-scale problems
    %   where forming and storing the full Jacobian matrix is impractical.
    %   The method only requires function evaluations of F(x). Based on the
    %   algorithm from JacobianFreeNewtonKrylov implementation.
    %
    % See also:
    %   JacobianFreeNewtonKrylov, gmres

    properties
        Opts            % Options structure for iterative solver
        Loss            % Function handle for the loss function F(x)
        krylovOpts      % Options structure for Krylov solvers
        lineSearchOpts  % Options structure for line search methods
    end

    methods
        function obj = JfnkOptimizer(options)
            % JFNKOPTIMIZER Construct an instance of JfnkOptimizer.
            %
            %   obj = JfnkOptimizer() creates JfnkOptimizer with default
            %   settings.
            %
            %   obj = JfnkOptimizer(Name=Value) creates JfnkOptimizer with
            %   specified options using @a Name and @a Value pairs.

            arguments
                options.maxIter{mustBeInteger, mustBePositive} = 30
                options.fTol{mustBeNumeric, mustBePositive} = 1e-6
                options.xTol{mustBeNumeric, mustBePositive} = 1e-10
                options.useMatrix logical = false
                options.epsilon{mustBeNumeric} = []
                options.mType string = "sparse"
                options.jacFun function_handle = function_handle.empty
                options.krylov string = "gmres"
                options.krylovOpts struct = struct()
                options.lineSearch string = "none"
                options.lineSearchOpts struct = struct()
                options.equilibrate logical = false
                options.reorderingMethod string = "dissect"
            end

            %< Initialize options structure with all JFNK parameters
            obj.Opts = struct();
            obj.Opts.maxIter = options.maxIter;
            obj.Opts.fTol = options.fTol;
            obj.Opts.xTol = options.xTol;
            obj.Opts.useMatrix = options.useMatrix;
            obj.Opts.epsilon = options.epsilon;
            obj.Opts.mType = options.mType;
            obj.Opts.jacFun = options.jacFun;
            obj.Opts.krylov = options.krylov;
            obj.Opts.lineSearch = options.lineSearch;
            obj.Opts.equilibrate = options.equilibrate;
            obj.Opts.reorderingMethod = options.reorderingMethod;

            %< Initialize separate options structures
            obj.krylovOpts = options.krylovOpts;
            obj.lineSearchOpts = options.lineSearchOpts;

            %< Initialize loss function
            obj.Loss = [];
        end

        function [x, stats] = solve(obj, f, x0)
            % SOLVE Solve the nonlinear system F(x) = 0 using JFNK method.
            %
            %   [x, stats] = solve(f, x0) solves the nonlinear system @a f(x) = 0
            %   starting from initial guess @a x0. Returns the solution @a x and
            %   statistics structure @a stats with convergence information.
            %
            %   The function @a f should accept a column vector and return a column
            %   vector of the same size representing the residual.

            arguments
                obj core.optim.JfnkOptimizer
                f function_handle
                x0 {mustBeNumeric, mustBeVector}
            end

            obj.Loss = f;
            x = x0(:);

            %< Initialize epsilon if not provided
            if isempty(obj.Opts.epsilon)
                obj.Opts.epsilon = sum(sqrt(eps) * (1 + abs(x0)));
            end

            %< Set up default jacFun if not provided
            if isempty(obj.Opts.jacFun)
                obj.Opts.jacFun = @(x) obj.jacmat(f, x);
            end

            %< Initialize statistics
            stats = struct();
            stats.nlIter = 0;
            stats.fNorm = zeros(obj.Opts.maxIter + 1, 1);
            stats.xNorm = zeros(obj.Opts.maxIter + 1, 1);
            stats.convergedId = 0;
            stats.convergedMessage = 'Not Started';
            stats.krylovStats = struct();
            stats.krylovStats.exitStatusId = zeros(obj.Opts.maxIter, 1);
            stats.krylovStats.linearRes = zeros(obj.Opts.maxIter, 1);
            stats.krylovStats.linearSteps = zeros(obj.Opts.maxIter, 1);
            stats.krylovStats.resVec = cell(obj.Opts.maxIter, 1);

            F = f(x);
            stats.fNorm(1) = norm(F);
            stats.xNorm(1) = norm(-x0);

            %< Main JFNK iteration loop following JacobianFreeNewtonKrylov algorithm
            for iter = 1:obj.Opts.maxIter
                F = f(x);

                %< Check convergence based on function tolerance
                if norm(F) <= obj.Opts.fTol
                    stats.convergedId = 1;
                    stats.convergedMessage = 'Residual below function tolerance';
                    break;
                end

                %< Solve linear system using matrix or matrix-free approach
                if obj.Opts.useMatrix
                    %< Use full Jacobian matrix
                    J = obj.Opts.jacFun(x);

                    if obj.Opts.equilibrate
                        %< Apply matrix equilibration and reordering
                        [P, R, C] = equilibrate(J);
                        JNew = R * P * J * C;
                        FNew = R * P * F;

                        %< Apply reordering
                        if strcmp(obj.Opts.reorderingMethod, 'dissect')
                            q = dissect(JNew);
                        elseif strcmp(obj.Opts.reorderingMethod, 'amd')
                            q = amd(JNew);
                        elseif strcmp(obj.Opts.reorderingMethod, 'symrcm')
                            q = symrcm(JNew);
                        else
                            q = 1:size(JNew, 1);
                        end

                        JNew = JNew(q, q);
                        FNew = FNew(q);

                        %< Solve reordered system
                        [deltaNew, linFlag, linRelRes, linIter, linResVec] = obj.krylov(JNew, -FNew, x);

                        %< Undo reordering and equilibration
                        delta = zeros(size(deltaNew));
                        delta(q) = deltaNew;
                        delta = C * delta;
                    else
                        %< Direct solve without equilibration
                        [delta, linFlag, linRelRes, linIter, linResVec] = obj.krylov(J, -F, x);
                    end
                else
                    %< Use Jacobian-free approach
                    jvProduct = @(v) obj.jacvec(v, x, F);
                    [delta, linFlag, linRelRes, linIter, linResVec] = obj.krylov(jvProduct, -F, x);
                end

                %< Store Krylov statistics
                stats.krylovStats.exitStatusId(iter) = linFlag;
                stats.krylovStats.linearRes(iter) = linRelRes;
                stats.krylovStats.linearSteps(iter) = max(linIter);
                stats.krylovStats.resVec{iter} = linResVec;

                %< Apply line search
                if strcmp(obj.Opts.lineSearch, 'none')
                    %< Standard Newton step
                    step = delta;
                elseif strcmp(obj.Opts.lineSearch, 'backtracking')
                    %< Backtracking line search
                    step = obj.backtrackLineSearch(f, x, F, delta);
                else
                    step = delta;
                end

                %< Update solution
                xOld = x;
                x = x + step;

                %< Update statistics
                stats.nlIter = stats.nlIter + 1;
                F = f(x);
                stats.fNorm(iter + 1) = norm(F);
                stats.xNorm(iter + 1) = norm(x - xOld);

                %< Check solution tolerance
                if stats.xNorm(iter + 1) <= obj.Opts.xTol
                    stats.convergedId = 2;
                    stats.convergedMessage = 'Difference in x below tolerance';
                    break;
                end
            end

            %< Trim statistics arrays
            stats.fNorm = stats.fNorm(1:stats.nlIter + 1);
            stats.xNorm = stats.xNorm(1:stats.nlIter + 1);
            stats.krylovStats.exitStatusId = stats.krylovStats.exitStatusId(1:stats.nlIter);
            stats.krylovStats.linearRes = stats.krylovStats.linearRes(1:stats.nlIter);
            stats.krylovStats.linearSteps = stats.krylovStats.linearSteps(1:stats.nlIter);
            stats.krylovStats.resVec(stats.nlIter + 1:end) = [];

            %< Final convergence check
            if stats.convergedId == 0
                stats.convergedMessage = sprintf('Newton method failed to converge after %d iterations', stats.nlIter);
                stats.convergedId = -1;
            end
        end
    end

    methods (Access = private)
        function [delta, linFlag, linRelRes, linIter, linResVec] = krylov(obj, A, b, x)
            % CALLKRYLOVSOLVER Call appropriate Krylov solver with preconditioners.

            %< Set up preconditioners
            M1 = [];
            M2 = [];

            if isfield(obj.krylovOpts, 'M1') && ~isempty(obj.krylovOpts.M1)
                if isa(obj.krylovOpts.M1, 'function_handle')
                    userData = struct();
                    if isfield(obj.krylovOpts, 'userData')
                        userData = obj.krylovOpts.userData;
                    end
                    M1 = @(b) obj.krylovOpts.M1(b, obj.Opts.jacFun, x, userData);
                else
                    M1 = obj.krylovOpts.M1;
                end
            end

            if isfield(obj.krylovOpts, 'M2') && ~isempty(obj.krylovOpts.M2)
                if isa(obj.krylovOpts.M2, 'function_handle')
                    userData = struct();
                    if isfield(obj.krylovOpts, 'userData')
                        userData = obj.krylovOpts.userData;
                    end
                    M2 = @(b) obj.krylovOpts.M2(b, obj.Opts.jacFun, x, userData);
                else
                    M2 = obj.krylovOpts.M2;
                end
            end

            %< Call appropriate Krylov solver
            switch lower(obj.Opts.krylov)
                case 'gmres'
                    obj.setDefaultKrylovOpts('gmres', length(x));
                    [delta, linFlag, linRelRes, linIter, linResVec] = gmres(A, b, ...
                        obj.krylovOpts.restart, obj.krylovOpts.tol, ...
                        obj.krylovOpts.maxIt, M1, M2, []);

                case 'bicgstab'
                    obj.setDefaultKrylovOpts('bicgstab', length(x));
                    [delta, linFlag, linRelRes, linIter, linResVec] = bicgstab(A, b, ...
                        obj.krylovOpts.tol, obj.krylovOpts.maxIt, M1, M2, []);

                case 'pcg'
                    obj.setDefaultKrylovOpts('pcg', length(x));
                    [delta, linFlag, linRelRes, linIter, linResVec] = pcg(A, b, ...
                        obj.krylovOpts.tol, obj.krylovOpts.maxIt, M1, M2, []);

                case 'minres'
                    obj.setDefaultKrylovOpts('minres', length(x));
                    [delta, linFlag, linRelRes, linIter, linResVec] = minres(A, b, ...
                        obj.krylovOpts.tol, obj.krylovOpts.maxIt, M1, M2, []);

                case 'symmlq'
                    obj.setDefaultKrylovOpts('symmlq', length(x));
                    [delta, linFlag, linRelRes, linIter, linResVec] = symmlq(A, b, ...
                        obj.krylovOpts.tol, obj.krylovOpts.maxIt, M1, M2, []);

                case 'bicgstabl'
                    obj.setDefaultKrylovOpts('bicgstabl', length(x));
                    [delta, linFlag, linRelRes, linIter, linResVec] = bicgstabl(A, b, ...
                        obj.krylovOpts.tol, obj.krylovOpts.maxIt, M1, M2, []);

                case 'cgs'
                    obj.setDefaultKrylovOpts('cgs', length(x));
                    [delta, linFlag, linRelRes, linIter, linResVec] = cgs(A, b, ...
                        obj.krylovOpts.tol, obj.krylovOpts.maxIt, M1, M2, []);

                case 'tfqmr'
                    obj.setDefaultKrylovOpts('tfqmr', length(x));
                    [delta, linFlag, linRelRes, linIter, linResVec] = tfqmr(A, b, ...
                        obj.krylovOpts.tol, obj.krylovOpts.maxIt, M1, M2, []);

                otherwise
                    % Default to GMRES
                    obj.setDefaultKrylovOpts('gmres', length(x));
                    [delta, linFlag, linRelRes, linIter, linResVec] = gmres(A, b, ...
                        obj.krylovOpts.restart, obj.krylovOpts.tol, ...
                        obj.krylovOpts.maxIt, M1, M2, []);
            end
        end

        function setDefaultKrylovOpts(obj, solver, n)
            % SETDEFAULTKRYLOVOPTS Set default options for Krylov solvers.
            switch lower(solver)
                case 'gmres'
                    if ~isfield(obj.krylovOpts, 'restart')
                        obj.krylovOpts.restart = min(30, n);
                    end
                    if ~isfield(obj.krylovOpts, 'tol')
                        obj.krylovOpts.tol = 1e-4;
                    end
                    if ~isfield(obj.krylovOpts, 'maxIt')
                        obj.krylovOpts.maxIt = min(n, 20);
                    end

                case {'bicgstab', 'pcg', 'minres', 'symmlq', 'bicgstabl', 'cgs', 'tfqmr'}
                    if ~isfield(obj.krylovOpts, 'tol')
                        obj.krylovOpts.tol = 1e-4;
                    end
                    if ~isfield(obj.krylovOpts, 'maxIt')
                        obj.krylovOpts.maxIt = min(2*n, 20);
                    end
            end
        end

        function jv = jacvec(obj, v, x, F)
            % JACOBIANVECTORPRODUCT Approximate Jacobian-vector product using finite differences.
            %
            %   jv = jacobianVectorProduct(v, x, F) computes the approximate product
            %   J*v where J is the Jacobian of the loss function at point @a x.
            %   Uses the algorithm from JacobianOperator implementation.

            n = length(x);

            %< Compute finite difference step size following JFNK algorithm
            if isempty(obj.Opts.epsilon)
                if norm(v, 2) > eps
                    s = sum(sqrt(eps) * (1 + abs(x)));
                    p = s / (n * norm(v, 2));
                else
                    s = sum(sqrt(eps) * (1 + abs(x)));
                    p = s / n;
                end
            else
                p = obj.Opts.epsilon;
            end

            %< Compute finite difference approximation using JacobianOperator algorithm
            xp = x + p * v;
            Fp = obj.Loss(xp);
            jv = (Fp - F) / p;
        end

        function J = jacmat(obj, f, x)
            % NUMERICALJACOBIAN Compute full Jacobian matrix using finite differences.
            %
            %   J = numericalJacobianMatrix(f, x) computes the full Jacobian matrix
            %   of function @a f at point @a x using column-by-column finite differences.
            %   Based on NumericalJacobianMatrix implementation.

            F = f(x);
            n = length(x);
            m = length(F);
            I = speye(m, n);

            %< Use epsilon from options or compute default
            if isempty(obj.Opts.epsilon)
                epsilon = sum(sqrt(eps) * (1 + abs(x)));
            else
                epsilon = obj.Opts.epsilon;
            end

            if strcmp(obj.Opts.mType, 'dense')
                %< Preallocate dense matrix
                J = zeros(m, n);
                for i = 1:n
                    J(:, i) = (f(x + I(:, i) * epsilon) - F) / epsilon;
                end
            else
                %< Preallocate sparse matrix
                J = spalloc(m, n, ceil(max([m*n/100, 10*m, 10*n])));
                for i = 1:n
                    J = J + sparse(1:m, i*ones(1,m), ((f(x + I(:, i) * epsilon) - F) / epsilon), m, n);
                end
            end
        end

        function s = backtrackLineSearch(obj, f, x, F, s)
            % BACKTRACKINESEARCH Perform backtracking line search.
            %
            %   s = backtrackLineSearch(f, x, F, s) performs backtracking line search
            %   to find an appropriate step size. Returns the scaled search direction.
            %   Based on BacktrackLinesearch implementation.

            %< Initial step size
            a = 1;

            %< Get line search options with defaults
            if isfield(obj.lineSearchOpts, 'minStep')
                aMin = obj.lineSearchOpts.minStep;
            else
                aMin = 1e-3;
            end

            if isfield(obj.lineSearchOpts, 'contractionRate')
                tau = obj.lineSearchOpts.contractionRate;
            else
                tau = 0.5;
            end

            %< Backtracking loop
            while (norm(f(x + a * s)) > norm(F) && a > aMin)
                a = a * tau;
            end

            %< Return scaled search direction
            s = s * a;
        end
    end
end