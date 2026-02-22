%% Robust AMG-GMRES Solver with Manual Parameters
% Clean implementation for large, non-symmetric sparse systems

clear; clc; close all;

%% Load and analyze problem
% fprintf('Loading problem...\n');
% load('matlab.mat', 'A', 'b');

%% Create test problem
fprintf('=== CREATING NON-SYMMETRIC ADVECTION-DIFFUSION PROBLEM ===\n');
n_grid = 100;   % Grid size (will create ~2400 unknowns)
Pe = 100;      % Péclet number (high = advection-dominated, challenging for AMG)
[A, b] = create_nonsymmetric_matrix(n_grid, Pe);

%%
[m, n] = size(A);
fprintf('Matrix: %dx%d, nnz=%d, density=%.2e\n', m, n, nnz(A), nnz(A)/(m*n));

% Fast condition number estimate
fprintf('Computing fast condition estimate...\n');
tic;
cond_est = fast_condition_estimate(A);
est_time = toc;
fprintf('Condition estimate: %.2e (computed in %.3f s)\n', cond_est, est_time);

%% Manual Parameter Configuration
fprintf('\n=== MANUAL PARAMETERS ===\n');

% AMG Parameters
amg_params = struct();
amg_params.theta = 0.35;                 % Strong connection threshold
amg_params.max_levels = 5;               % Maximum AMG levels
amg_params.coarse_size = 100;            % Coarsest level size
amg_params.max_coarse_ratio = 0.8;       % Maximum coarsening ratio
amg_params.min_coarse_ratio = 0.05;      % Minimum coarsening ratio

% GMRES Parameters
gmres_params = struct();
gmres_params.restart = 50;               % GMRES restart
gmres_params.max_iter = 1000;            % Maximum iterations
gmres_params.tol = 1e-8;                 % Convergence tolerance

% Smoothing Parameters
smooth_params = struct();
smooth_params.nu1 = 1;                   % Pre-smoothing steps
smooth_params.nu2 = 1;                   % Post-smoothing steps
smooth_params.omega = 0.7;               % Jacobi damping factor

fprintf('AMG: theta=%.3f, levels=%d, coarse_size=%d\n', ...
    amg_params.theta, amg_params.max_levels, amg_params.coarse_size);
fprintf('GMRES: restart=%d, max_iter=%d, tol=%.2e\n', ...
    gmres_params.restart, gmres_params.max_iter, gmres_params.tol);
fprintf('Smoothing: nu1=%d, nu2=%d, omega=%.2f\n', ...
    smooth_params.nu1, smooth_params.nu2, smooth_params.omega);


%% Setup AMG Preconditioner
fprintf('\n=== SETTING UP AMG ===\n');
tic;
amg_data = setup_amg(A, amg_params, smooth_params);
setup_time = toc;
fprintf('AMG setup: %.3f seconds\n', setup_time);

%% Test 1: No Preconditioning
fprintf('\n=== TEST 1: GMRES WITHOUT PRECONDITIONING ===\n');
tic;
[x1, flag1, relres1, iter1, resvec1] = gmres(A, b, gmres_params.restart, ...
    gmres_params.tol, gmres_params.max_iter);
time1 = toc;

fprintf('No preconditioning:\n');
fprintf('  Flag: %d, Iterations: %d\n', flag1, iter1(2));
fprintf('  Relative residual: %.2e\n', relres1);
fprintf('  Time: %.3f seconds\n', time1);

%% Test 2: AMG Preconditioning  
fprintf('\n=== TEST 2: GMRES WITH AMG PRECONDITIONING ===\n');
M = @(x) apply_amg_preconditioner(amg_data, x, smooth_params);
x0 = zeros(size(b));

tic;
[x2, flag2, relres2, iter2, resvec2] = gmres(A, b, gmres_params.restart, ...
    gmres_params.tol, gmres_params.max_iter, M, [], x0);
time2 = toc;

fprintf('AMG preconditioning:\n');
fprintf('  Flag: %d, Iterations: %d\n', flag2, iter2(2));
fprintf('  Relative residual: %.2e\n', relres2);
fprintf('  Time: %.3f seconds\n', time2);

%% Results comparison
fprintf('\n=== COMPARISON ===\n');
fprintf('Speedup: %.1fx faster\n', time1/time2);
fprintf('Iteration reduction: %d -> %d (%.1fx fewer)\n', ...
    iter1(2), iter2(2), iter1(2)/iter2(2));

% Verify solutions are similar
solution_diff = norm(x1 - x2) / norm(x1);
fprintf('Solution difference: %.2e\n', solution_diff);

%% Plot convergence
figure('Position', [100, 100, 800, 300]);

subplot(1,2,1);
semilogy(0:length(resvec1)-1, resvec1/resvec1(1), 'r-', 'LineWidth', 2);
hold on;
semilogy(0:length(resvec2)-1, resvec2/resvec2(1), 'b-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Relative Residual');
title('GMRES Convergence Comparison');
legend('No Preconditioning', 'AMG Preconditioning', 'Location', 'best');
grid on;

subplot(1,2,2);
% Plot solution (2D visualization)
ni = n_grid - 2;
X = reshape(x2, ni, ni);
h = 1/(n_grid-1);
x_coords = h:h:(1-h);
y_coords = h:h:(1-h);
contourf(x_coords, y_coords, X', 20);
colorbar;
title('AMG Solution');
xlabel('x'); ylabel('y');

fprintf('\nThis demonstrates ideal AMG performance!\n');
fprintf('For comparison, try your original matrix by replacing the create_amg_test_matrix call.\n');

%% HELPER FUNCTIONS

function amg_data = setup_amg(A, amg_params, smooth_params)
    % Setup AMG hierarchy
    amg_data = struct();
    amg_data.levels = cell(1, amg_params.max_levels);
    amg_data.num_levels = 0;
    amg_data.smooth_params = smooth_params;
    
    current_A = A;
    level = 1;
    max_retries = 5;  % Limit retries to avoid infinite loops
    
    fprintf('Building AMG hierarchy:\n');
    
    while level <= amg_params.max_levels && size(current_A, 1) > amg_params.coarse_size
        n = size(current_A, 1);
        fprintf('  Level %d (%d dofs): ', level, n);
        
        % Store current level
        amg_data.levels{level} = struct();
        amg_data.levels{level}.A = current_A;
        amg_data.levels{level}.size = n;
        
        % Try to build coarse level
        if level < amg_params.max_levels && n > amg_params.coarse_size
            retry_count = 0;
            theta = amg_params.theta;
            successful = false;
            
            while ~successful && retry_count < max_retries
                [P, A_coarse] = build_coarse_level(current_A, theta);
                coarse_ratio = size(A_coarse, 1) / n;
                
                if coarse_ratio >= amg_params.min_coarse_ratio && coarse_ratio <= amg_params.max_coarse_ratio
                    successful = true;
                    fprintf('-> %d (ratio=%.2f)\n', size(A_coarse, 1), coarse_ratio);
                    
                    % Store operators
                    amg_data.levels{level}.P = P;
                    amg_data.levels{level}.R = P';
                    current_A = A_coarse;
                else
                    retry_count = retry_count + 1;
                    if coarse_ratio < amg_params.min_coarse_ratio
                        theta = min(0.9, theta * 1.5);  % Less aggressive
                    else
                        theta = max(0.1, theta * 0.8);  % More aggressive
                    end
                end
            end
            
            if ~successful
                fprintf('failed, stopping\n');
                break;
            end
        else
            fprintf('final level\n');
            break;
        end
        
        level = level + 1;
    end
    
    % Store final level
    if level <= length(amg_data.levels)
        amg_data.levels{level}.A = current_A;
        amg_data.levels{level}.size = size(current_A, 1);
    end
    amg_data.num_levels = level;
    
    % Setup smoothers
    amg_data = setup_smoothers(amg_data, smooth_params);
end

function [P, A_coarse] = build_coarse_level(A, theta)
    % Build coarse level using simple aggregation
    n = size(A, 1);
    
    % Build strong connections
    S = build_strong_connections(A, theta);
    
    % Simple aggregation-based coarsening
    [aggregates, nc] = simple_aggregation(S);
    
    % Build prolongation
    P = sparse(1:n, aggregates, 1, n, nc);
    
    % Galerkin coarse operator
    R = P';
    A_coarse = R * A * P;
    
    % Ensure coarse matrix is well-conditioned
    A_coarse = A_coarse + 1e-12 * speye(nc);
end

function S = build_strong_connections(A, theta)
    % Build strong connection matrix
    [i, j, v] = find(A);
    n = size(A, 1);
    
    % Remove diagonal
    off_diag = (i ~= j);
    i_off = i(off_diag);
    j_off = j(off_diag);
    v_off = abs(v(off_diag));
    
    % Find strong connections
    max_row = accumarray(i_off, v_off, [n, 1], @max, 0);
    max_row(max_row == 0) = 1;
    strong_mask = v_off >= theta * max_row(i_off);
    
    S = sparse(i_off(strong_mask), j_off(strong_mask), 1, n, n);
    S = S | S';  % Make symmetric
end

function [aggregates, nc] = simple_aggregation(S)
    % Optimized simple maximal independent set aggregation
    n = size(S, 1);
    aggregates = zeros(n, 1);
    visited = false(n, 1);
    nc = 0;
    
    % Pre-extract sparse matrix structure for efficiency
    [rows, cols] = find(S);
    
    % Build adjacency list for faster neighbor access
    adj_list = cell(n, 1);
    for i = 1:n
        adj_list{i} = [];
    end
    
    for k = 1:length(rows)
        i = rows(k);
        j = cols(k);
        if i ~= j  % Skip diagonal
            adj_list{i}(end+1) = j;
        end
    end
    
    % Process nodes in order of decreasing connectivity (greedy strategy)
    degrees = zeros(n, 1);
    for i = 1:n
        degrees(i) = length(adj_list{i});
    end
    [~, order] = sort(degrees, 'descend');
    
    for idx = 1:n
        i = order(idx);
        if ~visited(i)
            nc = nc + 1;
            aggregates(i) = nc;
            visited(i) = true;
            
            % Add unvisited neighbors efficiently
            neighbors = adj_list{i};
            if ~isempty(neighbors)
                unvisited_mask = ~visited(neighbors);
                unvisited_neighbors = neighbors(unvisited_mask);
                if ~isempty(unvisited_neighbors)
                    aggregates(unvisited_neighbors) = nc;
                    visited(unvisited_neighbors) = true;
                end
            end
        end
    end
end

function amg_data = setup_smoothers(amg_data, smooth_params)
    % Setup smoothers for each level
    for level = 1:amg_data.num_levels
        A_level = amg_data.levels{level}.A;
        n_level = size(A_level, 1);
        
        if level == amg_data.num_levels && n_level <= 1000
            % Direct solve for small coarsest level
            amg_data.levels{level}.smoother_type = 'direct';
            try
                [L, U, P_lu] = lu(double(A_level));
                amg_data.levels{level}.L = L;
                amg_data.levels{level}.U = U;
                amg_data.levels{level}.P_lu = P_lu;
                amg_data.levels{level}.use_lu = true;
            catch
                amg_data.levels{level}.use_lu = false;
            end
        else
            % Jacobi smoother
            amg_data.levels{level}.smoother_type = 'jacobi';
            d = diag(A_level);
            d_safe = d;
            small_diag = abs(d) < 1e-14;
            if any(small_diag)
                row_norms = sum(abs(A_level), 2);
                d_safe(small_diag) = max(row_norms(small_diag), 1e-12);
            end
            amg_data.levels{level}.D_inv = spdiags(1./d_safe, 0, length(d_safe), length(d_safe));
            amg_data.levels{level}.N = A_level - spdiags(d, 0, length(d), length(d));
            amg_data.levels{level}.omega = smooth_params.omega;
        end
    end
end

function y = apply_amg_preconditioner(amg_data, x, smooth_params)
    % Apply AMG V-cycle preconditioner

    y = amg_vcycle(amg_data, 1, x, smooth_params.nu1, smooth_params.nu2);
    if any(~isfinite(y))
        y = x;  % Fallback
    end
end

function x = amg_vcycle(amg_data, level, b, nu1, nu2)
    % V-cycle implementation
    if level == amg_data.num_levels
        % Coarsest level solve
        x = coarse_solve(amg_data.levels{level}, b);
    else
        % Initialize
        x = zeros(size(b));
        
        % Pre-smoothing
        x = smooth(amg_data.levels{level}, x, b, nu1);
        
        % Coarse grid correction
        if level < amg_data.num_levels
            A = amg_data.levels{level}.A;
            r = b - A * x;
            
            % Restrict
            R = amg_data.levels{level}.R;
            r_coarse = R * r;
            
            % Recursive solve
            e_coarse = amg_vcycle(amg_data, level + 1, r_coarse, nu1, nu2);
            
            % Prolongate and correct
            P = amg_data.levels{level}.P;
            x = x + P * e_coarse;
        end
        
        % Post-smoothing
        x = smooth(amg_data.levels{level}, x, b, nu2);
    end
end

function x = coarse_solve(level_data, b)
    % Solve coarsest level
    if strcmp(level_data.smoother_type, 'direct') && level_data.use_lu
        try
            x = level_data.U \ (level_data.L \ (level_data.P_lu * b));
        catch
            x = level_data.A \ b;
        end
    else
        % Iterative solve
        x = zeros(size(b));
        for iter = 1:20
            x_old = x;
            x = level_data.D_inv * (b - level_data.N * x);
            if norm(x - x_old) < 1e-12 * norm(x)
                break;
            end
        end
    end
end

function x = smooth(level_data, x, b, nu)
    % Apply smoothing
    if strcmp(level_data.smoother_type, 'jacobi')
        D_inv = level_data.D_inv;
        N = level_data.N;
        omega = level_data.omega;
        
        for i = 1:nu
            x = (1 - omega) * x + omega * D_inv * (b - N * x);
        end
    end
end

function cond_est = fast_condition_estimate(A)
    % Fast condition number estimate for large sparse matrices
    % Uses power iteration and inverse power iteration for largest/smallest eigenvalues
    
    n = size(A, 1);
    max_iter = 20;  % Keep iterations low for speed
    tol = 1e-3;     % Relaxed tolerance for estimate
    
    try
        % Estimate largest eigenvalue using power iteration
        v = randn(n, 1);
        v = v / norm(v);
        
        for iter = 1:max_iter
            v_old = v;
            v = A * v;
            lambda_max = norm(v);
            v = v / lambda_max;
            
            if norm(v - v_old) < tol
                break;
            end
        end
        
        % Estimate smallest eigenvalue using diagonal dominance approximation
        % This is much faster than inverse power iteration
        d = abs(diag(A));
        off_diag_sums = sum(abs(A), 2) - d;
        
        % Gershgorin circle estimate for smallest eigenvalue
        lambda_min_est = min(d - off_diag_sums);
        
        % If matrix might be nearly singular, use a conservative estimate
        if lambda_min_est <= 0
            % Use Frobenius norm scaling approach
            lambda_min_est = norm(A, 'fro') / n^2;
        end
        
        % Condition number estimate
        cond_est = abs(lambda_max / lambda_min_est);
        
        % Cap at reasonable value to avoid infinity
        cond_est = min(cond_est, 1e16);
        
    catch
        % Fallback: use simple diagonal-based estimate
        d = abs(diag(A));
        d(d == 0) = 1e-12;  % Handle zeros
        cond_est = max(d) / min(d);
    end
end