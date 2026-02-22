%% Test AMG with 2D Poisson Matrix
% This demonstrates ideal AMG performance on a model problem

clear; clc; close all;

%% Create test problem
fprintf('=== CREATING 2D POISSON TEST PROBLEM ===\n');
n_grid = 50;  % Grid size (will create ~2400 unknowns)
[A, b] = create_amg_test_matrix(n_grid);

%% AMG Parameters (optimized for Poisson)
amg_params = struct();
amg_params.theta = 0.25;                 % Standard for Poisson
amg_params.max_levels = 10;              % Allow more levels
amg_params.coarse_size = 50;             % Smaller coarse size
amg_params.max_coarse_ratio = 0.5;       % Good coarsening
amg_params.min_coarse_ratio = 0.1;       % Reasonable bounds

% GMRES Parameters  
gmres_params = struct();
gmres_params.restart = 20;               % Small restart for Poisson
gmres_params.max_iter = 100;             % Should converge fast
gmres_params.tol = 1e-10;                % Tight tolerance

% Smoothing Parameters
smooth_params = struct();
smooth_params.nu1 = 1;                   % Standard pre-smoothing
smooth_params.nu2 = 1;                   % Standard post-smoothing  
smooth_params.omega = 2/3;               % Optimal Jacobi for Poisson

fprintf('\nAMG Parameters:\n');
fprintf('  theta=%.3f, max_levels=%d, coarse_size=%d\n', ...
    amg_params.theta, amg_params.max_levels, amg_params.coarse_size);

%% Setup AMG (using your existing function)
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