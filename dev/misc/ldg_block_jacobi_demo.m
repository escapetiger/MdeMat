%% Block Jacobi Preconditioner for 3D LDG Augmented Systems
% Test script for 3D Poisson with periodic BC using 3rd-order LDG method
% Handles augmented/saddle-point systems

clear; clc; close all;

%% Load LDG matrix
fprintf('=== LOADING LDG MATRIX FROM FILE ===\n');
load('matlab.mat', 'A', 'b');
[m, n] = size(A);
fprintf('Matrix: %dx%d, nnz=%d, density=%.2e\n', m, n, nnz(A), nnz(A)/(m*n));

% Check symmetry
symmetry_measure = norm(A - A', 'fro') / norm(A, 'fro');
fprintf('  Non-symmetry measure: %.2e\n', symmetry_measure);
fprintf('  Matrix type: %s\n', symmetry_measure < 1e-12 );

%% Block structure detection
fprintf('\n=== DETECTING BLOCK STRUCTURE ===\n');

% For 3rd-order LDG in 3D: 27 DOFs per element (could be augmented)
% Augmented systems may have additional constraint variables
possible_block_sizes = [9, 18, 27, 36, 54, 81];
best_block_size = 27;  % Default assumption for 3D LDG

% Check if this might be a saddle-point system
fprintf('Analyzing system structure:\n');
% Simple check: look for zero diagonal blocks (typical of saddle-point)
diag_A = diag(A);
zero_diag_count = sum(abs(diag_A) < 1e-14);

if zero_diag_count > 0.1 * length(diag_A)
    fprintf('  System type: Likely saddle-point/augmented\n');
    fprintf('  Recommendation: Use specialized saddle-point preconditioners\n');
else
    fprintf('  System type: Standard elliptic\n');
end

% Auto-detect block size by analyzing matrix structure
for block_size = possible_block_sizes
    if mod(m, block_size) == 0
        num_blocks = int(m / block_size);
        % Check if diagonal blocks are reasonably dense
        sample_block = A(1:block_size, 1:block_size);
        block_density = nnz(sample_block) / block_size^2;
        
        fprintf('  Block size %d: %d blocks, density=%.3f\n', ...
            block_size, num_blocks, block_density);
        
        if block_density > 0.3  % Reasonable density for a good block
            best_block_size = block_size;
        end
    end
end

block_size = best_block_size;
num_blocks = m / block_size;
fprintf('Selected block size: %d (%d blocks)\n', block_size, num_blocks);

%% GMRES Parameters
gmres_params = struct();
gmres_params.restart = 50;
gmres_params.max_iter = 500;
gmres_params.tol = 1e-8;

fprintf('\nGMRES parameters: restart=%d, max_iter=%d, tol=%.2e\n', ...
    gmres_params.restart, gmres_params.max_iter, gmres_params.tol);

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

%% Test 2: Block Jacobi Preconditioning
fprintf('\n=== TEST 2: GMRES WITH BLOCK JACOBI PRECONDITIONING ===\n');

% Setup Block Jacobi preconditioner
fprintf('Setting up Block Jacobi preconditioner...\n');
tic;
[M_block, block_data] = setup_block_jacobi(A, block_size);
setup_time = toc;
fprintf('Block Jacobi setup: %.3f seconds\n', setup_time);
fprintf('  Block factorizations: %d successful, %d failed\n', ...
    block_data.successful, block_data.failed);

% Solve with Block Jacobi preconditioning
x0 = zeros(size(b));
tic;
[x2, flag2, relres2, iter2, resvec2] = gmres(A, b, gmres_params.restart, ...
    gmres_params.tol, gmres_params.max_iter, M_block, [], x0);
time2 = toc;

fprintf('Block Jacobi preconditioning:\n');
fprintf('  Flag: %d, Iterations: %d\n', flag2, iter2(2));
fprintf('  Relative residual: %.2e\n', relres2);
fprintf('  Time: %.3f seconds (including %.3f setup)\n', time2, setup_time);

%% Test 3: ILU Preconditioning (for comparison)
fprintf('\n=== TEST 3: GMRES WITH ILU PRECONDITIONING ===\n');
try
    tic;
    setup_start = tic;
    [L, U] = ilu(A, struct('type', 'ilutp', 'droptol', 1e-4));
    ilu_setup_time = toc(setup_start);
    fprintf('ILU setup: %.3f seconds\n', ilu_setup_time);

    M1 = @(x) L \ x;
    M2 = @(x) U \ x;

    [x3, flag3, relres3, iter3, resvec3] = gmres(A, b, gmres_params.restart, ...
        gmres_params.tol, gmres_params.max_iter, M1, M2, x0);
    time3 = toc;

    fprintf('ILU preconditioning:\n');
    fprintf('  Flag: %d, Iterations: %d\n', flag3, iter3(2));
    fprintf('  Relative residual: %.2e\n', relres3);
    fprintf('  Time: %.3f seconds (including %.3f setup)\n', time3, ilu_setup_time);
    
    ilu_success = true;
catch ME
    fprintf('ILU failed: %s\n', ME.message);
    ilu_success = false;
end

%% Results comparison
fprintf('\n=== COMPARISON ===\n');
fprintf('Method           | Iter | Time   | Speedup | Residual\n');
fprintf('-----------------|------|--------|---------|----------\n');
fprintf('No Precond       | %4d | %6.3f |   1.0x  | %.2e\n', ...
    iter1(2), time1, relres1);
fprintf('Block Jacobi     | %4d | %6.3f |   %.1fx  | %.2e\n', ...
    iter2(2), time2, time1/time2, relres2);
if ilu_success
    fprintf('ILU              | %4d | %6.3f |   %.1fx  | %.2e\n', ...
        iter3(2), time3, time1/time3, relres3);
end

% Verify solutions are similar (if multiple converged)
if flag1 == 0 && flag2 == 0
    solution_diff = norm(x1 - x2) / norm(x1);
    fprintf('\nSolution difference (No Precond vs Block Jacobi): %.2e\n', solution_diff);
end

%% Plot convergence
figure('Position', [100, 100, 1000, 400]);

subplot(1,2,1);
semilogy(0:length(resvec1)-1, resvec1/resvec1(1), 'r-', 'LineWidth', 2);
hold on;
semilogy(0:length(resvec2)-1, resvec2/resvec2(1), 'b-', 'LineWidth', 2);
if ilu_success
    semilogy(0:length(resvec3)-1, resvec3/resvec3(1), 'g-', 'LineWidth', 2);
    legend('No Preconditioning', 'Block Jacobi', 'ILU', 'Location', 'best');
else
    legend('No Preconditioning', 'Block Jacobi', 'Location', 'best');
end
xlabel('Iteration');
ylabel('Relative Residual');
title('GMRES Convergence Comparison');
grid on;

subplot(1,2,2);
% Visualize block structure
spy(A(1:min(200,m), 1:min(200,n)));
title(sprintf('Matrix Sparsity Pattern (first %dx%d)', min(200,m), min(200,n)));
% Highlight block boundaries
hold on;
for i = block_size:block_size:min(200,m)
    plot([0.5, min(200,n)+0.5], [i+0.5, i+0.5], 'r-', 'LineWidth', 1);
    plot([i+0.5, i+0.5], [0.5, min(200,m)+0.5], 'r-', 'LineWidth', 1);
end

%% HELPER FUNCTIONS

function [M, block_data] = setup_block_jacobi(A, block_size)
    % Setup Block Jacobi preconditioner
    n = size(A, 1);
    num_blocks = floor(n / block_size);
    
    % Extract and factor diagonal blocks
    block_factors = cell(num_blocks, 1);
    successful = 0;
    failed = 0;
    
    for k = 1:num_blocks
        idx = (k-1)*block_size + (1:block_size);
        A_block = A(idx, idx);
        
        try
            % Try LU factorization with partial pivoting
            [L, U, P] = lu(A_block);
            % Check if factorization is reasonable
            if rcond(U) > 1e-14
                block_factors{k} = struct('L', L, 'U', U, 'P', P, 'method', 'lu');
                successful = successful + 1;
            else
                error('Poor conditioning in LU');
            end
        catch
            try
                % For saddle-point blocks, try regularized solve
                % Add small regularization to diagonal
                A_reg = A_block + 1e-12 * speye(size(A_block));
                [L, U, P] = lu(A_reg);
                block_factors{k} = struct('L', L, 'U', U, 'P', P, 'method', 'lu_reg');
                successful = successful + 1;
            catch
                try
                    % Fallback: QR for rank-deficient blocks
                    [Q, R, P] = qr(A_block);
                    block_factors{k} = struct('Q', Q, 'R', R, 'P', P, 'method', 'qr');
                    successful = successful + 1;
                catch
                    % Final fallback: pseudo-inverse
                    block_factors{k} = struct('Ainv', pinv(full(A_block)), 'method', 'pinv');
                    failed = failed + 1;
                end
            end
        end
    end
    
    % Create preconditioner function
    M = @(x) apply_block_jacobi(x, block_factors, block_size);
    
    block_data = struct('successful', successful, 'failed', failed, ...
                       'num_blocks', num_blocks, 'block_size', block_size);
end

function y = apply_block_jacobi(x, block_factors, block_size)
    % Apply Block Jacobi preconditioner
    n = length(x);
    num_blocks = n / block_size;
    y = zeros(size(x));
    
    for k = 1:num_blocks
        idx = (k-1)*block_size + (1:block_size);
        x_block = x(idx);
        
        switch block_factors{k}.method
            case {'lu', 'lu_reg'}
                y(idx) = block_factors{k}.U \ (block_factors{k}.L \ (block_factors{k}.P * x_block));
            case 'qr'
                % For QR: solve R*y = Q'*P*x
                y(idx) = block_factors{k}.P \ (block_factors{k}.R \ (block_factors{k}.Q' * x_block));
            case 'chol'
                y(idx) = block_factors{k}.R \ (block_factors{k}.R' \ x_block);
            case 'pinv'
                y(idx) = block_factors{k}.Ainv * x_block;
        end
    end
end