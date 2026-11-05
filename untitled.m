clear; clc; close all;

%% Comparison of CardinalKernelBasis vs RbfKernelBasis
fprintf('Comparing CardinalKernelBasis and RbfKernelBasis\n');
fprintf('================================================\n\n');

%% Test function definitions
f1 = @(x) exp(-2*sum(x.^2, 1));
f2 = @(x) sin(3*pi*x(1,:)) .* cos(2*pi*x(2,:));
f3 = @(x) 1./(1 + 5*sum((x - 0.3).^2, 1));

functions = {f1, f2, f3};
function_names = {'Gaussian Peak', 'Oscillatory', 'Rational Peak'};

%% Test configurations
test_configs = struct();
test_configs.nNodes = [9, 16, 25, 36];
test_configs.nEval = 50;
test_configs.epsilon_values = [0.5, 1.0, 2.0, 5.0];
test_configs.rbf_kernels = {'gaussian', 'multiquadric', 'inverse_multiquadric'};

%% Setup evaluation domain
[x_eval, y_eval] = meshgrid(linspace(0, 1, test_configs.nEval), ...
                            linspace(0, 1, test_configs.nEval));
X_eval = [x_eval(:)'; y_eval(:)'];

%% Main comparison loop
for f_idx = 1:length(functions)
    fprintf('Testing function: %s\n', function_names{f_idx});
    fprintf('------------------\n');
    
    f = functions{f_idx};
    f_true = f(X_eval);
    
    for n_idx = 1:length(test_configs.nNodes)
        nNodes = test_configs.nNodes(n_idx);
        nSide = sqrt(nNodes);
        
        [x_nodes, y_nodes] = meshgrid(linspace(0, 1, nSide), ...
                                      linspace(0, 1, nSide));
        X_nodes = [x_nodes(:)'; y_nodes(:)'];
        f_nodes = f(X_nodes);
        
        fprintf('  Nodes: %d (%dx%d grid)\n', nNodes, nSide, nSide);
        
        try
            cardinal_basis = approx.basis.CardinalKernelBasis(X_nodes);
            Y_cardinal = cardinal_basis.evaluate(X_eval);
            
            f_cardinal = Y_cardinal' * f_nodes(:);
            err_cardinal = abs(f_true(:) - f_cardinal);
            max_err_cardinal = max(err_cardinal);
            rms_err_cardinal = sqrt(mean(err_cardinal.^2));
            
            fprintf('    Cardinal Basis:\n');
            fprintf('      Max Error: %.3e\n', max_err_cardinal);
            fprintf('      RMS Error: %.3e\n', rms_err_cardinal);
            
        catch ME
            fprintf('    Cardinal Basis: Failed (%s)\n', ME.message);
            max_err_cardinal = Inf;
            rms_err_cardinal = Inf;
        end
        
        best_rbf_max = Inf;
        best_rbf_rms = Inf;
        best_config = '';
        
        for kernel_idx = 1:length(test_configs.rbf_kernels)
            kernel_type = test_configs.rbf_kernels{kernel_idx};
            
            for eps_idx = 1:length(test_configs.epsilon_values)
                epsilon = test_configs.epsilon_values(eps_idx);
                
                try
                    rbf_basis = approx.basis.RbfKernelBasis(X_nodes, kernel_type, epsilon);
                    Y_rbf = rbf_basis.evaluate(X_eval);
                    
                    f_rbf = Y_rbf' * f_nodes(:);
                    err_rbf = abs(f_true(:) - f_rbf);
                    max_err_rbf = max(err_rbf);
                    rms_err_rbf = sqrt(mean(err_rbf.^2));
                    
                    if max_err_rbf < best_rbf_max
                        best_rbf_max = max_err_rbf;
                        best_rbf_rms = rms_err_rbf;
                        best_config = sprintf('%s, ε=%.1f', kernel_type, epsilon);
                    end
                    
                catch ME
                    continue;
                end
            end
        end
        
        if ~isinf(best_rbf_max)
            fprintf('    Best RBF (%s):\n', best_config);
            fprintf('      Max Error: %.3e\n', best_rbf_max);
            fprintf('      RMS Error: %.3e\n', best_rbf_rms);
            
            if best_rbf_max < max_err_cardinal
                fprintf('      --> RBF wins (%.1fx better)\n', ...
                        max_err_cardinal / best_rbf_max);
            else
                fprintf('      --> Cardinal wins (%.1fx better)\n', ...
                        best_rbf_max / max_err_cardinal);
            end
        else
            fprintf('    Best RBF: All failed\n');
        end
        
        fprintf('\n');
    end
    
    fprintf('\n');
end

%% Detailed comparison for visualization
fprintf('Detailed Comparison for Visualization\n');
fprintf('====================================\n');

nNodes = 50;
nSide = 5;
[x_nodes, y_nodes] = meshgrid(linspace(0, 1, nSide), linspace(0, 1, nSide));
X_nodes = [x_nodes(:)'; y_nodes(:)'];

f_test = f2;
f_nodes = f_test(X_nodes);
f_true_grid = reshape(f_test(X_eval), test_configs.nEval, test_configs.nEval);

try
    cardinal_basis = approx.basis.CardinalKernelBasis(X_nodes);
    Y_cardinal = cardinal_basis.evaluate(X_eval);
    f_cardinal = Y_cardinal' * f_nodes(:);
    f_cardinal_grid = reshape(f_cardinal, test_configs.nEval, test_configs.nEval);
    cardinal_success = true;
catch ME
    fprintf('Cardinal basis failed: %s\n', ME.message);
    f_cardinal_grid = nan(test_configs.nEval, test_configs.nEval);
    cardinal_success = false;
end

rbf_basis = approx.basis.RbfKernelBasis(X_nodes, 'gaussian', 2.0);
Y_rbf = rbf_basis.evaluate(X_eval);
f_rbf = f_nodes(:).' * Y_rbf;
f_rbf_grid = reshape(f_rbf, test_configs.nEval, test_configs.nEval);

figure('Position', [100, 100, 1200, 800]);

subplot(2, 3, 1);
surf(x_eval, y_eval, f_true_grid);
title('True Function: sin(3πx)cos(2πy)');
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
shading interp; colorbar;

subplot(2, 3, 2);
if cardinal_success
    surf(x_eval, y_eval, f_cardinal_grid);
    title('Cardinal Basis Approximation');
else
    text(0.5, 0.5, 'Cardinal Basis\nFailed', 'HorizontalAlignment', 'center');
    title('Cardinal Basis (Failed)');
end
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
shading interp; colorbar;

subplot(2, 3, 3);
surf(x_eval, y_eval, f_rbf_grid);
title('RBF Approximation (Gaussian, ε=2)');
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
shading interp; colorbar;

subplot(2, 3, 4);
if cardinal_success
    contour(x_eval, y_eval, abs(f_true_grid - f_cardinal_grid), 20);
    title('Cardinal Error');
else
    text(0.5, 0.5, 'No Error\n(Basis Failed)', 'HorizontalAlignment', 'center');
    title('Cardinal Error (N/A)');
end
xlabel('x'); ylabel('y');
colorbar;
hold on;
plot(X_nodes(1,:), X_nodes(2,:), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

subplot(2, 3, 5);
contour(x_eval, y_eval, abs(f_true_grid - f_rbf_grid), 20);
title('RBF Error');
xlabel('x'); ylabel('y');
colorbar;
hold on;
plot(X_nodes(1,:), X_nodes(2,:), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

subplot(2, 3, 6);
hold on;
nEval_1d = 100;
x_slice = linspace(0, 1, nEval_1d);
X_slice = [x_slice; 0.5 * ones(1, nEval_1d)];

f_true_slice = f_test(X_slice);
plot(x_slice, f_true_slice, 'k-', 'LineWidth', 2, 'DisplayName', 'True');

if cardinal_success
    Y_cardinal_slice = cardinal_basis.evaluate(X_slice);
    f_cardinal_slice = Y_cardinal_slice' * f_nodes(:);
    plot(x_slice, f_cardinal_slice, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Cardinal');
end

Y_rbf_slice = rbf_basis.evaluate(X_slice);
f_rbf_slice = Y_rbf_slice' * f_nodes(:);
plot(x_slice, f_rbf_slice, 'r:', 'LineWidth', 1.5, 'DisplayName', 'RBF');

x_nodes_slice = X_nodes(1, X_nodes(2,:) == 0.2);
f_nodes_slice = f_test([x_nodes_slice; 0.5 * ones(size(x_nodes_slice))]);
plot(x_nodes_slice, f_nodes_slice, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Nodes');

title('Cross-section at y=0.5');
xlabel('x'); ylabel('f(x,0.5)');
legend('Location', 'best');
grid on;

sgtitle('Comparison: Cardinal vs RBF Kernel Basis', 'FontSize', 14, 'FontWeight', 'bold');

%% Performance comparison
fprintf('Performance Comparison\n');
fprintf('=====================\n');

nNodes_perf = 100;
X_nodes_perf = rand(2, nNodes_perf);
X_eval_perf = rand(2, 1000);

try
    tic;
    cardinal_basis_perf = approx.basis.CardinalKernelBasis(X_nodes_perf);
    Y_cardinal_perf = cardinal_basis_perf.evaluate(X_eval_perf);
    time_cardinal = toc;
    fprintf('Cardinal Basis (100 nodes, 1000 eval): %.3f seconds\n', time_cardinal);
catch ME
    fprintf('Cardinal Basis: Failed (%s)\n', ME.message);
    time_cardinal = Inf;
end

tic;
rbf_basis_perf = approx.basis.RbfKernelBasis(X_nodes_perf, 'gaussian', 1.0);
Y_rbf_perf = rbf_basis_perf.evaluate(X_eval_perf);
time_rbf = toc;
fprintf('RBF Basis (100 nodes, 1000 eval): %.3f seconds\n', time_rbf);

if ~isinf(time_cardinal)
    fprintf('Speed ratio (Cardinal/RBF): %.2f\n', time_cardinal / time_rbf);
end

%% Derivative comparison
fprintf('\nDerivative Comparison\n');
fprintf('====================\n');

X_test = [0.3; 0.7];

try
    cardinal_basis_deriv = approx.basis.CardinalKernelBasis(X_nodes);
    dY_cardinal_dx = cardinal_basis_deriv.differentiate(X_test, [1, 0]);
    dY_cardinal_dy = cardinal_basis_deriv.differentiate(X_test, [0, 1]);
    fprintf('Cardinal derivatives: available\n');
catch ME
    fprintf('Cardinal derivatives: not implemented (%s)\n', ME.message);
end

rbf_basis_deriv = approx.basis.RbfKernelBasis(X_nodes, 'gaussian', 1.0);
dY_rbf_dx = rbf_basis_deriv.differentiate(X_test, [1, 0]);
dY_rbf_dy = rbf_basis_deriv.differentiate(X_test, [0, 1]);
fprintf('RBF derivatives: available\n');

fprintf('\nComparison Summary:\n');
fprintf('==================\n');
fprintf('Cardinal Basis:\n');
fprintf('  + Exact interpolation at nodes\n');
fprintf('  + No shape parameter tuning\n');
fprintf('  + Partition of unity property\n');
fprintf('  - May be numerically unstable\n');
fprintf('  - Limited to specific node configurations\n');
fprintf('  - Derivatives not implemented\n\n');

fprintf('RBF Basis:\n');
fprintf('  + Numerically stable\n');
fprintf('  + Flexible for scattered data\n');
fprintf('  + Multiple kernel options\n');
fprintf('  + Derivatives available\n');
fprintf('  - Requires epsilon tuning\n');
fprintf('  - Not exact interpolation (without constraints)\n');