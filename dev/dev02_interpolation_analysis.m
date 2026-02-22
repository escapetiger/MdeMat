function dev02_interpolation_analysis()
% DEV02_INTERPOLATION_ANALYSIS Investigate RBF interpolation error on unit circle.
%
%   Analyzes the interpolation error distribution for unified RBF interpolation
%   in velocity space to understand ray effects and interpolation quality.

    clc; close all;

    fprintf('=== RBF Interpolation Error Analysis ===\n\n');

    % Setup test parameters
    params = setupTestParameters();

    % Create test functions with different characteristics
    test_functions = createTestFunctions();

    % Setup velocity grids and RBF centers
    [v_grids, rbf_centers] = setupVelocityGrids(params);

    % Analyze interpolation for each test function
    for i = 1:length(test_functions)
        fprintf('Analyzing test function %d: %s\n', i, test_functions{i}.name);
        analyzeInterpolationError(test_functions{i}, v_grids, rbf_centers, params, i);
    end

    % Investigate ray effect patterns
    investigateRayEffects(v_grids, rbf_centers, params);

    fprintf('\nAnalysis complete. Check figures for results.\n');

end

function params = setupTestParameters()
% SETUPTESTPARAMETERS Initialize parameters for interpolation analysis.

    % RBF parameters
    params.Nc = 32;           % Number of RBF centers
    params.epsilon = 1.0;     % RBF shape parameter
    params.Nv_coarse = 8;     % Coarse grid (discrete ordinate)
    params.Nv_fine = 64;      % Fine grid (for error analysis)

    % Analysis parameters
    params.poly_degree = 1;   % Linear polynomial reproduction
    params.plot_results = true;

end

function test_functions = createTestFunctions()
% CREATETESTFUNCTIONS Create various test functions to analyze interpolation.

    test_functions = {};

    % 1. Smooth isotropic function
    test_functions{1}.name = 'Smooth Isotropic';
    test_functions{1}.func = @(theta) 1 + 0.5 * cos(2*theta);
    test_functions{1}.description = 'Smooth function with low frequency variation';

    % 2. Sharp directional peak
    test_functions{2}.name = 'Sharp Directional Peak';
    test_functions{2}.func = @(theta) exp(-10*(theta - pi/4).^2);
    test_functions{2}.description = 'Sharp peak in specific direction (ray-like)';

    % 3. High frequency oscillation
    test_functions{3}.name = 'High Frequency';
    test_functions{3}.func = @(theta) 1 + 0.3 * sin(8*theta);
    test_functions{3}.description = 'High frequency angular variation';

    % 4. Discontinuous step function
    test_functions{4}.name = 'Step Function';
    test_functions{4}.func = @(theta) double(theta > pi & theta < 3*pi/2);
    test_functions{4}.description = 'Discontinuous angular distribution';

    % 5. Ray effect simulation
    test_functions{5}.name = 'Simulated Ray Effect';
    test_functions{5}.func = @(theta) createRayPattern(theta);
    test_functions{5}.description = 'Simulated discrete ordinate ray pattern';

end

function f = createRayPattern(theta)
% CREATERAYPATTERN Create a function mimicking ray effects from discrete ordinates.

    % Simulate 8 discrete rays with Gaussian spreading
    ray_angles = linspace(0, 2*pi, 9);
    ray_angles = ray_angles(1:end-1); % Remove duplicate

    f = zeros(size(theta));
    for i = 1:length(ray_angles)
        % Angular distance with periodic boundary
        angular_dist = min(abs(theta - ray_angles(i)), 2*pi - abs(theta - ray_angles(i)));
        f = f + exp(-angular_dist.^2 / (2 * (pi/16)^2));
    end

end

function [v_grids, rbf_centers] = setupVelocityGrids(params)
% SETUPVELOCITYGRIDS Create velocity grids and RBF centers.

    % Coarse grid (discrete ordinate resolution)
    theta_coarse = linspace(0, 2*pi, params.Nv_coarse+1);
    theta_coarse = theta_coarse(1:end-1);
    v_grids.coarse = [cos(theta_coarse)', sin(theta_coarse)'];
    v_grids.theta_coarse = theta_coarse;

    % Fine grid (for error analysis)
    theta_fine = linspace(0, 2*pi, params.Nv_fine+1);
    theta_fine = theta_fine(1:end-1);
    v_grids.fine = [cos(theta_fine)', sin(theta_fine)'];
    v_grids.theta_fine = theta_fine;

    % RBF centers
    theta_rbf = linspace(0, 2*pi, params.Nc+1);
    theta_rbf = theta_rbf(1:end-1);
    rbf_centers = [cos(theta_rbf)', sin(theta_rbf)'];

end

function analyzeInterpolationError(test_func, v_grids, rbf_centers, params, func_idx)
% ANALYZEINTERPOLATIONERROR Analyze interpolation error for a specific test function.

    % Evaluate test function on grids
    f_coarse = test_func.func(v_grids.theta_coarse);
    f_fine_exact = test_func.func(v_grids.theta_fine);

    % Apply unified RBF interpolation
    [f_fine_rbf, ~] = unifiedInterpolationAnalysis(v_grids.fine, rbf_centers, f_coarse, ...
                                                   v_grids.coarse, params.epsilon, params.poly_degree);

    % Compute interpolation error
    error_abs = abs(f_fine_rbf - f_fine_exact);
    error_rel = error_abs ./ (abs(f_fine_exact) + 1e-12);

    % Compute error metrics
    l2_error = sqrt(mean(error_abs.^2));
    max_error = max(error_abs);
    mean_rel_error = mean(error_rel);

    fprintf('  L2 Error: %.6f\n', l2_error);
    fprintf('  Max Error: %.6f\n', max_error);
    fprintf('  Mean Relative Error: %.6f\n', mean_rel_error);

    % Plot results
    if params.plot_results
        plotInterpolationResults(test_func, v_grids, f_coarse, f_fine_exact, f_fine_rbf, ...
                                error_abs, error_rel, func_idx);
    end

    % Analyze error patterns
    analyzeErrorPatterns(error_abs, v_grids.theta_fine, test_func.name);

end

function [s_eval, rho_analytic] = unifiedInterpolationAnalysis(eval_points, centers, values, data_points, epsilon, poly_degree)
% UNIFIEDINTERPOLATIONANALYSIS Unified interpolation for analysis (simplified).

    N_data = length(values);
    N_centers = size(centers, 1);
    M = size(eval_points, 1);

    % Build polynomial basis matrices
    P_data = buildPolynomialBasisAnalysis(data_points, poly_degree);
    P_centers = buildPolynomialBasisAnalysis(centers, poly_degree);
    P_eval = buildPolynomialBasisAnalysis(eval_points, poly_degree);
    Q = size(P_data, 2);

    % Build kernel matrices
    Phi_data_centers = zeros(N_data, N_centers);
    for i = 1:N_data
        for j = 1:N_centers
            r = norm(data_points(i,:) - centers(j,:));
            Phi_data_centers(i,j) = wendlandC2Analysis(r, epsilon);
        end
    end

    Phi_eval = zeros(M, N_centers);
    for i = 1:M
        for j = 1:N_centers
            r = norm(eval_points(i,:) - centers(j,:));
            Phi_eval(i,j) = wendlandC2Analysis(r, epsilon);
        end
    end

    % Solve unified system
    A = [Phi_data_centers, P_data; P_centers', zeros(Q, Q)];
    b = [values'; zeros(Q, 1)];
    coeffs = A \ b;
    c = coeffs(1:N_centers);
    d = coeffs(N_centers+1:N_centers+Q);

    % Evaluate interpolant
    s_eval = Phi_eval * c + P_eval * d;
    rho_analytic = sum(c) * (2*pi/(epsilon^2)) * (7/15); % Simplified analytical density

end

function P = buildPolynomialBasisAnalysis(points, degree)
% BUILDPOLYNOMIALBASISANALYSIS Build polynomial basis matrix.

    if size(points, 2) == 1
        % Convert angle to Cartesian
        x = cos(points);
        y = sin(points);
    else
        x = points(:, 1);
        y = points(:, 2);
    end

    N = length(x);

    if degree == 0
        P = ones(N, 1);
    elseif degree == 1
        P = [ones(N, 1), x, y];
    elseif degree == 2
        P = [ones(N, 1), x, y, x.^2, x.*y, y.^2];
    else
        error('Polynomial degree > 2 not implemented');
    end

end

function phi = wendlandC2Analysis(r, epsilon)
% WENDLANDC2ANALYSIS Wendland C2 compactly supported RBF.

    s = epsilon * r;
    if s >= 1
        phi = 0;
    else
        phi = (1 - s)^4 * (4*s + 1);
    end

end

function plotInterpolationResults(test_func, v_grids, f_coarse, f_fine_exact, f_fine_rbf, error_abs, error_rel, func_idx)
% PLOTINTERPOLATIONRESULTS Plot interpolation results and errors.

    figure(func_idx);

    subplot(2, 2, 1);
    plot(v_grids.theta_fine, f_fine_exact, 'k-', 'LineWidth', 2); hold on;
    plot(v_grids.theta_coarse, f_coarse, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    plot(v_grids.theta_fine, f_fine_rbf, 'b--', 'LineWidth', 1.5);
    xlabel('Angle (rad)'); ylabel('Function Value');
    title(sprintf('%s: Function Values', test_func.name));
    legend('Exact (Fine)', 'Data (Coarse)', 'RBF Interpolation', 'Location', 'best');
    grid on;

    subplot(2, 2, 2);
    semilogy(v_grids.theta_fine, error_abs, 'r-', 'LineWidth', 2);
    xlabel('Angle (rad)'); ylabel('Absolute Error (log scale)');
    title('Absolute Interpolation Error');
    grid on;

    subplot(2, 2, 3);
    plot(v_grids.theta_fine, error_rel, 'g-', 'LineWidth', 2);
    xlabel('Angle (rad)'); ylabel('Relative Error');
    title('Relative Interpolation Error');
    grid on;

    subplot(2, 2, 4);
    polarplot(v_grids.theta_fine, error_abs, 'r-', 'LineWidth', 2);
    title('Error Distribution (Polar)');

end

function analyzeErrorPatterns(error_abs, theta, func_name)
% ANALYZEERRORPATTERNS Analyze spatial patterns in interpolation error.

    % Find error peaks
    [peaks, peak_locs] = findpeaks(error_abs);
    peak_angles = theta(peak_locs);

    % Analyze periodicity
    fft_error = abs(fft(error_abs));
    frequencies = (0:length(fft_error)-1) / length(fft_error) * 2*pi;

    fprintf('  Error Analysis for %s:\n', func_name);
    fprintf('    Number of error peaks: %d\n', length(peaks));
    if ~isempty(peaks)
        fprintf('    Largest error peak: %.6f at angle %.3f\n', max(peaks), peak_angles(peaks == max(peaks)));
    end

    % Check for dominant frequencies (ray-like patterns)
    [~, dominant_freq_idx] = max(fft_error(2:end/2)); % Exclude DC component
    dominant_freq = frequencies(dominant_freq_idx + 1);
    fprintf('    Dominant error frequency: %.3f (corresponding to %d-fold symmetry)\n', ...
            dominant_freq, round(2*pi/dominant_freq));

end

function investigateRayEffects(v_grids, rbf_centers, params)
% INVESTIGATERAYEFFECTS Investigate ray effect patterns specifically.

    fprintf('\n=== Ray Effect Investigation ===\n');

    % Create test scenario: 8 discrete rays (typical DO method)
    ray_angles = linspace(0, 2*pi, 9);
    ray_angles = ray_angles(1:end-1);

    % Test different RBF center densities
    center_counts = [8, 16, 32, 64];

    figure(10);

    for i = 1:length(center_counts)
        Nc = center_counts(i);

        % Setup RBF centers
        theta_rbf = linspace(0, 2*pi, Nc+1);
        theta_rbf = theta_rbf(1:end-1);
        centers = [cos(theta_rbf)', sin(theta_rbf)'];

        % Create ray pattern
        f_rays = zeros(size(ray_angles));
        for j = 1:length(ray_angles)
            f_rays(j) = 1.0; % Unit intensity per ray
        end

        % Interpolate to fine grid
        ray_data_points = [cos(ray_angles)', sin(ray_angles)'];
        [f_interpolated, ~] = unifiedInterpolationAnalysis(v_grids.fine, centers, f_rays, ...
                                                          ray_data_points, params.epsilon, params.poly_degree);

        % Plot results
        subplot(2, 2, i);
        plot(v_grids.theta_fine, f_interpolated, 'b-', 'LineWidth', 2); hold on;
        stem(ray_angles, f_rays, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        xlabel('Angle (rad)'); ylabel('Interpolated Value');
        title(sprintf('Ray Effect Mitigation (%d RBF centers)', Nc));
        legend('RBF Interpolation', 'Original Rays', 'Location', 'best');
        grid on;

        % Analyze smoothing quality
        variation = std(f_interpolated);
        fprintf('  %d RBF centers: std deviation = %.6f\n', Nc, variation);
    end

    % Compare with ideal smooth distribution
    f_ideal = ones(size(v_grids.theta_fine)) * mean(f_rays);
    fprintf('  Ideal smooth distribution: std deviation = %.6f\n', std(f_ideal));

end