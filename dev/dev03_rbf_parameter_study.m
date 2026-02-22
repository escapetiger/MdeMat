function dev03_rbf_parameter_study()
% DEV03_RBF_PARAMETER_STUDY Study RBF interpolation parameter sensitivity.
%
%   Investigates how different parameters affect interpolation quality
%   and ray effect mitigation in unified RBF interpolation.

    clc; close all;

    fprintf('=== RBF Parameter Sensitivity Study ===\n\n');

    % Setup base parameters
    base_params = setupBaseParameters();

    % Study 1: Shape parameter (epsilon) sensitivity
    fprintf('Study 1: Shape parameter sensitivity\n');
    studyShapeParameter(base_params);

    % Study 2: Number of RBF centers
    fprintf('\nStudy 2: RBF center count sensitivity\n');
    studyRBFCenterCount(base_params);

    % Study 3: Polynomial degree effects
    fprintf('\nStudy 3: Polynomial degree effects\n');
    studyPolynomialDegree(base_params);

    % Study 4: Ray effect mitigation strategies
    fprintf('\nStudy 4: Ray effect mitigation comparison\n');
    compareRayMitigationStrategies(base_params);

    fprintf('\nParameter study complete. Check figures for results.\n');

end

function params = setupBaseParameters()
% SETUPBASEPARAMETERS Set up baseline parameters for studies.

    params.Nv_coarse = 8;     % Discrete ordinate resolution
    params.Nv_fine = 64;     % Analysis resolution
    params.Nc_base = 32;     % Base number of RBF centers
    params.epsilon_base = 1.0; % Base shape parameter
    params.poly_degree_base = 1; % Base polynomial degree

    % Test function (ray-like pattern)
    params.test_function = @(theta) createTestRayPattern(theta, 8);

end

function f = createTestRayPattern(theta, n_rays)
% CREATETESTRAYPATTERN Create test pattern with specified number of rays.

    ray_angles = linspace(0, 2*pi, n_rays+1);
    ray_angles = ray_angles(1:end-1);

    f = zeros(size(theta));
    for i = 1:length(ray_angles)
        angular_dist = min(abs(theta - ray_angles(i)), 2*pi - abs(theta - ray_angles(i)));
        f = f + exp(-angular_dist.^2 / (2 * (pi/32)^2));
    end

end

function studyShapeParameter(params)
% STUDYSHAPEPARAMETER Study effect of RBF shape parameter epsilon.

    epsilon_values = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0];
    errors = zeros(size(epsilon_values));
    smoothness = zeros(size(epsilon_values));

    % Setup grids
    theta_coarse = linspace(0, 2*pi, params.Nv_coarse+1);
    theta_coarse = theta_coarse(1:end-1);
    theta_fine = linspace(0, 2*pi, params.Nv_fine+1);
    theta_fine = theta_fine(1:end-1);

    v_coarse = [cos(theta_coarse)', sin(theta_coarse)'];
    v_fine = [cos(theta_fine)', sin(theta_fine)'];

    % RBF centers
    theta_rbf = linspace(0, 2*pi, params.Nc_base+1);
    theta_rbf = theta_rbf(1:end-1);
    rbf_centers = [cos(theta_rbf)', sin(theta_rbf)'];

    % Test function values
    f_coarse = params.test_function(theta_coarse);
    f_fine_exact = params.test_function(theta_fine);

    figure(20);

    for i = 1:length(epsilon_values)
        epsilon = epsilon_values(i);

        % Apply interpolation
        [f_interpolated, ~] = unifiedInterpolationParameterStudy(v_fine, rbf_centers, f_coarse, ...
                                                                v_coarse, epsilon, params.poly_degree_base);

        % Compute metrics
        errors(i) = sqrt(mean((f_interpolated - f_fine_exact).^2));
        smoothness(i) = computeSmoothness(f_interpolated);

        % Plot sample results
        if i <= 4
            subplot(2, 2, i);
            plot(theta_fine, f_fine_exact, 'k-', 'LineWidth', 2); hold on;
            plot(theta_coarse, f_coarse, 'ro', 'MarkerSize', 6);
            plot(theta_fine, f_interpolated, 'b--', 'LineWidth', 1.5);
            title(sprintf('\\epsilon = %.1f', epsilon));
            xlabel('Angle (rad)'); ylabel('Value');
            legend('Exact', 'Data', 'RBF', 'Location', 'best');
            grid on;
        end

        fprintf('  ε = %.1f: L2 error = %.6f, smoothness = %.6f\n', epsilon, errors(i), smoothness(i));
    end

    % Summary plot
    figure(21);
    subplot(1, 2, 1);
    semilogx(epsilon_values, errors, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Shape Parameter \\epsilon'); ylabel('L2 Error');
    title('Interpolation Error vs Shape Parameter');
    grid on;

    subplot(1, 2, 2);
    semilogx(epsilon_values, smoothness, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Shape Parameter \\epsilon'); ylabel('Smoothness Index');
    title('Solution Smoothness vs Shape Parameter');
    grid on;

end

function studyRBFCenterCount(params)
% STUDYRBFCENTERCOUNT Study effect of number of RBF centers.

    center_counts = [8, 16, 24, 32, 48, 64];
    errors = zeros(size(center_counts));
    condition_numbers = zeros(size(center_counts));

    % Setup grids
    theta_coarse = linspace(0, 2*pi, params.Nv_coarse+1);
    theta_coarse = theta_coarse(1:end-1);
    theta_fine = linspace(0, 2*pi, params.Nv_fine+1);
    theta_fine = theta_fine(1:end-1);

    v_coarse = [cos(theta_coarse)', sin(theta_coarse)'];
    v_fine = [cos(theta_fine)', sin(theta_fine)'];

    f_coarse = params.test_function(theta_coarse);
    f_fine_exact = params.test_function(theta_fine);

    figure(22);

    for i = 1:length(center_counts)
        Nc = center_counts(i);

        % Setup RBF centers
        theta_rbf = linspace(0, 2*pi, Nc+1);
        theta_rbf = theta_rbf(1:end-1);
        rbf_centers = [cos(theta_rbf)', sin(theta_rbf)'];

        % Apply interpolation and compute condition number
        [f_interpolated, cond_num] = unifiedInterpolationParameterStudy(v_fine, rbf_centers, f_coarse, ...
                                                                        v_coarse, params.epsilon_base, params.poly_degree_base);

        errors(i) = sqrt(mean((f_interpolated - f_fine_exact).^2));
        condition_numbers(i) = cond_num;

        % Plot sample results
        if i <= 4
            subplot(2, 2, i);
            plot(theta_fine, f_fine_exact, 'k-', 'LineWidth', 2); hold on;
            plot(theta_coarse, f_coarse, 'ro', 'MarkerSize', 6);
            plot(theta_fine, f_interpolated, 'b--', 'LineWidth', 1.5);
            title(sprintf('N_c = %d', Nc));
            xlabel('Angle (rad)'); ylabel('Value');
            legend('Exact', 'Data', 'RBF', 'Location', 'best');
            grid on;
        end

        fprintf('  Nc = %d: L2 error = %.6f, condition = %.2e\n', Nc, errors(i), condition_numbers(i));
    end

    % Summary plot
    figure(23);
    subplot(1, 2, 1);
    plot(center_counts, errors, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Number of RBF Centers'); ylabel('L2 Error');
    title('Interpolation Error vs RBF Center Count');
    grid on;

    subplot(1, 2, 2);
    semilogy(center_counts, condition_numbers, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Number of RBF Centers'); ylabel('Condition Number (log scale)');
    title('System Condition vs RBF Center Count');
    grid on;

end

function studyPolynomialDegree(params)
% STUDYPOLYNOMIALDEGREE Study effect of polynomial reproduction degree.

    poly_degrees = [0, 1, 2];
    degree_names = {'Constant', 'Linear', 'Quadratic'};

    % Setup grids
    theta_coarse = linspace(0, 2*pi, params.Nv_coarse+1);
    theta_coarse = theta_coarse(1:end-1);
    theta_fine = linspace(0, 2*pi, params.Nv_fine+1);
    theta_fine = theta_fine(1:end-1);

    v_coarse = [cos(theta_coarse)', sin(theta_coarse)'];
    v_fine = [cos(theta_fine)', sin(theta_fine)'];

    % RBF centers
    theta_rbf = linspace(0, 2*pi, params.Nc_base+1);
    theta_rbf = theta_rbf(1:end-1);
    rbf_centers = [cos(theta_rbf)', sin(theta_rbf)'];

    f_coarse = params.test_function(theta_coarse);
    f_fine_exact = params.test_function(theta_fine);

    figure(24);

    for i = 1:length(poly_degrees)
        degree = poly_degrees(i);

        try
            [f_interpolated, ~] = unifiedInterpolationParameterStudy(v_fine, rbf_centers, f_coarse, ...
                                                                    v_coarse, params.epsilon_base, degree);

            error = sqrt(mean((f_interpolated - f_fine_exact).^2));
            smoothness = computeSmoothness(f_interpolated);

            subplot(1, 3, i);
            plot(theta_fine, f_fine_exact, 'k-', 'LineWidth', 2); hold on;
            plot(theta_coarse, f_coarse, 'ro', 'MarkerSize', 6);
            plot(theta_fine, f_interpolated, 'b--', 'LineWidth', 1.5);
            title(sprintf('%s (degree %d)', degree_names{i}, degree));
            xlabel('Angle (rad)'); ylabel('Value');
            legend('Exact', 'Data', 'RBF', 'Location', 'best');
            grid on;

            fprintf('  Degree %d (%s): L2 error = %.6f, smoothness = %.6f\n', ...
                    degree, degree_names{i}, error, smoothness);

        catch ME
            fprintf('  Degree %d failed: %s\n', degree, ME.message);
        end
    end

end

function compareRayMitigationStrategies(params)
% COMPARERAYMITIGATIONSTRATEGIES Compare different approaches to ray effect mitigation.

    strategies = {'None', 'Wider Kernels', 'More Centers', 'Higher Polynomial', 'Combined'};

    % Setup grids
    theta_coarse = linspace(0, 2*pi, params.Nv_coarse+1);
    theta_coarse = theta_coarse(1:end-1);
    theta_fine = linspace(0, 2*pi, params.Nv_fine+1);
    theta_fine = theta_fine(1:end-1);

    v_coarse = [cos(theta_coarse)', sin(theta_coarse)'];
    v_fine = [cos(theta_fine)', sin(theta_fine)'];

    f_coarse = params.test_function(theta_coarse);
    f_fine_exact = ones(size(theta_fine)) * mean(f_coarse); % Target: smooth isotropic

    figure(25);

    for i = 1:length(strategies)
        strategy = strategies{i};

        % Configure parameters for each strategy
        switch strategy
            case 'None'
                epsilon = params.epsilon_base;
                Nc = params.Nc_base;
                poly_degree = 0;
            case 'Wider Kernels'
                epsilon = params.epsilon_base * 0.5; % Wider support
                Nc = params.Nc_base;
                poly_degree = params.poly_degree_base;
            case 'More Centers'
                epsilon = params.epsilon_base;
                Nc = params.Nc_base * 2; % More centers
                poly_degree = params.poly_degree_base;
            case 'Higher Polynomial'
                epsilon = params.epsilon_base;
                Nc = params.Nc_base;
                poly_degree = 2; % Quadratic reproduction
            case 'Combined'
                epsilon = params.epsilon_base * 0.5;
                Nc = params.Nc_base * 1.5;
                poly_degree = 1;
        end

        % Setup RBF centers
        theta_rbf = linspace(0, 2*pi, Nc+1);
        theta_rbf = theta_rbf(1:end-1);
        rbf_centers = [cos(theta_rbf)', sin(theta_rbf)'];

        try
            [f_interpolated, ~] = unifiedInterpolationParameterStudy(v_fine, rbf_centers, f_coarse, ...
                                                                    v_coarse, epsilon, poly_degree);

            % Compute ray effect metric (deviation from isotropy)
            ray_effect = std(f_interpolated) / mean(f_interpolated);
            error_to_smooth = sqrt(mean((f_interpolated - f_fine_exact).^2));

            subplot(2, 3, i);
            plot(theta_fine, f_fine_exact, 'k-', 'LineWidth', 2); hold on;
            plot(theta_coarse, f_coarse, 'ro', 'MarkerSize', 6);
            plot(theta_fine, f_interpolated, 'b--', 'LineWidth', 1.5);
            title(strategy);
            xlabel('Angle (rad)'); ylabel('Value');
            legend('Target Smooth', 'Ray Data', 'RBF Result', 'Location', 'best');
            grid on;

            fprintf('  %s: Ray effect = %.6f, Error to smooth = %.6f\n', ...
                    strategy, ray_effect, error_to_smooth);

        catch ME
            fprintf('  %s failed: %s\n', strategy, ME.message);
        end
    end

end

function [s_eval, condition_number] = unifiedInterpolationParameterStudy(eval_points, centers, values, data_points, epsilon, poly_degree)
% UNIFIEDINTERPOLATIONPARAMETERSTUDY Unified interpolation with condition number monitoring.

    N_data = length(values);
    N_centers = size(centers, 1);
    M = size(eval_points, 1);

    % Build polynomial basis matrices
    P_data = buildPolynomialBasisStudy(data_points, poly_degree);
    P_centers = buildPolynomialBasisStudy(centers, poly_degree);
    P_eval = buildPolynomialBasisStudy(eval_points, poly_degree);
    Q = size(P_data, 2);

    % Build kernel matrices
    Phi_data_centers = zeros(N_data, N_centers);
    for i = 1:N_data
        for j = 1:N_centers
            r = norm(data_points(i,:) - centers(j,:));
            Phi_data_centers(i,j) = wendlandC2Study(r, epsilon);
        end
    end

    Phi_eval = zeros(M, N_centers);
    for i = 1:M
        for j = 1:N_centers
            r = norm(eval_points(i,:) - centers(j,:));
            Phi_eval(i,j) = wendlandC2Study(r, epsilon);
        end
    end

    % Solve unified system and compute condition number
    A = [Phi_data_centers, P_data; P_centers', zeros(Q, Q)];
    condition_number = cond(A);

    b = [values'; zeros(Q, 1)];
    coeffs = A \ b;
    c = coeffs(1:N_centers);
    d = coeffs(N_centers+1:N_centers+Q);

    % Evaluate interpolant
    s_eval = Phi_eval * c + P_eval * d;

end

function P = buildPolynomialBasisStudy(points, degree)
% BUILDPOLYNOMIALBASISTUDY Build polynomial basis for parameter study.

    if size(points, 2) == 1
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

function phi = wendlandC2Study(r, epsilon)
% WENDLANDC2STUDY Wendland C2 RBF for parameter study.

    s = epsilon * r;
    if s >= 1
        phi = 0;
    else
        phi = (1 - s)^4 * (4*s + 1);
    end

end

function smoothness = computeSmoothness(f)
% COMPUTESMOOTHNESS Compute smoothness metric (inverse of total variation).

    % Compute discrete derivative
    df = diff(f);

    % Handle periodic boundary
    df = [df; f(1) - f(end)];

    % Total variation
    total_variation = sum(abs(df));

    % Smoothness index (lower total variation = higher smoothness)
    smoothness = 1 / (1 + total_variation);

end