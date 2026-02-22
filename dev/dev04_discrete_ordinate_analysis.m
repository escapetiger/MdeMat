function dev04_discrete_ordinate_analysis()
% DEV04_DISCRETE_ORDINATE_ANALYSIS Analyze failure modes of discrete ordinate method.
%
%   Investigates why the classical discrete ordinate method fails and
%   demonstrates the fundamental sources of ray effects and angular discretization errors.

    clc; close all;

    fprintf('=== Discrete Ordinate Method Failure Analysis ===\n\n');

    % Setup analysis parameters
    params = setupAnalysisParameters();

    % Analysis 1: Angular discretization effects
    fprintf('Analysis 1: Angular discretization effects\n');
    analyzeAngularDiscretization(params);

    % Analysis 2: Ray effect demonstration
    fprintf('\nAnalysis 2: Ray effect demonstration\n');
    demonstrateRayEffects(params);

    % Analysis 3: Information propagation analysis
    fprintf('\nAnalysis 3: Information propagation patterns\n');
    analyzeInformationPropagation(params);

    % Analysis 4: Comparison with exact solution characteristics
    fprintf('\nAnalysis 4: Comparison with exact characteristics\n');
    compareWithExactCharacteristics(params);

    % Analysis 5: Angular coupling investigation
    fprintf('\nAnalysis 5: Angular coupling investigation\n');
    investigateAngularCoupling(params);

    fprintf('\nDiscrete ordinate analysis complete. Check figures for detailed results.\n');

end

function params = setupAnalysisParameters()
% SETUPANALYSISPARAMETERS Setup parameters for discrete ordinate analysis.

    % Spatial domain
    params.x_min = -0.5; params.x_max = 0.5;
    params.y_min = -0.5; params.y_max = 0.5;
    params.Nx = 64; params.Ny = 64;

    % Time parameters
    params.T_final = 0.2;
    params.Nt = 20;
    params.dt = params.T_final / params.Nt;

    % Angular discretizations to compare
    params.Nv_values = [4, 8, 16, 32];

    % Initial condition parameters
    params.sigma = 0.02;  % Gaussian width
    params.x0 = 0; params.y0 = 0;  % Initial center

    % Analysis parameters
    params.output_times = [0.05, 0.1, 0.15, 0.2];

end

function analyzeAngularDiscretization(params)
% ANALYZEANGULARDISCRETIZATION Analyze effects of angular discretization.

    fprintf('  Investigating angular resolution effects...\n');

    % Setup spatial grid
    x_grid = linspace(params.x_min, params.x_max, params.Nx);
    y_grid = linspace(params.y_min, params.y_max, params.Ny);
    [X, Y] = meshgrid(x_grid, y_grid);

    figure(1);

    for i = 1:length(params.Nv_values)
        Nv = params.Nv_values(i);

        fprintf('    Nv = %d directions\n', Nv);

        % Setup velocity grid
        theta = linspace(0, 2*pi, Nv+1);
        theta = theta(1:end-1);
        v_grid = [cos(theta)', sin(theta)'];

        % Setup initial condition
        f0 = setupInitialCondition(X, Y, params, Nv);

        % Time array
        times = linspace(0, params.T_final, params.Nt+1);

        % Solve with discrete ordinate method
        f_all = solveDiscreteOrdinateAnalysis(f0, params, x_grid, y_grid, v_grid, times);

        % Compute mass density at final time
        rho_final = sum(f_all(:,:,:,end), 3) * (2*pi/Nv);

        % Plot result
        subplot(2, 2, i);
        imagesc(x_grid, y_grid, rho_final, [0, max(rho_final(:))]);
        colorbar;
        title(sprintf('N_v = %d', Nv));
        xlabel('x'); ylabel('y');
        axis equal tight;

        % Analyze angular errors
        analyzeAngularErrors(rho_final, x_grid, y_grid, Nv, params);
    end

    sgtitle('Angular Discretization Effects on Mass Density');

end

function analyzeAngularErrors(rho, x_grid, y_grid, Nv, params)
% ANALYZEANGULARERRORS Analyze specific angular discretization errors.

    % Find center of mass
    [X, Y] = meshgrid(x_grid, y_grid);
    total_mass = sum(rho(:)) * (x_grid(2)-x_grid(1)) * (y_grid(2)-y_grid(1));
    x_cm = sum(sum(rho .* X)) / sum(rho(:));
    y_cm = sum(sum(rho .* Y)) / sum(rho(:));

    % Analyze symmetry breaking
    % Check if the solution maintains rotational symmetry
    r = sqrt((X - x_cm).^2 + (Y - y_cm).^2);
    theta_grid = atan2(Y - y_cm, X - x_cm);

    % Compute angular average and deviations
    theta_bins = linspace(-pi, pi, Nv*4);
    r_max = min(abs([params.x_min, params.x_max, params.y_min, params.y_max])) * 0.8;

    % Sample at fixed radius
    if r_max > 0.1
        mask = (r > r_max*0.8) & (r < r_max*1.2);
        if sum(mask(:)) > 10
            rho_angular = rho(mask);
            theta_angular = theta_grid(mask);

            angular_variation = std(rho_angular) / mean(rho_angular);
            fprintf('      Angular variation (std/mean) at r=%.2f: %.6f\n', r_max, angular_variation);
        end
    end

    % Count discrete rays in the solution
    countDiscreteRays(rho, x_grid, y_grid, Nv);

end

function countDiscreteRays(rho, x_grid, y_grid, Nv)
% COUNTDISCRETERAYS Count and identify discrete ray structures.

    [Ny, Nx] = size(rho);
    center_i = round(Ny/2);
    center_j = round(Nx/2);

    % Sample along radial lines from center
    n_angles = 64;
    angles = linspace(0, 2*pi, n_angles+1);
    angles = angles(1:end-1);

    radius = min(center_i, center_j) * 0.8;

    radial_profile = zeros(size(angles));
    for k = 1:length(angles)
        % Compute end point of radial line
        end_i = center_i + radius * sin(angles(k));
        end_j = center_j + radius * cos(angles(k));

        % Interpolate along line
        if end_i >= 1 && end_i <= Ny && end_j >= 1 && end_j <= Nx
            radial_profile(k) = interp2(1:Nx, 1:Ny, rho, end_j, end_i, 'linear', 0);
        end
    end

    % Find peaks in radial profile (discrete rays)
    [peaks, peak_locs] = findpeaks(radial_profile, 'MinPeakHeight', max(radial_profile)*0.1);

    fprintf('      Detected %d ray-like peaks (expected %d)\n', length(peaks), Nv);

    if length(peaks) > 0
        peak_angles = angles(peak_locs);
        angular_spacing = diff(sort(peak_angles));
        if ~isempty(angular_spacing)
            fprintf('      Mean angular spacing: %.3f rad (expected %.3f)\n', ...
                    mean(angular_spacing), 2*pi/Nv);
        end
    end

end

function demonstrateRayEffects(params)
% DEMONSTRATERAYEFFECTS Demonstrate ray effects with detailed visualization.

    fprintf('  Demonstrating ray effect mechanisms...\n');

    % Use 8 directions for clear ray effect demonstration
    Nv = 8;
    theta = linspace(0, 2*pi, Nv+1);
    theta = theta(1:end-1);
    v_grid = [cos(theta)', sin(theta)'];

    % Setup grids
    x_grid = linspace(params.x_min, params.x_max, params.Nx);
    y_grid = linspace(params.y_min, params.y_max, params.Ny);
    [X, Y] = meshgrid(x_grid, y_grid);

    % Initial condition
    f0 = setupInitialCondition(X, Y, params, Nv);

    % Time array
    times = linspace(0, params.T_final, params.Nt+1);

    % Solve
    f_all = solveDiscreteOrdinateAnalysis(f0, params, x_grid, y_grid, v_grid, times);

    % Visualize ray effects at different times
    figure(2);

    time_indices = [1, round(params.Nt/3), round(2*params.Nt/3), params.Nt+1];

    for i = 1:length(time_indices)
        ti = time_indices(i);
        t = times(ti);

        % Compute mass density
        rho = sum(f_all(:,:,:,ti), 3) * (2*pi/Nv);

        subplot(2, 2, i);
        imagesc(x_grid, y_grid, rho);
        colorbar;
        title(sprintf('t = %.3f', t));
        xlabel('x'); ylabel('y');
        axis equal tight;

        % Overlay characteristic lines
        hold on;
        plotCharacteristicLines(x_grid, y_grid, v_grid, t, params);
        hold off;
    end

    sgtitle('Ray Effect Development Over Time (8 Directions)');

    % Detailed ray structure analysis
    analyzeRayStructure(f_all, x_grid, y_grid, v_grid, times, params);

end

function plotCharacteristicLines(x_grid, y_grid, v_grid, t, params)
% PLOTCHARACTERISTICLINES Plot characteristic lines for discrete directions.

    % Plot characteristics passing through initial center
    x0 = params.x0;
    y0 = params.y0;

    for k = 1:size(v_grid, 1)
        vx = v_grid(k, 1);
        vy = v_grid(k, 2);

        % Characteristic line: x = x0 + vx*t, y = y0 + vy*t
        x_char = x0 + vx * t;
        y_char = y0 + vy * t;

        % Plot line from center to current position
        plot([x0, x_char], [y0, y_char], 'w--', 'LineWidth', 1.5);

        % Mark direction
        if abs(x_char) < max(abs(x_grid)) && abs(y_char) < max(abs(y_grid))
            plot(x_char, y_char, 'wo', 'MarkerSize', 6, 'MarkerFaceColor', 'red');
        end
    end

end

function analyzeRayStructure(f_all, x_grid, y_grid, v_grid, times, params)
% ANALYZERAYSTRUCTURE Detailed analysis of ray structure formation.

    fprintf('    Analyzing ray structure formation...\n');

    figure(3);

    % Analyze evolution of angular distribution at a fixed point
    % Choose a point off-center to see directional effects
    test_i = round(3*params.Ny/4);
    test_j = round(3*params.Nx/4);
    test_x = x_grid(test_j);
    test_y = y_grid(test_i);

    fprintf('      Analysis point: (%.3f, %.3f)\n', test_x, test_y);

    subplot(2, 2, 1);
    % Plot angular distribution evolution at test point
    for ti = 1:4:length(times)
        f_angular = squeeze(f_all(test_i, test_j, :, ti));
        theta = atan2(v_grid(:,2), v_grid(:,1));
        [theta_sorted, idx] = sort(theta);
        f_sorted = f_angular(idx);

        plot(theta_sorted, f_sorted, '-o', 'DisplayName', sprintf('t=%.3f', times(ti)));
        hold on;
    end
    xlabel('Angle (rad)'); ylabel('f(v)');
    title('Angular Distribution Evolution');
    legend('Location', 'best');
    grid on;

    % Analyze information propagation speed
    subplot(2, 2, 2);
    analyzeInformationSpeed(f_all, x_grid, y_grid, times, params);

    % Show directional flux patterns
    subplot(2, 2, 3);
    showDirectionalFlux(f_all, x_grid, y_grid, v_grid, times, params);

    % Demonstrate missing angular information
    subplot(2, 2, 4);
    demonstrateMissingAngularInfo(v_grid, params);

end

function analyzeInformationSpeed(f_all, x_grid, y_grid, times, params)
% ANALYZEINFORMATIONSPEED Analyze how information propagates.

    % Track the spreading of mass from initial center
    [X, Y] = meshgrid(x_grid, y_grid);

    max_radius = zeros(size(times));
    for ti = 1:length(times)
        rho = sum(f_all(:,:,:,ti), 3);

        % Find radius containing 95% of mass
        r = sqrt((X - params.x0).^2 + (Y - params.y0).^2);
        total_mass = sum(rho(:));

        r_values = linspace(0, max(r(:)), 100);
        cumulative_mass = zeros(size(r_values));

        for i = 1:length(r_values)
            mask = r <= r_values(i);
            cumulative_mass(i) = sum(rho(mask));
        end

        cumulative_fraction = cumulative_mass / total_mass;
        idx_95 = find(cumulative_fraction >= 0.95, 1, 'first');

        if ~isempty(idx_95)
            max_radius(ti) = r_values(idx_95);
        else
            max_radius(ti) = max(r(:));
        end
    end

    plot(times, max_radius, 'b-o', 'LineWidth', 2);
    hold on;
    plot(times, times, 'r--', 'LineWidth', 2); % Expected for |v|=1
    xlabel('Time'); ylabel('95% Mass Radius');
    title('Information Propagation Speed');
    legend('Discrete Ordinate', 'Expected (|v|=1)', 'Location', 'best');
    grid on;

    % Check if propagation is isotropic
    speed_discrete = max_radius(end) / times(end);
    fprintf('      Apparent propagation speed: %.3f (expected: 1.0)\n', speed_discrete);

end

function showDirectionalFlux(f_all, x_grid, y_grid, v_grid, times, params)
% SHOWDIRECTIONALFLUX Show directional flux patterns.

    % Compute flux in each direction at final time
    ti = length(times);

    % Sample flux along a line through center
    center_j = round(params.Nx/2);
    y_line = y_grid;

    flux_total = zeros(size(y_line));

    for k = 1:size(v_grid, 1)
        f_k = f_all(:, center_j, k, ti);
        vx = v_grid(k, 1);
        vy = v_grid(k, 2);

        % Directional flux component
        flux_k = f_k * vx; % x-component of flux
        flux_total = flux_total + flux_k;
    end

    plot(y_line, flux_total, 'b-', 'LineWidth', 2);
    xlabel('y'); ylabel('Net x-flux');
    title('Directional Flux Pattern');
    grid on;

    % Show discrete nature
    fprintf('      Flux computed from %d discrete directions only\n', size(v_grid, 1));

end

function demonstrateMissingAngularInfo(v_grid, params)
% DEMONSTRATEMISSINGANGULARINFO Show what angular information is missing.

    % Show discrete directions vs continuous coverage
    theta_discrete = atan2(v_grid(:,2), v_grid(:,1));
    theta_continuous = linspace(0, 2*pi, 1000);

    % Plot coverage
    polarplot(theta_discrete, ones(size(theta_discrete)), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'red');
    hold on;
    polarplot(theta_continuous, 0.5*ones(size(theta_continuous)), 'b-', 'LineWidth', 1);

    title('Angular Coverage: Discrete vs Continuous');
    legend('Discrete Ordinates', 'Continuous Coverage', 'Location', 'best');

    % Calculate angular gaps
    theta_sorted = sort(theta_discrete);
    gaps = diff([theta_sorted; theta_sorted(1) + 2*pi]);
    max_gap = max(gaps);

    fprintf('      Maximum angular gap: %.3f rad (%.1f degrees)\n', max_gap, max_gap*180/pi);
    fprintf('      Number of directions: %d\n', length(theta_discrete));

end

function analyzeInformationPropagation(params)
% ANALYZEINFORMATIONPROPAGATION Analyze how information propagates in DO method.

    fprintf('  Analyzing information propagation patterns...\n');

    % Create a localized initial condition (point source)
    x_grid = linspace(params.x_min, params.x_max, params.Nx);
    y_grid = linspace(params.y_min, params.y_max, params.Ny);
    [X, Y] = meshgrid(x_grid, y_grid);

    % Very localized initial condition
    params_local = params;
    params_local.sigma = 0.01; % Very sharp

    Nv = 8;
    theta = linspace(0, 2*pi, Nv+1);
    theta = theta(1:end-1);
    v_grid = [cos(theta)', sin(theta)'];

    f0 = setupInitialCondition(X, Y, params_local, Nv);

    % Short time evolution to see initial propagation
    times_short = linspace(0, 0.05, 6);
    f_all = solveDiscreteOrdinateAnalysis(f0, params_local, x_grid, y_grid, v_grid, times_short);

    figure(4);

    for i = 1:length(times_short)
        subplot(2, 3, i);
        rho = sum(f_all(:,:,:,i), 3) * (2*pi/Nv);

        % Use logarithmic scale to see propagation better
        rho_log = log10(rho + 1e-10);
        imagesc(x_grid, y_grid, rho_log);
        colorbar;
        title(sprintf('log(ρ), t = %.3f', times_short(i)));
        xlabel('x'); ylabel('y');
        axis equal tight;

        % Overlay characteristic rays
        if i > 1
            hold on;
            for k = 1:Nv
                vx = v_grid(k,1);
                vy = v_grid(k,2);
                t = times_short(i);

                % Ray from origin
                x_ray = [0, vx*t];
                y_ray = [0, vy*t];
                plot(x_ray, y_ray, 'w--', 'LineWidth', 2);
            end
            hold off;
        end
    end

    sgtitle('Information Propagation Along Discrete Rays');

    % Quantify directional bias
    quantifyDirectionalBias(f_all, x_grid, y_grid, v_grid, times_short);

end

function quantifyDirectionalBias(f_all, x_grid, y_grid, v_grid, times)
% QUANTIFYDIRECTIONALBIAS Quantify how much solution is biased toward discrete directions.

    fprintf('    Quantifying directional bias...\n');

    % Analyze final time
    ti = length(times);
    rho = sum(f_all(:,:,:,ti), 3);

    % Compute moments to detect directional bias
    [X, Y] = meshgrid(x_grid, y_grid);

    % Normalize coordinates relative to center of mass
    total_mass = sum(rho(:));
    x_cm = sum(sum(rho .* X)) / total_mass;
    y_cm = sum(sum(rho .* Y)) / total_mass;

    X_rel = X - x_cm;
    Y_rel = Y - y_cm;

    % Compute directional moments
    Mxx = sum(sum(rho .* X_rel.^2)) / total_mass;
    Myy = sum(sum(rho .* Y_rel.^2)) / total_mass;
    Mxy = sum(sum(rho .* X_rel .* Y_rel)) / total_mass;

    % Compute eccentricity and principal directions
    M = [Mxx, Mxy; Mxy, Myy];
    [eigvec, eigval] = eig(M);

    eccentricity = sqrt(1 - min(diag(eigval))/max(diag(eigval)));

    fprintf('      Eccentricity: %.6f (0=isotropic, ~1=highly directional)\n', eccentricity);

    % Principal direction
    [~, max_idx] = max(diag(eigval));
    principal_dir = atan2(eigvec(2, max_idx), eigvec(1, max_idx));

    fprintf('      Principal direction: %.3f rad (%.1f degrees)\n', principal_dir, principal_dir*180/pi);

    % Check alignment with discrete ordinate directions
    do_directions = atan2(v_grid(:,2), v_grid(:,1));
    angular_diffs = abs(principal_dir - do_directions);
    angular_diffs = min(angular_diffs, 2*pi - angular_diffs); % Handle periodicity

    [min_diff, closest_idx] = min(angular_diffs);

    fprintf('      Closest DO direction: %.3f rad, difference: %.3f rad (%.1f degrees)\n', ...
            do_directions(closest_idx), min_diff, min_diff*180/pi);

end

function compareWithExactCharacteristics(params)
% COMPAREWITHEXACTCHARACTERISTICS Compare DO method with exact characteristics.

    fprintf('  Comparing with exact characteristic solution...\n');

    % Setup for comparison
    Nv = 8;
    theta = linspace(0, 2*pi, Nv+1);
    theta = theta(1:end-1);
    v_grid = [cos(theta)', sin(theta)'];

    x_grid = linspace(params.x_min, params.x_max, params.Nx);
    y_grid = linspace(params.y_min, params.y_max, params.Ny);
    [X, Y] = meshgrid(x_grid, y_grid);

    % Initial condition
    f0_do = setupInitialCondition(X, Y, params, Nv);
    f0_exact = setupInitialCondition(X, Y, params, 64); % High resolution for exact

    % Solve with both methods
    times = linspace(0, params.T_final, params.Nt+1);

    % Discrete ordinate
    f_do = solveDiscreteOrdinateAnalysis(f0_do, params, x_grid, y_grid, v_grid, times);

    % Exact characteristics (high resolution)
    theta_exact = linspace(0, 2*pi, 65);
    theta_exact = theta_exact(1:end-1);
    v_grid_exact = [cos(theta_exact)', sin(theta_exact)'];
    f_exact = solveExactCharacteristics(f0_exact, params, x_grid, y_grid, v_grid_exact, times);

    % Compare at final time
    figure(5);

    ti = length(times);

    % Discrete ordinate result
    subplot(1, 3, 1);
    rho_do = sum(f_do(:,:,:,ti), 3) * (2*pi/Nv);
    imagesc(x_grid, y_grid, rho_do);
    colorbar; title('Discrete Ordinate');
    xlabel('x'); ylabel('y'); axis equal tight;

    % Exact result
    subplot(1, 3, 2);
    rho_exact = sum(f_exact(:,:,:,ti), 3) * (2*pi/64);
    imagesc(x_grid, y_grid, rho_exact);
    colorbar; title('Exact Characteristics');
    xlabel('x'); ylabel('y'); axis equal tight;

    % Difference
    subplot(1, 3, 3);
    % Interpolate exact to same resolution for comparison
    rho_exact_interp = rho_exact; % Same grid already
    diff = rho_do - rho_exact_interp;
    imagesc(x_grid, y_grid, diff);
    colorbar; title('Difference (DO - Exact)');
    xlabel('x'); ylabel('y'); axis equal tight;

    sgtitle(sprintf('Comparison at t = %.3f', times(ti)));

    % Quantify error
    l2_error = sqrt(mean(diff(:).^2));
    max_error = max(abs(diff(:)));

    fprintf('    L2 error: %.6f\n', l2_error);
    fprintf('    Max error: %.6f\n', max_error);

end

function investigateAngularCoupling(params)
% INVESTIGATEANGULARCOUPLING Investigate lack of angular coupling in DO method.

    fprintf('  Investigating angular coupling mechanisms...\n');

    % Setup test case with multiple direction initial condition
    Nv = 8;
    theta = linspace(0, 2*pi, Nv+1);
    theta = theta(1:end-1);
    v_grid = [cos(theta)', sin(theta)'];

    x_grid = linspace(params.x_min, params.x_max, params.Nx);
    y_grid = linspace(params.y_min, params.y_max, params.Ny);
    [X, Y] = meshgrid(x_grid, y_grid);

    % Create initial condition with specific directional structure
    f0 = setupDirectionalInitialCondition(X, Y, v_grid, params);

    times = linspace(0, params.T_final, params.Nt+1);
    f_all = solveDiscreteOrdinateAnalysis(f0, params, x_grid, y_grid, v_grid, times);

    figure(6);

    % Show evolution of directional distribution at center
    center_i = round(params.Ny/2);
    center_j = round(params.Nx/2);

    subplot(2, 2, 1);
    for ti = 1:5:length(times)
        f_center = squeeze(f_all(center_i, center_j, :, ti));
        plot(theta, f_center, '-o', 'DisplayName', sprintf('t=%.3f', times(ti)));
        hold on;
    end
    xlabel('Direction (rad)'); ylabel('f');
    title('Angular Distribution at Center');
    legend('Location', 'best'); grid on;

    % Show that directions evolve independently
    subplot(2, 2, 2);
    % Track two specific directions
    dir1 = 1; dir2 = 3;
    f_dir1 = squeeze(f_all(center_i, center_j, dir1, :));
    f_dir2 = squeeze(f_all(center_i, center_j, dir2, :));

    plot(times, f_dir1, 'r-o', 'DisplayName', sprintf('Direction %d (%.1f°)', dir1, theta(dir1)*180/pi));
    hold on;
    plot(times, f_dir2, 'b-s', 'DisplayName', sprintf('Direction %d (%.1f°)', dir2, theta(dir2)*180/pi));
    xlabel('Time'); ylabel('f');
    title('Independent Evolution of Directions');
    legend('Location', 'best'); grid on;

    % Demonstrate lack of angular diffusion
    subplot(2, 2, 3);
    demonstrateNoAngularDiffusion(f_all, theta, times);

    % Show angular correlation analysis
    subplot(2, 2, 4);
    analyzeAngularCorrelations(f_all, theta, times, params);

end

function f0 = setupDirectionalInitialCondition(X, Y, v_grid, params)
% SETUPDIRECTIONALINITIALCONDITION Create initial condition with directional structure.

    Nv = size(v_grid, 1);
    [Ny, Nx] = size(X);
    f0 = zeros(Ny, Nx, Nv);

    % Spatial Gaussian
    r2 = (X - params.x0).^2 + (Y - params.y0).^2;
    spatial_part = exp(-r2 / (2 * params.sigma^2));

    % Different directional structure: two peaks
    theta = atan2(v_grid(:,2), v_grid(:,1));
    directional_part = exp(-2*(theta - pi/4).^2) + exp(-2*(theta - 5*pi/4).^2);
    directional_part = directional_part / sum(directional_part);

    for k = 1:Nv
        f0(:,:,k) = spatial_part * directional_part(k);
    end

end

function demonstrateNoAngularDiffusion(f_all, theta, times)
% DEMONSTRATENOANGULARDIFFUSION Show that angular distribution doesn't diffuse.

    % Sample at a fixed spatial point away from center
    test_i = round(size(f_all,1) * 0.7);
    test_j = round(size(f_all,2) * 0.7);

    % Compute angular variance over time
    angular_variance = zeros(size(times));

    for ti = 1:length(times)
        f_angular = squeeze(f_all(test_i, test_j, :, ti));

        if sum(f_angular) > 1e-10
            % Compute angular mean and variance
            f_normalized = f_angular / sum(f_angular);

            % Ensure consistent vector orientations
            f_normalized = f_normalized(:);
            theta_vec = theta(:);

            % Convert to complex representation for circular statistics
            z = sum(f_normalized .* exp(1i * theta_vec));
            angular_variance(ti) = 1 - abs(z); % 0 = concentrated, 1 = uniform
        else
            angular_variance(ti) = 0;
        end
    end

    plot(times, angular_variance, 'g-o', 'LineWidth', 2);
    xlabel('Time'); ylabel('Angular Variance');
    title('Angular Diffusion (Should Increase for Physical Diffusion)');
    grid on;

    fprintf('      Final angular variance: %.6f (no increase = no angular diffusion)\n', ...
            angular_variance(end));

end

function analyzeAngularCorrelations(f_all, theta, times, params)
% ANALYZEANGULARCORRELATIONS Analyze correlations between angular directions.

    % Sample multiple spatial points
    [Ny, Nx, Nv, Nt] = size(f_all);

    % Choose several points
    points_i = round(linspace(Ny*0.3, Ny*0.7, 5));
    points_j = round(linspace(Nx*0.3, Nx*0.7, 5));

    % Compute cross-correlations between directions
    ti = length(times); % Final time

    correlations = zeros(Nv, Nv);

    for pi = 1:length(points_i)
        for pj = 1:length(points_j)
            i = points_i(pi);
            j = points_j(pj);

            f_point = squeeze(f_all(i, j, :, ti));

            if sum(f_point) > 1e-10
                % Compute correlation matrix
                for v1 = 1:Nv
                    for v2 = 1:Nv
                        if v1 == v2
                            correlations(v1, v2) = correlations(v1, v2) + 1;
                        else
                            % Simple correlation: product of normalized values
                            corr = f_point(v1) * f_point(v2) / (sum(f_point)^2);
                            correlations(v1, v2) = correlations(v1, v2) + corr;
                        end
                    end
                end
            end
        end
    end

    % Normalize
    correlations = correlations / (length(points_i) * length(points_j));

    % Plot correlation matrix
    imagesc(correlations);
    colorbar;
    xlabel('Direction Index'); ylabel('Direction Index');
    title('Angular Direction Correlations');

    % Check off-diagonal correlations
    off_diag = correlations - diag(diag(correlations));
    mean_off_diag = mean(abs(off_diag(:)));

    fprintf('      Mean off-diagonal correlation: %.6f (should be ~0 for independent directions)\n', ...
            mean_off_diag);

end

% Helper functions

function f0 = setupInitialCondition(X, Y, params, Nv)
% SETUPINITIALCONDITION Setup initial condition for analysis.

    r2 = (X - params.x0).^2 + (Y - params.y0).^2;
    f_spatial = exp(-r2 / (2 * params.sigma^2));
    f0 = repmat(f_spatial, [1, 1, Nv]);

end

function f_all = solveDiscreteOrdinateAnalysis(f0, params, x_grid, y_grid, v_grid, times)
% SOLVEDISCRETEORDINATEANALYSIS Solve discrete ordinate for analysis.

    Nv = size(v_grid, 1);
    Nt = length(times);
    f_all = zeros(params.Ny, params.Nx, Nv, Nt);
    f_all(:,:,:,1) = f0;

    dx = x_grid(2) - x_grid(1);
    dy = y_grid(2) - y_grid(1);

    for n = 1:Nt-1
        f_current = f_all(:,:,:,n);
        f_new = f_current;

        for k = 1:Nv
            vx = v_grid(k,1);
            vy = v_grid(k,2);

            % Simple upwind scheme
            f_temp = f_current(:,:,k);

            % x-direction
            if abs(vx) > 1e-12
                if vx > 0
                    for i = 2:size(f_temp,2)
                        f_temp(:,i) = f_temp(:,i) - params.dt * vx * (f_temp(:,i) - f_temp(:,i-1)) / dx;
                    end
                    f_temp(:,1) = f_temp(:,1) - params.dt * vx * f_temp(:,1) / dx;
                else
                    for i = size(f_temp,2)-1:-1:1
                        f_temp(:,i) = f_temp(:,i) - params.dt * vx * (f_temp(:,i+1) - f_temp(:,i)) / dx;
                    end
                    f_temp(:,end) = f_temp(:,end) + params.dt * abs(vx) * f_temp(:,end) / dx;
                end
            end

            % y-direction
            if abs(vy) > 1e-12
                if vy > 0
                    for j = 2:size(f_temp,1)
                        f_temp(j,:) = f_temp(j,:) - params.dt * vy * (f_temp(j,:) - f_temp(j-1,:)) / dy;
                    end
                    f_temp(1,:) = f_temp(1,:) - params.dt * vy * f_temp(1,:) / dy;
                else
                    for j = size(f_temp,1)-1:-1:1
                        f_temp(j,:) = f_temp(j,:) - params.dt * vy * (f_temp(j+1,:) - f_temp(j,:)) / dy;
                    end
                    f_temp(end,:) = f_temp(end,:) + params.dt * abs(vy) * f_temp(end,:) / dy;
                end
            end

            f_new(:,:,k) = f_temp;
        end

        f_all(:,:,:,n+1) = f_new;
    end

end

function f_exact = solveExactCharacteristics(f0, params, x_grid, y_grid, v_grid, times)
% SOLVEEXACTCHARACTERISTICS Solve using exact method of characteristics.

    Nv = size(v_grid, 1);
    Nt = length(times);
    f_exact = zeros(params.Ny, params.Nx, Nv, Nt);
    f_exact(:,:,:,1) = f0;

    [X, Y] = meshgrid(x_grid, y_grid);

    for n = 2:Nt
        t = times(n);

        for k = 1:Nv
            vx = v_grid(k,1);
            vy = v_grid(k,2);

            % Characteristic foot
            X_char = X - vx * t;
            Y_char = Y - vy * t;

            % Evaluate initial condition at characteristic foot
            r2 = (X_char - params.x0).^2 + (Y_char - params.y0).^2;
            f_char = exp(-r2 / (2 * params.sigma^2));

            f_exact(:,:,k,n) = f_char;
        end
    end

end