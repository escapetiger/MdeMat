%% ex04_van_der_pol.m
% Solve the Van der Pol oscillator using multiple integrator types
%
% Problem: d²u/dt² - μ(1-u²)du/dt + u = 0
% Rewritten as first-order system:
%   du₁/dt = u₂
%   du₂/dt = μ(1-u₁²)u₂ - u₁
%
% This is a classic nonlinear stiff oscillator problem.
% Tests implicit, explicit, and IMEX integrators.

clear; close all; clc;

%% Problem Setup
mu = 5;                  % Stiffness parameter (try 1, 5, 10)
final_time = 20;         % Long integration to see limit cycle

% Initial conditions
u0 = [2; 0];             % [position; velocity]

% For implicit methods: treat everything implicitly
L_implicit = sparse([0, 1; -1, 0]);  % Linear part
S_implicit = @(u) [0; mu * (1 - u(1)^2) * u(2)];  % Nonlinear source

% For explicit methods: combine everything
F_explicit = @(u) [u(2); mu * (1 - u(1)^2) * u(2) - u(1)];

% For IMEX methods: split linear and nonlinear parts
L_imex = sparse([0, 1; -1, 0]);      % Linear part (implicit)
F_imex = @(u) [0; mu * (1 - u(1)^2) * u(2)];  % Nonlinear part (explicit)

fprintf('=== VAN DER POL OSCILLATOR ===\n');
fprintf('Parameter μ = %g\n', mu);
fprintf('Integration time: 0 to %g\n', final_time);

%% Setup integrators
% Different time steps based on stability requirements
dt_implicit = 0.05;      % Larger step for implicit methods
dt_explicit = 0.01;      % Smaller step for explicit methods
dt_imex = 0.02;          % Intermediate step for IMEX methods

% Implicit integrators
implicit_integrators = {
    approx.odeint.BeIntegrator(final_time)
    approx.odeint.Bdf2Integrator(final_time)
    approx.odeint.Bdf3Integrator(final_time)
    approx.odeint.Sdirk2Integrator(final_time)
    approx.odeint.Sdirk3Integrator(final_time)
};

% Explicit integrators
explicit_integrators = {
    approx.odeint.FeIntegrator(final_time)
    approx.odeint.HeunIntegrator(final_time)
    approx.odeint.Exrk4Integrator(final_time)
    approx.odeint.Ssprk3Integrator(final_time)
};

% IMEX integrators
imex_integrators = {
    approx.odeint.Ars111Integrator(final_time)
    approx.odeint.Ars222Integrator(final_time)
    approx.odeint.Ars443Integrator(final_time)
};

% Combine all integrators
all_integrators = [implicit_integrators; explicit_integrators; imex_integrators];
n_methods = length(all_integrators);

% Extract integrator names and types
integrator_names = cell(n_methods, 1);
integrator_types = cell(n_methods, 1);

for i = 1:n_methods
    class_name = class(all_integrators{i});
    name_parts = strsplit(class_name, '.');
    base_name = name_parts{end};
    integrator_names{i} = strrep(base_name, 'Integrator', '');
    
    % Determine type
    if i <= length(implicit_integrators)
        integrator_types{i} = 'Implicit';
    elseif i <= length(implicit_integrators) + length(explicit_integrators)
        integrator_types{i} = 'Explicit';
    else
        integrator_types{i} = 'IMEX';
    end
end

%% Solve with all integrators
fprintf('\nSolving with multiple integrator types...\n');
fprintf('Implicit methods dt: %.3f\n', dt_implicit);
fprintf('Explicit methods dt: %.3f\n', dt_explicit);
fprintf('IMEX methods dt: %.3f\n\n', dt_imex);

% Storage for solutions
U_solutions = cell(n_methods, 1);
t_solutions = cell(n_methods, 1);
computational_times = zeros(n_methods, 1);

for i = 1:n_methods
    fprintf('Solving with %s (%s)...\n', integrator_names{i}, integrator_types{i});
    
    integrator = all_integrators{i};
    
    % Choose appropriate time step based on method type
    if strcmp(integrator_types{i}, 'Explicit')
        dt = dt_explicit;
    elseif strcmp(integrator_types{i}, 'IMEX')
        dt = dt_imex;
    else
        dt = dt_implicit;
    end
    
    n_steps = ceil(final_time / dt);
    t_values = (0:n_steps) * dt;
    
    % Initialize solution history
    integrator.update(u0);
    
    % Storage for this method
    U_current_method = zeros(2, n_steps + 1);
    U_current_method(:, 1) = u0;
    
    % Time stepping
    tic;
    for k = 1:n_steps
        % Get current solution for nonlinear term evaluation
        U_current = integrator.U0{1, 1};
        
        integrator.setTimeStep(dt=dt);
        
        % Set coefficients for BDF methods
        if isa(integrator, 'approx.odeint.BdfIntegrator')
            integrator.setCoefficients();
        end
        
        % Call appropriate step method based on integrator type
        if strcmp(integrator_types{i}, 'Implicit')
            % Treat nonlinear term as source (simplified implicit treatment)
            S_current = S_implicit(U_current);
            U_new = integrator.step(L=L_implicit, S=S_current, M=[]);
        elseif strcmp(integrator_types{i}, 'Explicit')
            % Treat everything explicitly
            F_current = F_explicit(U_current);
            U_new = integrator.step(L=[], F=F_current, S=[], M=[]);
        else % IMEX
            % Split: linear part implicit, nonlinear part explicit
            F_current = F_imex(U_current);
            U_new = integrator.step(L=L_imex, F=F_current, S=[], M=[]);
        end
        
        integrator.update(U_new);
        integrator.advance();
        U_current_method(:, k + 1) = U_new;
        
        % Progress indicator for longer runs
        if mod(k, round(n_steps/10)) == 0
            fprintf('  Progress: %d%%\n', round(100*k/n_steps));
        end
    end
    computational_times(i) = toc;
    
    % Store results
    U_solutions{i} = U_current_method;
    t_solutions{i} = t_values;
    
    fprintf('  Completed in %.3f seconds (%d steps)\n', computational_times(i), n_steps);
end

%% Plot results
figure('Position', [100, 100, 1800, 1000]);

% Define colors for different types
colors_implicit = lines(length(implicit_integrators));
colors_explicit = autumn(length(explicit_integrators));
colors_imex = spring(length(imex_integrators));

% Time series - Position
subplot(3, 4, 1);
color_idx = [1, 1, 1]; % Index for each type
for i = 1:n_methods
    if strcmp(integrator_types{i}, 'Implicit')
        color = colors_implicit(color_idx(1), :);
        style = '-';
        color_idx(1) = color_idx(1) + 1;
    elseif strcmp(integrator_types{i}, 'Explicit')
        color = colors_explicit(color_idx(2), :);
        style = '--';
        color_idx(2) = color_idx(2) + 1;
    else % IMEX
        color = colors_imex(color_idx(3), :);
        style = '-.';
        color_idx(3) = color_idx(3) + 1;
    end
    
    plot(t_solutions{i}, U_solutions{i}(1, :), style, 'Color', color, ...
        'LineWidth', 1.2, 'DisplayName', [integrator_names{i} ' (' integrator_types{i} ')']);
    hold on;
end
xlabel('Time');
ylabel('u₁ (position)');
title('Position vs Time');
legend('Location', 'best', 'FontSize', 7);
grid on;

% Time series - Velocity
subplot(3, 4, 2);
color_idx = [1, 1, 1]; % Reset index
for i = 1:n_methods
    if strcmp(integrator_types{i}, 'Implicit')
        color = colors_implicit(color_idx(1), :);
        style = '-';
        color_idx(1) = color_idx(1) + 1;
    elseif strcmp(integrator_types{i}, 'Explicit')
        color = colors_explicit(color_idx(2), :);
        style = '--';
        color_idx(2) = color_idx(2) + 1;
    else % IMEX
        color = colors_imex(color_idx(3), :);
        style = '-.';
        color_idx(3) = color_idx(3) + 1;
    end
    
    plot(t_solutions{i}, U_solutions{i}(2, :), style, 'Color', color, ...
        'LineWidth', 1.2, 'DisplayName', [integrator_names{i} ' (' integrator_types{i} ')']);
    hold on;
end
xlabel('Time');
ylabel('u₂ (velocity)');
title('Velocity vs Time');
legend('Location', 'best', 'FontSize', 7);
grid on;

% Phase portraits by method type
method_types = {'Implicit', 'Explicit', 'IMEX'};
colors_by_type = {colors_implicit, colors_explicit, colors_imex};

for type_idx = 1:3
    subplot(3, 4, 2 + type_idx);
    
    type_name = method_types{type_idx};
    type_indices = find(strcmp(integrator_types, type_name));
    colors = colors_by_type{type_idx};
    
    color_idx = 1;
    for i = type_indices'
        plot(U_solutions{i}(1, :), U_solutions{i}(2, :), '-', ...
            'Color', colors(color_idx, :), 'LineWidth', 1.5, ...
            'DisplayName', integrator_names{i});
        hold on;
        color_idx = color_idx + 1;
    end
    
    plot(u0(1), u0(2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', ...
        'DisplayName', 'Initial point');
    xlabel('u₁ (position)');
    ylabel('u₂ (velocity)');
    title([type_name ' Methods - Phase Portrait']);
    legend('Location', 'best', 'FontSize', 7);
    grid on;
    axis equal;
end

% Energy-like quantity evolution
subplot(3, 4, 6);
color_idx = [1, 1, 1]; % Reset index
for i = 1:n_methods
    if strcmp(integrator_types{i}, 'Implicit')
        color = colors_implicit(color_idx(1), :);
        style = '-';
        color_idx(1) = color_idx(1) + 1;
    elseif strcmp(integrator_types{i}, 'Explicit')
        color = colors_explicit(color_idx(2), :);
        style = '--';
        color_idx(2) = color_idx(2) + 1;
    else % IMEX
        color = colors_imex(color_idx(3), :);
        style = '-.';
        color_idx(3) = color_idx(3) + 1;
    end
    
    energy = 0.5 * (U_solutions{i}(1, :).^2 + U_solutions{i}(2, :).^2);
    plot(t_solutions{i}, energy, style, 'Color', color, ...
        'LineWidth', 1.2, 'DisplayName', [integrator_names{i} ' (' integrator_types{i} ')']);
    hold on;
end
xlabel('Time');
ylabel('Energy-like quantity');
title('Pseudo-energy Evolution');
legend('Location', 'best', 'FontSize', 7);
grid on;

% Computational cost comparison
subplot(3, 4, 7);
color_idx = [1, 1, 1]; % Reset index
bar_colors = zeros(n_methods, 3);
for i = 1:n_methods
    if strcmp(integrator_types{i}, 'Implicit')
        bar_colors(i, :) = colors_implicit(color_idx(1), :);
        color_idx(1) = color_idx(1) + 1;
    elseif strcmp(integrator_types{i}, 'Explicit')
        bar_colors(i, :) = colors_explicit(color_idx(2), :);
        color_idx(2) = color_idx(2) + 1;
    else % IMEX
        bar_colors(i, :) = colors_imex(color_idx(3), :);
        color_idx(3) = color_idx(3) + 1;
    end
end

bar_handle = bar(1:n_methods, computational_times, 'FaceColor', 'flat');
bar_handle.CData = bar_colors;
set(gca, 'XTick', 1:n_methods);
set(gca, 'XTickLabel', integrator_names);
xtickangle(45);
ylabel('Computational Time (seconds)');
title('Computational Cost');
grid on;

% Method comparison - orbit overlay
subplot(3, 4, 8);
% Show comparison of all methods on same phase portrait
color_idx = [1, 1, 1];
for i = 1:min(n_methods, 6)  % Limit to first 6 for readability
    if strcmp(integrator_types{i}, 'Implicit')
        color = colors_implicit(color_idx(1), :);
        style = '-';
        if color_idx(1) < size(colors_implicit, 1)
            color_idx(1) = color_idx(1) + 1;
        end
    elseif strcmp(integrator_types{i}, 'Explicit')
        color = colors_explicit(color_idx(2), :);
        style = '--';
        if color_idx(2) < size(colors_explicit, 1)
            color_idx(2) = color_idx(2) + 1;
        end
    else % IMEX
        color = colors_imex(color_idx(3), :);
        style = '-.';
        if color_idx(3) < size(colors_imex, 1)
            color_idx(3) = color_idx(3) + 1;
        end
    end
    
    plot(U_solutions{i}(1, :), U_solutions{i}(2, :), style, ...
        'Color', color, 'LineWidth', 1.2, ...
        'DisplayName', [integrator_names{i} ' (' integrator_types{i} ')']);
    hold on;
end
plot(u0(1), u0(2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', ...
    'DisplayName', 'Initial point');
xlabel('u₁ (position)');
ylabel('u₂ (velocity)');
title('Method Comparison - Phase Portrait');
legend('Location', 'best', 'FontSize', 7);
grid on;
axis equal;

% Amplitude evolution (envelope detection)
subplot(3, 4, 9);
color_idx = [1, 1, 1];
for i = 1:n_methods
    if strcmp(integrator_types{i}, 'Implicit')
        color = colors_implicit(color_idx(1), :);
        style = '-';
        color_idx(1) = color_idx(1) + 1;
    elseif strcmp(integrator_types{i}, 'Explicit')
        color = colors_explicit(color_idx(2), :);
        style = '--';
        color_idx(2) = color_idx(2) + 1;
    else % IMEX
        color = colors_imex(color_idx(3), :);
        style = '-.';
        color_idx(3) = color_idx(3) + 1;
    end
    
    % Compute amplitude (distance from origin)
    amplitude = sqrt(U_solutions{i}(1, :).^2 + U_solutions{i}(2, :).^2);
    plot(t_solutions{i}, amplitude, style, 'Color', color, ...
        'LineWidth', 1.2, 'DisplayName', [integrator_names{i} ' (' integrator_types{i} ')']);
    hold on;
end
xlabel('Time');
ylabel('Amplitude');
title('Amplitude Evolution');
legend('Location', 'best', 'FontSize', 7);
grid on;

% Step size comparison
subplot(3, 4, 10);
step_sizes = [dt_implicit * ones(length(implicit_integrators), 1); ...
              dt_explicit * ones(length(explicit_integrators), 1); ...
              dt_imex * ones(length(imex_integrators), 1)];

bar_handle = bar(1:n_methods, step_sizes, 'FaceColor', 'flat');
bar_handle.CData = bar_colors;
set(gca, 'XTick', 1:n_methods);
set(gca, 'XTickLabel', integrator_names);
xtickangle(45);
ylabel('Time Step Size');
title('Step Size Used');
grid on;

% Final orbit comparison - zoomed
subplot(3, 4, 11);
% Show final portion of orbits to check steady state
final_portion = 0.8;  % Last 20% of simulation
color_idx = [1, 1, 1];
for i = 1:n_methods
    if strcmp(integrator_types{i}, 'Implicit')
        color = colors_implicit(color_idx(1), :);
        style = '-';
        color_idx(1) = color_idx(1) + 1;
    elseif strcmp(integrator_types{i}, 'Explicit')
        color = colors_explicit(color_idx(2), :);
        style = '--';
        color_idx(2) = color_idx(2) + 1;
    else % IMEX
        color = colors_imex(color_idx(3), :);
        style = '-.';
        color_idx(3) = color_idx(3) + 1;
    end
    
    n_pts = size(U_solutions{i}, 2);
    start_idx = round(final_portion * n_pts);
    plot(U_solutions{i}(1, start_idx:end), U_solutions{i}(2, start_idx:end), ...
        style, 'Color', color, 'LineWidth', 1.2, ...
        'DisplayName', [integrator_names{i} ' (' integrator_types{i} ')']);
    hold on;
end
xlabel('u₁ (position)');
ylabel('u₂ (velocity)');
title('Final Orbits (Steady State)');
legend('Location', 'best', 'FontSize', 7);
grid on;
axis equal;

% Efficiency comparison (accuracy vs cost)
subplot(3, 4, 12);
% Use final amplitude as a measure of accuracy consistency
final_amplitudes = zeros(n_methods, 1);
color_idx = [1, 1, 1];

for i = 1:n_methods
    final_amplitudes(i) = sqrt(U_solutions{i}(1, end)^2 + U_solutions{i}(2, end)^2);
    
    if strcmp(integrator_types{i}, 'Implicit')
        color = colors_implicit(color_idx(1), :);
        marker = 'o';
        color_idx(1) = color_idx(1) + 1;
    elseif strcmp(integrator_types{i}, 'Explicit')
        color = colors_explicit(color_idx(2), :);
        marker = 's';
        color_idx(2) = color_idx(2) + 1;
    else % IMEX
        color = colors_imex(color_idx(3), :);
        marker = '^';
        color_idx(3) = color_idx(3) + 1;
    end
    
    scatter(computational_times(i), final_amplitudes(i), 100, color, ...
        marker, 'filled', 'DisplayName', [integrator_names{i} ' (' integrator_types{i} ')']);
    hold on;
end

xlabel('Computational Time (seconds)');
ylabel('Final Amplitude');
title('Efficiency Comparison');
legend('Location', 'best', 'FontSize', 7);
grid on;

%% Analysis
fprintf('\n=== SOLUTION ANALYSIS ===\n');

% Check if solution has reached limit cycle (steady oscillation)
period_start = round(0.75 * length(t_solutions{1}));  % Look at last quarter

method_types = {'Implicit', 'Explicit', 'IMEX'};
fprintf('\nLimit cycle analysis by method type:\n');

for type_idx = 1:3
    type_name = method_types{type_idx};
    type_indices = find(strcmp(integrator_types, type_name));
    
    if isempty(type_indices)
        continue;
    end
    
    fprintf('\n--- %s METHODS ---\n', upper(type_name));
    
    for i = type_indices'
        u1_final = U_solutions{i}(1, period_start:end);
        
        % Find approximate amplitude and period
        [peaks, max_indices] = findpeaks(u1_final, 'MinPeakDistance', round(0.1 * length(u1_final)));
        if length(max_indices) >= 2
            dt_used = t_solutions{i}(2) - t_solutions{i}(1);
            estimated_period = mean(diff(max_indices)) * dt_used;
            estimated_amplitude = mean(peaks);
            
            fprintf('%s:\n', integrator_names{i});
            fprintf('  Estimated period: %.2f\n', estimated_period);
            fprintf('  Estimated amplitude: %.2f\n', estimated_amplitude);
            fprintf('  Final amplitude: %.2f\n', final_amplitudes(i));
            fprintf('  Computational time: %.3f s\n', computational_times(i));
        else
            fprintf('%s:\n', integrator_names{i});
            fprintf('  Solution may not have reached steady limit cycle yet.\n');
            fprintf('  Final amplitude: %.2f\n', final_amplitudes(i));
            fprintf('  Computational time: %.3f s\n', computational_times(i));
        end
    end
end

% Compare methods across types
fprintf('\n=== METHOD COMPARISON ACROSS TYPES ===\n');

% Find most accurate method (closest to expected limit cycle)
target_amplitude = 2.0;  % Approximate expected amplitude for this μ
amplitude_errors = abs(final_amplitudes - target_amplitude);
[~, most_accurate_idx] = min(amplitude_errors);

fprintf('Most accurate method: %s (%s)\n', ...
    integrator_names{most_accurate_idx}, integrator_types{most_accurate_idx});
fprintf('  Final amplitude error: %.4f\n', amplitude_errors(most_accurate_idx));

% Find most efficient method (best accuracy per unit time)
efficiency = 1 ./ (amplitude_errors .* computational_times);
[~, most_efficient_idx] = max(efficiency);

fprintf('Most efficient method: %s (%s)\n', ...
    integrator_names{most_efficient_idx}, integrator_types{most_efficient_idx});
fprintf('  Efficiency metric: %.2f\n', efficiency(most_efficient_idx));

% Summary by type
fprintf('\nSummary by method type:\n');
for type_idx = 1:3
    type_name = method_types{type_idx};
    type_indices = find(strcmp(integrator_types, type_name));
    
    if ~isempty(type_indices)
        avg_time = mean(computational_times(type_indices));
        avg_amplitude = mean(final_amplitudes(type_indices));
        
        if strcmp(type_name, 'Implicit')
            dt_used = dt_implicit;
        elseif strcmp(type_name, 'Explicit')
            dt_used = dt_explicit;
        else
            dt_used = dt_imex;
        end
        
        fprintf('%s methods:\n', type_name);
        fprintf('  Average computation time: %.3f s\n', avg_time);
        fprintf('  Average final amplitude: %.3f\n', avg_amplitude);
        fprintf('  Time step used: %.3f\n', dt_used);
        fprintf('  Number of methods: %d\n', length(type_indices));
    end
end

fprintf('\nVan der Pol oscillator analysis completed with μ = %g.\n', mu);