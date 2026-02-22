%% ex05_stability_regions.m
% Analyze and visualize stability regions for all integrator types
%
% Tests A-stability using the Dahlquist test equation: du/dt = λu
% A method is A-stable if |R(z)| ≤ 1 for all z with Re(z) ≤ 0
% where R(z) is the stability function and z = λ*dt
% Includes explicit, implicit, and IMEX integrators.

clear; close all; clc;

%% Setup for stability analysis
fprintf('=== STABILITY REGION ANALYSIS ===\n\n');

% Complex plane grid
lambda_real = linspace(-8, 2, 200);
lambda_imag = linspace(-5, 5, 200);
[LR, LI] = meshgrid(lambda_real, lambda_imag);
Z = LR + 1i * LI;

%% Define integrators to analyze
final_time = 5.0;
% Implicit integrators
implicit_integrators = {
    approx.odeint.BeIntegrator(final_time)
    approx.odeint.Bdf2Integrator(final_time)
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

% IMEX integrators (for A-stability, we consider their implicit part)
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

%% Compute stability functions
R = cell(n_methods, 1);

% Implicit methods
% Backward Euler: R(z) = 1/(1-z)
R{1} = 1 ./ (1 - Z);

% BDF2: R(z) = (3-2z)/(1 - 2z/3)^2 (corrected BDF2 stability function)
R{2} = (3 - 2*Z) ./ (1 - (2/3)*Z).^2;

% SDIRK2: R(z) = (1 + (1-2γ)z) / (1 - γz)^2
gamma2 = 1 - sqrt(2)/2;
R{3} = (1 + (1-2*gamma2)*Z) ./ (1 - gamma2*Z).^2;

% SDIRK3: Approximate stability function
gamma3 = 0.435866521508459;
R{4} = 1 ./ (1 - gamma3*Z).^3;  % Simplified approximation

% Explicit methods
% Forward Euler: R(z) = 1 + z
R{5} = 1 + Z;

% Heun (RK2): R(z) = 1 + z + z^2/2
R{6} = 1 + Z + Z.^2/2;

% RK4: R(z) = 1 + z + z^2/2 + z^3/6 + z^4/24
R{7} = 1 + Z + Z.^2/2 + Z.^3/6 + Z.^4/24;

% SSPRK3: Approximate stability function (3rd order)
R{8} = 1 + Z + Z.^2/2 + Z.^3/6;

% IMEX methods (analyze implicit part for A-stability)
% ARS111: Similar to Backward Euler for implicit part
R{9} = 1 ./ (1 - Z);

% ARS222: Similar to SDIRK2 for implicit part
gamma_ars2 = 1 - sqrt(2)/2;
R{10} = (1 + (1-2*gamma_ars2)*Z) ./ (1 - gamma_ars2*Z).^2;

% ARS443: More complex, approximate as 3rd order DIRK
R{11} = 1 ./ (1 - 0.5*Z).^3;  % Simplified

%% Plot stability regions
figure('Position', [100, 100, 1600, 1200]);

% Define colors for different types
colors_implicit = lines(length(implicit_integrators));
colors_explicit = autumn(length(explicit_integrators));
colors_imex = spring(length(imex_integrators));

% Plot all stability regions
subplot_positions = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];

% Check A-stability for all methods
is_A_stable = false(n_methods, 1);
left_half_plane = (LR <= 0);

for i = 1:n_methods
    subplot(3, 4, subplot_positions(i));
    
    % Compute magnitude of stability function
    R_mag = abs(R{i});
    
    % Plot stability region (where |R(z)| ≤ 1)
    contourf(LR, LI, R_mag, [0:0.1:1, 1.5:0.5:5], 'ShowText', 'off');
    hold on;
    
    % Highlight the stability boundary |R(z)| = 1
    contour(LR, LI, R_mag, [1, 1], 'k-', 'LineWidth', 3);
    
    % Highlight left half-plane (where A-stability is required)
    fill([-8, 0, 0, -8], [-5, -5, 5, 5], 'white', 'FaceAlpha', 0.1, ...
        'EdgeColor', 'red', 'LineWidth', 2, 'LineStyle', '--');
    
    % Color scheme: stable region in blue, unstable in red
    colormap(flipud(hot));
    colorbar;
    caxis([0, 3]);
    
    xlabel('Re(λ·dt)');
    ylabel('Im(λ·dt)');
    title(sprintf('%s (%s)', integrator_names{i}, integrator_types{i}));
    grid on;
    axis equal;
    xlim([-8, 2]);
    ylim([-5, 5]);
    
    % Check A-stability (only meaningful for implicit and IMEX methods)
    if strcmp(integrator_types{i}, 'Explicit')
        % Explicit methods are never A-stable
        is_A_stable(i) = false;
        text(-6, 4, 'NOT A-STABLE', 'FontSize', 10, 'FontWeight', 'bold', ...
            'Color', 'red', 'BackgroundColor', 'white');
        text(-6, 3, '(Explicit)', 'FontSize', 8, ...
            'Color', 'red', 'BackgroundColor', 'white');
    else
        % Check A-stability for implicit and IMEX methods
        is_A_stable(i) = all(R_mag(left_half_plane) <= 1.05);  % Small tolerance
        
        if is_A_stable(i)
            text(-6, 4, 'A-STABLE', 'FontSize', 10, 'FontWeight', 'bold', ...
                'Color', 'green', 'BackgroundColor', 'white');
        else
            text(-6, 4, 'NOT A-STABLE', 'FontSize', 10, 'FontWeight', 'bold', ...
                'Color', 'red', 'BackgroundColor', 'white');
        end
    end
    
    % Print result
    if strcmp(integrator_types{i}, 'Explicit')
        fprintf('%s (%s): Not A-stable (explicit method)\n', ...
            integrator_names{i}, integrator_types{i});
    else
        stability_str = is_A_stable(i);
        if stability_str
            result_str = "A-stable";
        else
            result_str = "Not A-stable";
        end
        fprintf('%s (%s): %s\n', ...
            integrator_names{i}, integrator_types{i}, result_str);
    end
end

%% Stability region comparison by method type
figure('Position', [200, 200, 1500, 500]);

method_types = {'Implicit', 'Explicit', 'IMEX'};
colors_by_type = {colors_implicit, colors_explicit, colors_imex};

for type_idx = 1:3
    subplot(1, 3, type_idx);
    
    type_name = method_types{type_idx};
    type_indices = find(strcmp(integrator_types, type_name));
    colors = colors_by_type{type_idx};
    
    if isempty(type_indices)
        continue;
    end
    
    % Plot multiple stability boundaries on same plot
    color_idx = 1;
    for i = type_indices'
        R_mag = abs(R{i});
        contour(LR, LI, R_mag, [1, 1], '-', 'Color', colors(color_idx, :), ...
            'LineWidth', 2, 'DisplayName', integrator_names{i});
        hold on;
        color_idx = color_idx + 1;
    end
    
    % Highlight left half-plane
    fill([-8, 0, 0, -8], [-5, -5, 5, 5], 'white', 'FaceAlpha', 0.1, ...
        'EdgeColor', 'red', 'LineWidth', 2, 'LineStyle', '--', ...
        'DisplayName', 'Left half-plane');
    
    xlabel('Re(λ·dt)');
    ylabel('Im(λ·dt)');
    title([type_name ' Methods - Stability Boundaries']);
    legend('Location', 'eastoutside', 'FontSize', 8);
    grid on;
    axis equal;
    xlim([-8, 2]);
    ylim([-5, 5]);
    
    % Add text about A-stability for this type
    if strcmp(type_name, 'Explicit')
        text(-6, 4, 'Explicit methods:', 'FontSize', 10, 'FontWeight', 'bold', ...
            'BackgroundColor', 'white');
        text(-6, 3.5, 'Never A-stable', 'FontSize', 9, ...
            'BackgroundColor', 'white');
        text(-6, 3, 'Stability limited', 'FontSize', 9, ...
            'BackgroundColor', 'white');
    else
        % Count A-stable methods in this type
        type_a_stable = sum(is_A_stable(type_indices));
        text(-6, 4, sprintf('%d/%d A-stable', type_a_stable, length(type_indices)), ...
            'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'white');
    end
end

%% Practical stability test with Dahlquist equation
fprintf('\n=== PRACTICAL STABILITY TEST ===\n\n');

% Test parameters
lambda_test = -10;  % Stiff parameter
dt_test = 0.5;      % Large time step to test stability
u0 = 1;

% Exact solution: u(t) = exp(λt)
exact_solution = @(t) exp(lambda_test * t);

fprintf('Testing with λ = %g, dt = %g\n', lambda_test, dt_test);
fprintf('Theoretical z = λ·dt = %g\n\n', lambda_test * dt_test);

% Test with subset of methods (representative from each type)
test_indices = [1, 3, 5, 7, 9, 10];  % Be, Sdirk2, Fe, Rk4, Ars111, Ars222
test_integrators = all_integrators(test_indices);
test_names = integrator_names(test_indices);
test_types = integrator_types(test_indices);

% Plot practical test results
figure('Position', [300, 300, 1200, 800]);

colors_test = lines(length(test_integrators));
n_steps = ceil(final_time / dt_test);
t_values = (0:n_steps) * dt_test;

% Storage for numerical solutions
U_solutions = zeros(length(test_integrators), n_steps + 1);

% Solution evolution plot
subplot(2, 2, 1);
for i = 1:length(test_integrators)
    integrator = test_integrators{i};
    
    % Initialize
    integrator.update(u0);
    
    % Store initial condition
    U_solutions(i, 1) = u0;
    
    % Time stepping
    stable = true;
    for k = 1:n_steps
        integrator.setTimeStep(dt=dt_test);
        
        % Set coefficients for BDF methods
        if isa(integrator, 'approx.odeint.BdfIntegrator')
            integrator.setCoefficients();
        end
        
        % Call appropriate step method based on integrator type
        if strcmp(test_types{i}, 'Implicit')
            U_new = integrator.step(L=lambda_test, S=[], M=[]);
        elseif strcmp(test_types{i}, 'Explicit')
            % For explicit methods, treat λu as nonlinear term
            F_eval = lambda_test * integrator.U0{1};
            U_new = integrator.step(L=[], F=F_eval, S=[], M=[]);
        else % IMEX
            U_new = integrator.step(L=lambda_test, F=[], S=[], M=[]);
        end
        
        integrator.update(U_new);
        integrator.advance();
        U_solutions(i, k + 1) = U_new;
        
        % Check for instability
        if abs(U_new) > 1e10
            stable = false;
            U_solutions(i, k+1:end) = NaN;
            break;
        end
    end
    
    % Plot solution
    valid_indices = ~isnan(U_solutions(i, :));
    if any(valid_indices)
        semilogy(t_values(valid_indices), abs(U_solutions(i, valid_indices)), ...
            '-', 'Color', colors_test(i, :), 'LineWidth', 2, ...
            'DisplayName', [test_names{i} ' (' test_types{i} ')']);
        hold on;
    end
    
    % Report stability
    final_value = U_solutions(i, end);
    exact_final = exact_solution(final_time);
    
    fprintf('%s (%s):\n', test_names{i}, test_types{i});
    if stable && ~isnan(final_value)
        fprintf('  Final value: %.2e (exact: %.2e)\n', final_value, exact_final);
        fprintf('  Stable: Yes\n');
        fprintf('  Error: %.2e\n\n', abs(final_value - exact_final));
    else
        fprintf('  Unstable (solution blew up)\n\n');
    end
end

% Add exact solution
U_exact = exact_solution(t_values);
semilogy(t_values, abs(U_exact), 'k--', 'LineWidth', 2, 'DisplayName', 'Exact');

xlabel('Time');
ylabel('|u(t)|');
title('Dahlquist Test: Solution Evolution');
legend('Location', 'best', 'FontSize', 8);
grid on;

% Error evolution plot
subplot(2, 2, 2);
for i = 1:length(test_integrators)
    U_exact_vals = exact_solution(t_values);
    error = abs(U_solutions(i, :) - U_exact_vals);
    
    valid_indices = ~isnan(error) & error > 0;
    if any(valid_indices)
        semilogy(t_values(valid_indices), error(valid_indices), ...
            '-', 'Color', colors_test(i, :), 'LineWidth', 2, ...
            'DisplayName', [test_names{i} ' (' test_types{i} ')']);
        hold on;
    end
end

xlabel('Time');
ylabel('Absolute Error');
title('Error Evolution');
legend('Location', 'best', 'FontSize', 8);
grid on;

% Stability summary by method type
subplot(2, 2, 3);
type_stability = zeros(3, 1);  % [Implicit, Explicit, IMEX]
type_counts = zeros(3, 1);

for i = 1:length(test_integrators)
    if strcmp(test_types{i}, 'Implicit')
        type_idx = 1;
    elseif strcmp(test_types{i}, 'Explicit')
        type_idx = 2;
    else
        type_idx = 3;
    end
    
    type_counts(type_idx) = type_counts(type_idx) + 1;
    
    final_value = U_solutions(i, end);
    if ~isnan(final_value) && abs(final_value) < 1e10
        type_stability(type_idx) = type_stability(type_idx) + 1;
    end
end

type_names = {'Implicit', 'Explicit', 'IMEX'};
bar_colors = [0.3, 0.6, 0.9; 0.9, 0.6, 0.3; 0.6, 0.9, 0.3];

bar_handle = bar(1:3, type_stability, 'FaceColor', 'flat');
bar_handle.CData = bar_colors;

set(gca, 'XTick', 1:3);
set(gca, 'XTickLabel', type_names);
ylabel('Number of Stable Methods');
title('Stability by Method Type');
grid on;

% Add text annotations
for i = 1:3
    if type_counts(i) > 0
        text(i, type_stability(i) + 0.1, sprintf('%d/%d', type_stability(i), type_counts(i)), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
end

% Theoretical vs practical stability comparison
subplot(2, 2, 4);
theoretical_stable = is_A_stable(test_indices);
practical_stable = ~isnan(U_solutions(:, end)) & abs(U_solutions(:, end)) < 1e10;

% Create scatter plot
for i = 1:length(test_integrators)
    if strcmp(test_types{i}, 'Implicit')
        color = [0.3, 0.6, 0.9];
        marker = 'o';
    elseif strcmp(test_types{i}, 'Explicit')
        color = [0.9, 0.6, 0.3];
        marker = 's';
    else
        color = [0.6, 0.9, 0.3];
        marker = '^';
    end
    
    scatter(theoretical_stable(i), practical_stable(i), 100, color, ...
        marker, 'filled', 'DisplayName', [test_names{i} ' (' test_types{i} ')']);
    hold on;
end

set(gca, 'XTick', [0, 1]);
set(gca, 'XTickLabel', {'Not A-stable', 'A-stable'});
set(gca, 'YTick', [0, 1]);
set(gca, 'YTickLabel', {'Unstable', 'Stable'});
xlabel('Theoretical A-stability');
ylabel('Practical Stability (λdt = -5)');
title('Theory vs Practice');
legend('Location', 'best', 'FontSize', 8);
grid on;

%% Summary analysis
fprintf('=== STABILITY ANALYSIS SUMMARY ===\n\n');

fprintf('A-stability analysis:\n');
for type_idx = 1:3
    type_name = method_types{type_idx};
    type_indices = find(strcmp(integrator_types, type_name));
    
    if isempty(type_indices)
        continue;
    end
    
    if strcmp(type_name, 'Explicit')
        fprintf('%s methods: Never A-stable (inherent limitation)\n', type_name);
    else
        n_a_stable = sum(is_A_stable(type_indices));
        n_total = length(type_indices);
        fprintf('%s methods: %d/%d are A-stable\n', type_name, n_a_stable, n_total);
        
        % List which ones are A-stable
        a_stable_indices = type_indices(is_A_stable(type_indices));
        if ~isempty(a_stable_indices)
            fprintf('  A-stable: ');
            for i = 1:length(a_stable_indices)
                fprintf('%s', integrator_names{a_stable_indices(i)});
                if i < length(a_stable_indices)
                    fprintf(', ');
                end
            end
            fprintf('\n');
        end
    end
end

fprintf('\nPractical stability test (λdt = %g):\n', lambda_test * dt_test);
fprintf('Results show whether methods can handle large time steps for stiff problems.\n');
fprintf('A-stable methods should remain stable, while explicit methods typically fail.\n');

% Compare theoretical and practical results
fprintf('\nTheory vs Practice correlation:\n');
correlation_correct = 0;
for i = 1:length(test_integrators)
    theoretical = is_A_stable(test_indices(i));
    practical = ~isnan(U_solutions(i, end)) && abs(U_solutions(i, end)) < 1e10;
    
    if (theoretical && practical) || (~theoretical && ~practical)
        correlation_correct = correlation_correct + 1;
    end
end

fprintf('Methods where theory matches practice: %d/%d\n', ...
    correlation_correct, length(test_integrators));

fprintf('\nNote: Explicit methods are inherently limited by stability constraints.\n');
fprintf('      For stiff problems, implicit or IMEX methods are preferred.\n');