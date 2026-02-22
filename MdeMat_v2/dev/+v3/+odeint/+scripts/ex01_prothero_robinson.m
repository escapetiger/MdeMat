%% ex01_prothero_robinson.m
% Solve the Prothero-Robinson equation using multiple integrator types
%
% Problem: du/dt = λ(u - φ(t)) + φ'(t)
% where φ(t) = exp(-t), λ = -100 (stiff parameter)
% Exact solution: u(t) = exp(-t) when u(0) = 1
%
% This is a standard test problem for stiff ODE solvers.
% Tests explicit, implicit, and IMEX integrators.

clear; close all; clc;

%% Problem Setup
lambda = -100;           % Stiffness parameter
phi = @(t) exp(-t);      % Exact solution
phi_prime = @(t) -exp(-t);
u0 = 1;                  % Initial condition
final_time = 2.0;

% For IMEX methods: split into stiff linear and non-stiff parts
L = lambda;              % Stiff linear operator (implicit)
F = [];                  % No nonlinear term
S = @(t) -lambda * phi(t) + phi_prime(t);  % Source term

% For explicit methods: combine everything into one term
F_explicit = @(u, t) lambda * (u - phi(t)) + phi_prime(t);

% Exact solution for comparison
exact_solution = @(t) exp(-t);

%% Setup integrators
dt_implicit = 0.1;      % Larger step for implicit methods
dt_imex = 0.05;         % Moderate step for implicit-explicit methods
dt_explicit = 0.01;     % Smaller step for explicit methods (stability)

% Create integrators with different time steps
implicit_integrators = {
    approx.odeint.BeIntegrator(final_time)
    approx.odeint.Bdf2Integrator(final_time)
    approx.odeint.Bdf3Integrator(final_time)
    approx.odeint.Sdirk2Integrator(final_time)
    approx.odeint.Sdirk3Integrator(final_time)
    approx.odeint.Sdirk4Integrator(final_time)
    approx.odeint.Esdirk3Integrator(final_time)
};

explicit_integrators = {
    approx.odeint.FeIntegrator(final_time)
    approx.odeint.HeunIntegrator(final_time)
    approx.odeint.Exrk4Integrator(final_time)
    approx.odeint.Ssprk3Integrator(final_time)
};

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
fprintf('Solving Prothero-Robinson equation with multiple integrator types...\n');
fprintf('Implicit methods dt: %.3f\n', dt_implicit);
fprintf('Explicit methods dt: %.3f\n', dt_explicit);
fprintf('IMEX methods dt: %.3f\n\n', dt_imex);

% Storage for solutions and time values
U_solutions = cell(n_methods, 1);
t_solutions = cell(n_methods, 1);

for i = 1:n_methods
    fprintf('Solving with %s (%s)...\n', integrator_names{i}, integrator_types{i});
    
    integrator = all_integrators{i};
    
    % Choose appropriate time step based on method type
    if strcmp(integrator_types{i}, 'Implicit')
        dt = dt_implicit;
    elseif strcmp(integrator_types{i}, 'Explicit')
        dt = dt_explicit;
    else
        dt = dt_imex;
    end
    
    n_steps = ceil(final_time / dt);
    t_values = (0:n_steps) * dt;
    
    % Initialize solution history
    integrator.update(u0);
    
    % Storage for this method
    U_current_method = zeros(1, n_steps + 1);
    U_current_method(1) = u0;
    
    % Time stepping
    for k = 1:n_steps
        integrator.timeline.setTimeStep(dt, 1, 1);
        
        % Set coefficients for BDF methods
        if isa(integrator, 'approx.odeint.BdfIntegrator')
            integrator.setCoefficients();
        end
        
        % Call appropriate step method based on integrator type
        if strcmp(integrator_types{i}, 'Implicit')
            U_new = integrator.step(L, S, []);
        elseif strcmp(integrator_types{i}, 'Explicit')
            % For explicit methods, treat everything as nonlinear term
            t_current = (k-1) * dt;
            F_eval = F_explicit(integrator.U0{1}, t_current);
            U_new = integrator.step([], F_eval, [], []);
        else % IMEX
            U_new = integrator.step(L, [], S, []);
        end
        
        integrator.update(U_new);
        integrator.advance();
        U_current_method(k + 1) = U_new;
    end
    
    % Store results
    U_solutions{i} = U_current_method;
    t_solutions{i} = t_values;
end

%% Plot results
figure('Position', [100, 100, 1500, 600]);

% Define colors for different types
colors_implicit = lines(length(implicit_integrators));
colors_explicit = autumn(length(explicit_integrators));
colors_imex = spring(length(imex_integrators));

% Solution comparison
subplot(1, 3, 1);
% Plot exact solution
t_exact = linspace(0, final_time, 1000);
U_exact_fine = exact_solution(t_exact);
plot(t_exact, U_exact_fine, 'k-', 'LineWidth', 3, 'DisplayName', 'Exact');
hold on;

% Plot numerical solutions
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
    
    plot(t_solutions{i}, U_solutions{i}, style, 'Color', color, ...
        'LineWidth', 1.5, 'DisplayName', [integrator_names{i} ' (' integrator_types{i} ')']);
end

xlabel('Time');
ylabel('u(t)');
title('Prothero-Robinson Equation Solution');
legend('Location', 'best', 'FontSize', 8);
grid on;

% Error comparison
subplot(1, 3, 2);
color_idx = [1, 1, 1]; % Reset index
for i = 1:n_methods
    % Interpolate exact solution to match numerical time points
    U_exact_interp = exact_solution(t_solutions{i});
    error = abs(U_solutions{i} - U_exact_interp);
    
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
    
    semilogy(t_solutions{i}, error, style, 'Color', color, ...
        'LineWidth', 1.5, 'DisplayName', [integrator_names{i} ' (' integrator_types{i} ')']);
    hold on;
end

xlabel('Time');
ylabel('Absolute Error');
title('Error Comparison');
legend('Location', 'best', 'FontSize', 8);
grid on;

% Work-precision diagram (error vs computational cost)
subplot(1, 3, 3);
final_errors = zeros(n_methods, 1);
computational_costs = zeros(n_methods, 1);

color_idx = [1, 1, 1]; % Reset index
for i = 1:n_methods
    % Compute final error
    U_exact_final = exact_solution(final_time);
    final_errors(i) = abs(U_solutions{i}(end) - U_exact_final);
    
    % Estimate computational cost (number of function evaluations)
    n_steps_used = length(t_solutions{i}) - 1;
    if strcmp(integrator_types{i}, 'Implicit')
        % Implicit methods: assume one linear solve per step/stage
        computational_costs(i) = n_steps_used * all_integrators{i}.nStages;
    elseif strcmp(integrator_types{i}, 'Explicit')
        % Explicit methods: one function evaluation per stage
        computational_costs(i) = n_steps_used * all_integrators{i}.nStages;
    else % IMEX
        % IMEX methods: mix of explicit and implicit evaluations
        computational_costs(i) = n_steps_used * all_integrators{i}.nStages * 1.5;
    end
    
    % Plot point
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
    
    loglog(computational_costs(i), final_errors(i), marker, 'Color', color, ...
        'MarkerSize', 8, 'MarkerFaceColor', color, 'LineWidth', 1.5, ...
        'DisplayName', [integrator_names{i} ' (' integrator_types{i} ')']);
    hold on;
end

xlabel('Computational Cost (Function Evaluations)');
ylabel('Final Error');
title('Work-Precision Diagram');
legend('Location', 'best', 'FontSize', 8);
grid on;

%% Error analysis
fprintf('\n=== ERROR ANALYSIS ===\n');
fprintf('Results at t = %.1f:\n', final_time);
fprintf('Exact solution: %.6f\n\n', exact_solution(final_time));

fprintf('%-15s | %-8s | %-12s | %-12s | %-8s\n', ...
    'Method', 'Type', 'Final Value', 'Error', 'Steps');
fprintf('%s\n', repmat('-', 1, 70));

for i = 1:n_methods
    final_value = U_solutions{i}(end);
    final_error = abs(final_value - exact_solution(final_time));
    n_steps_used = length(t_solutions{i}) - 1;
    
    fprintf('%-15s | %-8s | %12.6f | %12.2e | %8d\n', ...
        integrator_names{i}, integrator_types{i}, final_value, final_error, n_steps_used);
end