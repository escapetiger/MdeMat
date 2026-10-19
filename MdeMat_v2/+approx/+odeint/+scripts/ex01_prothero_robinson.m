%% ex01_prothero_robinson.m
% Solve the Prothero-Robinson equation using multiple integrator types
%
% Problem: du/dt = λ(u - φ(t)) + φ'(t)
% where φ(t) = exp(-t), λ = -100 (stiff parameter)
% Exact solution: u(t) = exp(-t) when u(0) = 1
%
% This is a standard test problem for stiff ODE solvers.
% Tests explicit, implicit, IMEX, and splitting integrators.

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
S = @(u, t) -lambda * phi(t) + phi_prime(t);  % Source term

% For explicit methods: combine everything into one term
F_explicit = @(u, t) lambda * (u - phi(t)) + phi_prime(t);

% For splitting methods: decompose the equation
% du/dt = λu + (-λφ(t) + φ'(t))
% Split 1: du/dt = λu (stiff part - treat implicitly)
% Split 2: du/dt = -λφ(t) + φ'(t) (non-stiff source - treat explicitly)
L_split = {lambda, []};
F_split = {[], []};
S_split = {@(u, t) -lambda * phi(t), @(u, t) phi_prime(t)};

% Alternative splitting for demonstration:
% Split 1: du/dt = λ(u - φ(t))
% Split 2: du/dt = φ'(t)
L_split_alt = {[], []};
F_split_alt = {@(u, t) lambda * (u - phi(t)), []};
S_split_alt = {[], @(u, t) phi_prime(t)};

% Exact solution for comparison
exact_solution = @(t) exp(-t);

%% Setup integrators
dt_implicit = 0.1;      % Larger step for implicit methods
dt_imex = 0.05;         % Moderate step for implicit-explicit methods
dt_explicit = 0.01;     % Smaller step for explicit methods (stability)
dt_splitting = 0.01;    % Moderate step for splitting methods

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
    approx.odeint.Ssprk2Integrator(final_time)
    approx.odeint.Ssprk3Integrator(final_time)
    approx.odeint.Exrk4Integrator(final_time)
};

imex_integrators = {
    approx.odeint.Ars111Integrator(final_time)
    approx.odeint.Ars222Integrator(final_time)
    approx.odeint.Ars443Integrator(final_time)
};

% Create splitting integrators
% Lie-Trotter splitting: implicit stiff part + explicit source
factor_integrators_lt1 = {
    approx.odeint.BeIntegrator(final_time)
    approx.odeint.FeIntegrator(final_time)
};
lie_trotter_1 = approx.odeint.LieTrotterIntegrator(factor_integrators_lt1);

% Strang splitting: implicit stiff part + explicit source
factor_integrators_st1 = {
    approx.odeint.BeIntegrator(final_time)
    approx.odeint.FeIntegrator(final_time)
};
strang_1 = approx.odeint.StrangIntegrator(factor_integrators_st1);

% Alternative Lie-Trotter: explicit nonlinear + explicit source
factor_integrators_lt2 = {
    approx.odeint.FeIntegrator(final_time)
    approx.odeint.FeIntegrator(final_time)
};
lie_trotter_2 = approx.odeint.LieTrotterIntegrator(factor_integrators_lt2);

% Alternative Strang: explicit nonlinear + explicit source
factor_integrators_st2 = {
    approx.odeint.FeIntegrator(final_time)
    approx.odeint.FeIntegrator(final_time)
};
strang_2 = approx.odeint.StrangIntegrator(factor_integrators_st2);

splitting_integrators = {
    lie_trotter_1
    strang_1
    lie_trotter_2
    strang_2
};

% Combine all integrators
all_integrators = [implicit_integrators; explicit_integrators; imex_integrators; splitting_integrators];
n_methods = length(all_integrators);

% Extract integrator names and types
integrator_names = cell(n_methods, 1);
integrator_types = cell(n_methods, 1);

for i = 1:n_methods
    if i <= length(implicit_integrators)
        class_name = class(all_integrators{i});
        name_parts = strsplit(class_name, '.');
        base_name = name_parts{end};
        integrator_names{i} = strrep(base_name, 'Integrator', '');
        integrator_types{i} = 'Implicit';
    elseif i <= length(implicit_integrators) + length(explicit_integrators)
        class_name = class(all_integrators{i});
        name_parts = strsplit(class_name, '.');
        base_name = name_parts{end};
        integrator_names{i} = strrep(base_name, 'Integrator', '');
        integrator_types{i} = 'Explicit';
    elseif i <= length(implicit_integrators) + length(explicit_integrators) + length(imex_integrators)
        class_name = class(all_integrators{i});
        name_parts = strsplit(class_name, '.');
        base_name = name_parts{end};
        integrator_names{i} = strrep(base_name, 'Integrator', '');
        integrator_types{i} = 'IMEX';
    else
        % Splitting methods
        splitting_idx = i - length(implicit_integrators) - length(explicit_integrators) - length(imex_integrators);
        switch splitting_idx
            case 1
                integrator_names{i} = 'LT-ImplEx';
                integrator_types{i} = 'Splitting';
            case 2
                integrator_names{i} = 'Strang-ImplEx';
                integrator_types{i} = 'Splitting';
            case 3
                integrator_names{i} = 'LT-ExplExpl';
                integrator_types{i} = 'Splitting';
            case 4
                integrator_names{i} = 'Strang-ExplExpl';
                integrator_types{i} = 'Splitting';
        end
    end
end

%% Solve with all integrators
fprintf('Solving Prothero-Robinson equation with multiple integrator types...\n');
fprintf('Implicit methods dt: %.3f\n', dt_implicit);
fprintf('Explicit methods dt: %.3f\n', dt_explicit);
fprintf('IMEX methods dt: %.3f\n', dt_imex);
fprintf('Splitting methods dt: %.3f\n\n', dt_splitting);

% Storage for solutions and time values
U_solutions = cell(n_methods, 1);
t_solutions = cell(n_methods, 1);

for i = 1:n_methods
    fprintf('Solving with %s (%s)...\n', integrator_names{i}, integrator_types{i});
    
    integrator = all_integrators{i};
    integrator.reset();
    
    % Choose appropriate time step based on method type
    if strcmp(integrator_types{i}, 'Implicit')
        dt = dt_implicit;
    elseif strcmp(integrator_types{i}, 'Explicit')
        dt = dt_explicit;
    elseif strcmp(integrator_types{i}, 'IMEX')
        dt = dt_imex;
    else % Splitting
        dt = dt_splitting;
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
        % Set timeline for integrator
        integrator.setTimeStep(dt, 1, 1);
        
        % Call appropriate step method based on integrator type
        if strcmp(integrator_types{i}, 'Implicit')
            U_new = integrator.step(L, S, []);
        elseif strcmp(integrator_types{i}, 'Explicit')
            % For explicit methods, treat everything as nonlinear term
            t_current = (k-1) * dt;
            F_eval = F_explicit(integrator.U0{1}, t_current);
            U_new = integrator.step([], F_eval, [], []);
        elseif strcmp(integrator_types{i}, 'IMEX')
            U_new = integrator.step(L, [], S, []);
        else % Splitting
            splitting_idx = i - length(implicit_integrators) - length(explicit_integrators) - length(imex_integrators);
            switch splitting_idx
                case {1, 2} % ImplicitExplicit splitting
                    U_new = integrator.step(L_split, F_split, S_split, []);
                case {3, 4} % ExplicitExplicit splitting
                    U_new = integrator.step(L_split_alt, F_split_alt, S_split_alt, []);
            end
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
colors_splitting = cool(length(splitting_integrators));

% Solution comparison
subplot(1, 3, 1);
% Plot exact solution
t_exact = linspace(0, final_time, 1000);
U_exact_fine = exact_solution(t_exact);
plot(t_exact, U_exact_fine, 'k-', 'LineWidth', 3, 'DisplayName', 'Exact');
hold on;

% Plot numerical solutions
color_idx = [1, 1, 1, 1]; % Index for each type
for i = 1:n_methods
    if strcmp(integrator_types{i}, 'Implicit')
        color = colors_implicit(color_idx(1), :);
        style = '-';
        color_idx(1) = color_idx(1) + 1;
    elseif strcmp(integrator_types{i}, 'Explicit')
        color = colors_explicit(color_idx(2), :);
        style = '--';
        color_idx(2) = color_idx(2) + 1;
    elseif strcmp(integrator_types{i}, 'IMEX')
        color = colors_imex(color_idx(3), :);
        style = '-.';
        color_idx(3) = color_idx(3) + 1;
    else % Splitting
        color = colors_splitting(color_idx(4), :);
        style = ':';
        color_idx(4) = color_idx(4) + 1;
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
color_idx = [1, 1, 1, 1]; % Reset index
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
    elseif strcmp(integrator_types{i}, 'IMEX')
        color = colors_imex(color_idx(3), :);
        style = '-.';
        color_idx(3) = color_idx(3) + 1;
    else % Splitting
        color = colors_splitting(color_idx(4), :);
        style = ':';
        color_idx(4) = color_idx(4) + 1;
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

color_idx = [1, 1, 1, 1]; % Reset index
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
    elseif strcmp(integrator_types{i}, 'IMEX')
        % IMEX methods: mix of explicit and implicit evaluations
        computational_costs(i) = n_steps_used * all_integrators{i}.nStages * 1.5;
    else % Splitting
        % Splitting methods: cost depends on factor integrators
        splitting_idx = i - length(implicit_integrators) - length(explicit_integrators) - length(imex_integrators);
        integrator = all_integrators{i};
        total_cost = 0;
        for j = 1:integrator.nFactors
            factor = integrator.getFactor(j);
            if isa(factor, 'approx.odeint.BeIntegrator')
                % Implicit factor
                factor_cost = factor.nStages * 2;
            else
                % Explicit factor
                factor_cost = factor.nStages;
            end
            total_cost = total_cost + factor_cost;
        end
        if mod(splitting_idx, 2) == 0 % Strang methods (even indices)
            computational_costs(i) = n_steps_used * total_cost * 1.5; % More sub-steps
        else % Lie-Trotter methods (odd indices)
            computational_costs(i) = n_steps_used * total_cost;
        end
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
    elseif strcmp(integrator_types{i}, 'IMEX')
        color = colors_imex(color_idx(3), :);
        marker = '^';
        color_idx(3) = color_idx(3) + 1;
    else % Splitting
        color = colors_splitting(color_idx(4), :);
        marker = 'd';
        color_idx(4) = color_idx(4) + 1;
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