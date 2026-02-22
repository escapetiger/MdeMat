%% ex02_heat_equation.m
% Solve the 1D heat equation using method of lines with multiple integrator types
%
% Problem: ∂u/∂t = α ∂²u/∂x², x ∈ [0,1], t > 0
% Boundary conditions: u(0,t) = u(1,t) = 0
% Initial condition: u(x,0) = sin(πx)
% Exact solution: u(x,t) = exp(-π²αt) sin(πx)
%
% Tests explicit, implicit, and IMEX integrators for parabolic PDEs.

clear; close all; clc;

%% Problem Setup
alpha = 1.0;             % Thermal diffusivity
n = 100;                 % Number of interior grid points
dx = 1 / (n + 1);        % Spatial step size
x = linspace(dx, 1 - dx, n)';
final_time = 0.2;

% Initial condition
u0 = sin(pi * x);

% Spatial discretization using finite differences
e = ones(n, 1);
L = alpha / dx^2 * spdiags([e, -2*e, e], -1:1, n, n);

% For IMEX methods: split the diffusion operator
% (Here we use implicit treatment for the entire operator)
F = [];                  % No nonlinear term
S = [];                  % No source term

% Exact solution for comparison
exact_solution = @(t) exp(-pi^2 * alpha * t) * sin(pi * x);

%% Setup integrators
% Different time step sizes based on stability requirements
dt_implicit = 0.001;     % Larger step for implicit methods
dt_explicit = 0.0001;    % Much smaller for explicit (CFL condition)
dt_imex = 0.001;         % Similar to implicit for IMEX

% Stability analysis for explicit methods
% CFL condition: dt ≤ dx²/(2α) for stability
dt_cfl = dx^2 / (2 * alpha);
fprintf('CFL stability limit for explicit methods: dt ≤ %.6f\n', dt_cfl);
if dt_explicit > dt_cfl
    fprintf('Warning: Explicit time step (%.6f) exceeds CFL limit!\n', dt_explicit);
    dt_explicit = 0.8 * dt_cfl;  % Use 80% of CFL limit for safety
    fprintf('Adjusted explicit time step to: %.6f\n', dt_explicit);
end

% Create integrators
implicit_integrators = {
    approx.odeint.BeIntegrator(final_time)
    approx.odeint.Bdf2Integrator(final_time)
    approx.odeint.Bdf3Integrator(final_time)
    approx.odeint.Sdirk2Integrator(final_time)
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
fprintf('\nSolving 1D heat equation with multiple integrator types...\n');
fprintf('Grid points: %d, dx = %.4f\n', n, dx);
fprintf('Implicit methods dt: %.5f\n', dt_implicit);
fprintf('Explicit methods dt: %.5f\n', dt_explicit);
fprintf('IMEX methods dt: %.5f\n\n', dt_imex);

% Storage for solutions at final time
U_solutions = zeros(n_methods, n);
computational_costs = zeros(n_methods, 1);

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
    
    % Initialize solution history
    integrator.update(u0);
    
    % Time stepping
    U_current = u0;
    tic;
    for k = 1:n_steps
        integrator.setTimeStep(dt=dt);
        
        % Set coefficients for BDF methods
        if isa(integrator, 'approx.odeint.BdfIntegrator')
            integrator.setCoefficients();
        end
        
        % Call appropriate step method based on integrator type
        if strcmp(integrator_types{i}, 'Implicit')
            U_new = integrator.step(L=L, S=S, M=[]);
        elseif strcmp(integrator_types{i}, 'Explicit')
            % For explicit methods, treat diffusion as nonlinear operator
            F_eval = L * integrator.U0{1};
            U_new = integrator.step(L=[], F=F_eval, S=[], M=[]);
        else % IMEX
            % For IMEX methods, treat diffusion implicitly
            U_new = integrator.step(L=L, F=F, S=S, M=[]);
        end
        
        integrator.update(U_new);
        integrator.advance();
        U_current = U_new;
        
        % Progress indicator for explicit methods (many steps)
        if strcmp(integrator_types{i}, 'Explicit') && mod(k, round(n_steps/10)) == 0
            fprintf('  Progress: %d%%\n', round(100*k/n_steps));
        end
    end
    computational_time = toc;
    
    % Store final solution and computational cost
    U_solutions(i, :) = U_current;
    computational_costs(i) = computational_time;
    
    fprintf('  Completed in %.3f seconds (%d steps)\n', computational_time, n_steps);
end

%% Compute exact solution at final time
U_exact = exact_solution(final_time);

%% Plot results
figure('Position', [100, 100, 1500, 800]);

% Define colors for different types
colors_implicit = lines(length(implicit_integrators));
colors_explicit = autumn(length(explicit_integrators));
colors_imex = spring(length(imex_integrators));

% Solution comparison at final time
subplot(2, 3, 1);
plot(x, U_exact, 'k-', 'LineWidth', 3, 'DisplayName', 'Exact');
hold on;

color_idx = [1, 1, 1]; % Index for each type
for i = 1:n_methods
    if strcmp(integrator_types{i}, 'Implicit')
        color = colors_implicit(color_idx(1), :);
        style = 'o';
        color_idx(1) = color_idx(1) + 1;
    elseif strcmp(integrator_types{i}, 'Explicit')
        color = colors_explicit(color_idx(2), :);
        style = 's';
        color_idx(2) = color_idx(2) + 1;
    else % IMEX
        color = colors_imex(color_idx(3), :);
        style = '^';
        color_idx(3) = color_idx(3) + 1;
    end
    
    plot(x, U_solutions(i, :), style, 'Color', color, 'MarkerSize', 4, ...
        'DisplayName', [integrator_names{i} ' (' integrator_types{i} ')']);
end

xlabel('x');
ylabel('u(x,t)');
title(sprintf('Heat Equation Solution at t = %.3f', final_time));
legend('Location', 'northeast', 'FontSize', 8);
grid on;

% Error comparison
subplot(2, 3, 2);
color_idx = [1, 1, 1]; % Reset index
for i = 1:n_methods
    error = abs(U_solutions(i, :)' - U_exact);
    
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
    
    semilogy(x, error, style, 'Color', color, 'LineWidth', 1.5, ...
        'DisplayName', [integrator_names{i} ' (' integrator_types{i} ')']);
    hold on;
end

xlabel('x');
ylabel('Absolute Error');
title('Error Distribution');
legend('Location', 'best', 'FontSize', 8);
grid on;

% Solution evolution for selected methods (implicit and explicit comparison)
subplot(2, 3, 3);
t_evolution = [0, final_time/4, final_time/2, final_time];
colors_evolution = {'k-', 'b-', 'r-', 'g-'};

for t_idx = 1:length(t_evolution)
    U_t = exact_solution(t_evolution(t_idx));
    plot(x, U_t, colors_evolution{t_idx}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('t = %.3f', t_evolution(t_idx)));
    hold on;
end

xlabel('x');
ylabel('u(x,t)');
title('Solution Evolution (Exact)');
legend('Location', 'northeast');
grid on;

% Computational cost comparison
subplot(2, 3, 4);
color_idx = [1, 1, 1]; % Reset index
for i = 1:n_methods
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
    
    semilogy(i, computational_costs(i), marker, 'Color', color, ...
        'MarkerSize', 10, 'MarkerFaceColor', color, 'LineWidth', 1.5);
    hold on;
end

set(gca, 'XTick', 1:n_methods);
set(gca, 'XTickLabel', integrator_names);
xtickangle(45);
ylabel('Computational Time (seconds)');
title('Computational Cost Comparison');
grid on;

% Work-precision diagram
subplot(2, 3, 5);
max_errors = zeros(n_methods, 1);
for i = 1:n_methods
    error = abs(U_solutions(i, :)' - U_exact);
    max_errors(i) = max(error);
end

color_idx = [1, 1, 1]; % Reset index
for i = 1:n_methods
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
    
    loglog(computational_costs(i), max_errors(i), marker, 'Color', color, ...
        'MarkerSize', 8, 'MarkerFaceColor', color, 'LineWidth', 1.5, ...
        'DisplayName', [integrator_names{i} ' (' integrator_types{i} ')']);
    hold on;
end

xlabel('Computational Time (seconds)');
ylabel('Maximum Error');
title('Work-Precision Diagram');
legend('Location', 'best', 'FontSize', 8);
grid on;

% L2 error comparison
subplot(2, 3, 6);
l2_errors = zeros(n_methods, 1);
color_idx = [1, 1, 1]; % Reset index

for i = 1:n_methods
    error = abs(U_solutions(i, :)' - U_exact);
    l2_errors(i) = sqrt(dx * sum(error.^2));
    
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
    
    semilogy(i, l2_errors(i), marker, 'Color', color, ...
        'MarkerSize', 10, 'MarkerFaceColor', color, 'LineWidth', 1.5);
    hold on;
end

set(gca, 'XTick', 1:n_methods);
set(gca, 'XTickLabel', integrator_names);
xtickangle(45);
ylabel('L2 Error');
title('L2 Error Comparison');