%% ex03_convergence_analysis.m
% Convergence order analysis for all integrator types
%
% Tests the theoretical convergence orders of different methods
% using the Prothero-Robinson equation with known exact solution.
% Includes explicit, implicit, IMEX, and splitting integrators.

clear; close all; clc;

%% Problem Setup
lambda = -50;            % Moderate stiffness
phi = @(t) exp(-t);
phi_prime = @(t) -exp(-t);
u0 = 1;
final_time = 2.0;

% For implicit and IMEX methods
L = lambda;
S = @(u, t) -lambda * phi(t) + phi_prime(t);

% For explicit methods - combined into nonlinear term
F = @(u, t) lambda * (u - phi(t)) + phi_prime(t);

% Exact solution
exact_solution = @(t) exp(-t);

%% Setup integrators
% Implicit integrators
implicit_integrators = {
    approx.odeint.BeIntegrator(final_time)
    approx.odeint.Bdf2Integrator(final_time)
    approx.odeint.Bdf3Integrator(final_time)
    approx.odeint.ImpIntegrator(final_time)
    approx.odeint.Sdirk2Integrator(final_time)
    approx.odeint.Sdirk3Integrator(final_time)
    approx.odeint.Sdirk4Integrator(final_time)
    approx.odeint.Esdirk3Integrator(final_time)
};

% Explicit integrators
explicit_integrators = {
    approx.odeint.FeIntegrator(final_time)
    approx.odeint.EmpIntegrator(final_time)
    approx.odeint.HeunIntegrator(final_time)
    approx.odeint.Ssprk2Integrator(final_time)
    approx.odeint.Ssprk3Integrator(final_time)
    approx.odeint.Exrk4Integrator(final_time)
};

% IMEX integrators
imex_integrators = {
    approx.odeint.Ars111Integrator(final_time)
    approx.odeint.Ars222Integrator(final_time)
    approx.odeint.Ars443Integrator(final_time)
};

% Splitting integrators (test with same problem structure)
splitting_integrators = {
    approx.odeint.LieTrotterIntegrator(factors={...
        approx.odeint.BeIntegrator(final_time), ...
        approx.odeint.FeIntegrator(final_time)})
    approx.odeint.StrangIntegrator(factors={...
        approx.odeint.Sdirk2Integrator(final_time), ...
        approx.odeint.HeunIntegrator(final_time)})
};

% Combine all integrators
all_integrators = [implicit_integrators; explicit_integrators; imex_integrators; splitting_integrators];
n_methods = length(all_integrators);

% Extract integrator names and types
integrator_names = cell(n_methods, 1);
integrator_types = cell(n_methods, 1);

for i = 1:n_methods
    class_name = class(all_integrators{i});
    name_parts = strsplit(class_name, '.');
    base_name = name_parts{end};
    
    % Determine type
    if i <= length(implicit_integrators)
        integrator_types{i} = 'Implicit';
        integrator_names{i} = strrep(base_name, 'Integrator', '');
    elseif i <= length(implicit_integrators) + length(explicit_integrators)
        integrator_types{i} = 'Explicit';
        integrator_names{i} = strrep(base_name, 'Integrator', '');
    elseif i <= length(implicit_integrators) + length(explicit_integrators) + length(imex_integrators)
        integrator_types{i} = 'IMEX';
        integrator_names{i} = strrep(base_name, 'Integrator', '');
    else % Splitting methods
        integrator_types{i} = 'Splitting';
        if contains(base_name, 'LieTrotter')
            integrator_names{i} = 'LieTrotter';
        else
            integrator_names{i} = 'Strang';
        end
    end
end

%% Step sizes for convergence study
% Different ranges for different method types due to stability constraints
dt_implicit = [0.1, 0.05, 0.025, 0.0125, 0.00625];
dt_explicit = [0.01, 0.005, 0.0025, 0.00125, 0.000625];  % Smaller due to stability
dt_imex = [0.05, 0.025, 0.0125, 0.00625, 0.003125];     % Intermediate
dt_splitting = [0.05, 0.025, 0.0125, 0.00625, 0.003125]; % Similar to IMEX

%% Convergence Study
fprintf('=== CONVERGENCE ORDER ANALYSIS ===\n\n');
fprintf('Note: Different step size ranges used due to stability constraints:\n');
fprintf('  Implicit: dt ∈ [%.5f, %.2f]\n', min(dt_implicit), max(dt_implicit));
fprintf('  Explicit: dt ∈ [%.6f, %.2f]\n', min(dt_explicit), max(dt_explicit));
fprintf('  IMEX:     dt ∈ [%.6f, %.2f]\n', min(dt_imex), max(dt_imex));
fprintf('  Splitting: dt ∈ [%.6f, %.2f]\n\n', min(dt_splitting), max(dt_splitting));

% Storage for errors
errors = cell(n_methods, 1);
dt_values_used = cell(n_methods, 1);

for i = 1:n_methods
    fprintf('Testing %s (%s, theoretical order %d)...\n', ...
        integrator_names{i}, integrator_types{i}, all_integrators{i}.Order);
    
    integrator = all_integrators{i};
    
    % Choose appropriate step sizes based on method type
    if strcmp(integrator_types{i}, 'Explicit')
        dt_values = dt_explicit;
    elseif strcmp(integrator_types{i}, 'IMEX')
        dt_values = dt_imex;
    elseif strcmp(integrator_types{i}, 'Splitting')
        dt_values = dt_splitting;
    else
        dt_values = dt_implicit;
    end
    
    n_dt = length(dt_values);
    dt_values_used{i} = dt_values;
    errors{i} = zeros(1, n_dt);

    for j = 1:n_dt
        dt = dt_values(j);
        
        integrator.reset();
        
        % Initialize
        integrator.update(u0);
        
        % Time stepping
        n_steps = ceil(final_time / dt);
        U_current = u0;
        
        try
            for k = 1:n_steps
                % Set time step
                if strcmp(integrator_types{i}, 'Splitting')
                    integrator.setTimeStep(dt=dt);
                else
                    integrator.setTimeStep(dt=dt);
                end
                
                % Set coefficients for BDF methods
                if isa(integrator, 'approx.odeint.BdfIntegrator')
                    integrator.setCoefficients();
                end
                
                % Call appropriate step method based on integrator type
                if strcmp(integrator_types{i}, 'Implicit')
                    U_new = integrator.step(L=L, S=S, M=[]);
                elseif strcmp(integrator_types{i}, 'Explicit')
                    % For explicit methods, treat everything as nonlinear term
                    U_new = integrator.step(L=[], F=F, S=[], M=[]);
                elseif strcmp(integrator_types{i}, 'IMEX')
                    U_new = integrator.step(L=L, F=[], S=S, M=[]);
                else % Splitting
                    % Split the problem: linear part + source part
                    L_split = {L, []};
                    F_split = {[], []};
                    S_split = {[], S};
                    M_split = {[], []};
                    U_new = integrator.step(L=L_split, F=F_split, S=S_split, M=M_split);
                end
                
                integrator.update(U_new);
                integrator.advance();
                U_current = U_new;
            end
            
            % Compute error
            exact_final = exact_solution(final_time);
            errors{i}(j) = abs(U_current - exact_final);
            
        catch ME
            fprintf('    Error at dt=%.6f: %s\n', dt, ME.message);
            errors{i}(j) = NaN;
        end
    end
    
    % Compute observed convergence order
    log_dt = log(dt_values);
    log_err = log(errors{i});
    
    % Use last 4 points for slope calculation (more stable)
    valid_indices = ~isinf(log_err) & ~isnan(log_err) & errors{i} > 0;
    if sum(valid_indices) >= 3
        valid_log_dt = log_dt(valid_indices);
        valid_log_err = log_err(valid_indices);
        
        if length(valid_log_dt) >= 4
            last_indices = (length(valid_log_dt)-3):length(valid_log_dt);
        else
            last_indices = 1:length(valid_log_dt);
        end
        
        p = polyfit(valid_log_dt(last_indices), valid_log_err(last_indices), 1);
        observed_order = -p(1);
    else
        observed_order = NaN;
    end
    
    if ~isnan(observed_order)
        fprintf('  Observed order: %.2f (theoretical: %d)\n', ...
            observed_order, all_integrators{i}.Order);
    else
        fprintf('  Observed order: N/A (insufficient valid data)\n');
    end
end

%% Plot convergence results
figure('Position', [100, 100, 1400, 1000]);

% Define colors for different types
colors_implicit = lines(length(implicit_integrators));
colors_explicit = autumn(length(explicit_integrators));
colors_imex = spring(length(imex_integrators));
colors_splitting = winter(length(splitting_integrators));

% Main convergence plot
subplot(2, 3, 1);
color_idx = [1, 1, 1, 1]; % Index for each type

for i = 1:n_methods
    if strcmp(integrator_types{i}, 'Implicit')
        color = colors_implicit(color_idx(1), :);
        marker = 'o';
        style = '-';
        color_idx(1) = color_idx(1) + 1;
    elseif strcmp(integrator_types{i}, 'Explicit')
        color = colors_explicit(color_idx(2), :);
        marker = 's';
        style = '--';
        color_idx(2) = color_idx(2) + 1;
    elseif strcmp(integrator_types{i}, 'IMEX')
        color = colors_imex(color_idx(3), :);
        marker = '^';
        style = '-.';
        color_idx(3) = color_idx(3) + 1;
    else % Splitting
        color = colors_splitting(color_idx(4), :);
        marker = 'd';
        style = ':';
        color_idx(4) = color_idx(4) + 1;
    end
    
    dt_vals = dt_values_used{i};
    err_vals = errors{i};
    valid_indices = err_vals > 0 & ~isnan(err_vals) & ~isinf(err_vals);
    
    if sum(valid_indices) > 0
        loglog(dt_vals(valid_indices), err_vals(valid_indices), ...
            [marker style], 'Color', color, 'LineWidth', 2, 'MarkerSize', 6, ...
            'DisplayName', [integrator_names{i} ' (' integrator_types{i} ')']);
        hold on;
    end
end

% Add reference lines for theoretical orders
dt_ref = logspace(log10(0.001), log10(0.1), 50);
for order = 1:4
    ref_slope = 0.01 * dt_ref.^order;
    loglog(dt_ref, ref_slope, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
    text(dt_ref(25), ref_slope(25)*1.5, sprintf('%d', order), ...
        'FontSize', 8, 'HorizontalAlignment', 'center');
end

xlabel('Time Step Size (dt)');
ylabel('Global Error at t = 2');
title('Convergence Analysis: All Method Types');
legend('Location', 'southeast', 'FontSize', 8);
grid on;
set(gca, 'FontSize', 10);

% Separate plots for each method type
method_types = {'Implicit', 'Explicit', 'IMEX', 'Splitting'};
colors_by_type = {colors_implicit, colors_explicit, colors_imex, colors_splitting};

for type_idx = 1:4
    subplot(2, 3, type_idx + 1);
    
    type_name = method_types{type_idx};
    type_indices = find(strcmp(integrator_types, type_name));
    
    if isempty(type_indices)
        continue;
    end
    
    colors = colors_by_type{type_idx};
    
    color_idx = 1;
    for i = type_indices'
        dt_vals = dt_values_used{i};
        err_vals = errors{i};
        valid_indices = err_vals > 0 & ~isnan(err_vals) & ~isinf(err_vals);
        
        if sum(valid_indices) > 0
            % Different markers for different method types
            if strcmp(type_name, 'Splitting')
                if contains(integrator_names{i}, 'LieTrotter')
                    marker_style = 'd--';
                else
                    marker_style = 'd-';
                end
            else
                marker_style = 'o-';
            end
            
            loglog(dt_vals(valid_indices), err_vals(valid_indices), ...
                marker_style, 'Color', colors(color_idx, :), 'LineWidth', 2, 'MarkerSize', 8, ...
                'DisplayName', integrator_names{i});
            hold on;
        end
        color_idx = color_idx + 1;
    end
    
    % Add reference lines
    if strcmp(type_name, 'Explicit')
        dt_ref = dt_explicit(1:3);
    elseif strcmp(type_name, 'Splitting')
        dt_ref = dt_splitting(1:3);
    elseif strcmp(type_name, 'IMEX')
        dt_ref = dt_imex(1:3);
    else
        dt_ref = dt_implicit(1:3);
    end
    
    max_order = 4;
    if strcmp(type_name, 'Splitting')
        max_order = 2; % Focus on orders 1 and 2 for splitting
    end
    
    for order = 1:max_order
        ref_slope = 0.01 * dt_ref.^order;
        loglog(dt_ref, ref_slope, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    
    xlabel('Time Step Size (dt)');
    ylabel('Global Error');
    title([type_name ' Methods']);
    legend('Location', 'southeast', 'FontSize', 8);
    grid on;
    set(gca, 'FontSize', 10);
end

% Additional subplot for splitting method detailed comparison
subplot(2, 3, 6);
splitting_indices = find(strcmp(integrator_types, 'Splitting'));

if ~isempty(splitting_indices)
    color_idx = 1;
    for i = splitting_indices'
        dt_vals = dt_values_used{i};
        err_vals = errors{i};
        valid_indices = err_vals > 0 & ~isnan(err_vals) & ~isinf(err_vals);
        
        if sum(valid_indices) > 0
            % Different line styles for Lie-Trotter vs Strang
            if contains(integrator_names{i}, 'LieTrotter')
                style = 'd--';
                lwidth = 2;
            else
                style = 'd-';
                lwidth = 3;
            end
            
            loglog(dt_vals(valid_indices), err_vals(valid_indices), ...
                style, 'Color', colors_splitting(color_idx, :), 'LineWidth', lwidth, 'MarkerSize', 10, ...
                'DisplayName', integrator_names{i});
            hold on;
        end
        color_idx = color_idx + 1;
    end
    
    % Reference lines for order 1 and 2
    dt_ref = dt_splitting(1:4);
    ref1 = 0.01 * dt_ref.^1;
    ref2 = 0.01 * dt_ref.^2;
    loglog(dt_ref, ref1, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
    loglog(dt_ref, ref2, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
    text(dt_ref(2), ref1(2)*2, '1', 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(dt_ref(2), ref2(2)*0.5, '2', 'FontSize', 10, 'HorizontalAlignment', 'center');
    
    xlabel('Time Step Size (dt)');
    ylabel('Global Error');
    title('Splitting Methods Detail');
    legend('Location', 'southeast', 'FontSize', 10);
    grid on;
    set(gca, 'FontSize', 10);
else
    text(0.5, 0.5, 'No Splitting Methods', 'HorizontalAlignment', 'center');
    title('Splitting Methods Detail');
end

%% Detailed analysis table by method type
fprintf('\n=== DETAILED CONVERGENCE ANALYSIS BY TYPE ===\n');

for type_idx = 1:4
    if type_idx <= 3
        type_name = method_types{type_idx};
    else
        type_name = 'Splitting';
    end
    
    type_indices = find(strcmp(integrator_types, type_name));
    
    if isempty(type_indices)
        continue;
    end
    
    fprintf('\n--- %s METHODS ---\n', upper(type_name));
    
    % Find common step size range for this type
    if strcmp(type_name, 'Explicit')
        common_dt = dt_explicit;
    elseif strcmp(type_name, 'Splitting')
        common_dt = dt_splitting;
    elseif strcmp(type_name, 'IMEX')
        common_dt = dt_imex;
    else
        common_dt = dt_implicit;
    end
    
    fprintf('%-15s', 'dt');
    for i = type_indices'
        fprintf(' | %-12s', integrator_names{i});
    end
    fprintf('\n');
    fprintf('%s\n', repmat('-', 1, 15 + 15 * length(type_indices)));
    
    for j = 1:length(common_dt)
        fprintf('%-15.6f', common_dt(j));
        for i = type_indices'
            err_val = errors{i}(j);
            if err_val > 0 && ~isnan(err_val) && ~isinf(err_val)
                fprintf(' | %12.2e', err_val);
            else
                fprintf(' | %12s', 'N/A');
            end
        end
        fprintf('\n');
    end
    
    % Convergence rates for this type
    fprintf('\nConvergence rates (consecutive step ratios):\n');
    fprintf('%-15s', 'dt ratio');
    for i = type_indices'
        fprintf(' | %-12s', integrator_names{i});
    end
    fprintf('\n');
    fprintf('%s\n', repmat('-', 1, 15 + 15 * length(type_indices)));
    
    for j = 2:length(common_dt)
        dt_ratio = common_dt(j-1) / common_dt(j);
        fprintf('%-15.1f', dt_ratio);
        
        for i = type_indices'
            err_prev = errors{i}(j-1);
            err_curr = errors{i}(j);
            
            if err_prev > 0 && err_curr > 0 && ~isnan(err_prev) && ~isnan(err_curr)
                rate = log(err_prev / err_curr) / log(dt_ratio);
                fprintf(' | %12.2f', rate);
            else
                fprintf(' | %12s', 'N/A');
            end
        end
        fprintf('\n');
    end
end

%% Summary comparison
fprintf('\n=== SUMMARY COMPARISON ===\n');
fprintf('Best performers by type (smallest error at finest resolution):\n\n');

for type_idx = 1:4
    if type_idx <= 3
        type_name = method_types{type_idx};
    else
        type_name = 'Splitting';
    end
    
    type_indices = find(strcmp(integrator_types, type_name));
    
    if isempty(type_indices)
        continue;
    end
    
    % Find best method in this category
    min_error = inf;
    best_idx = 0;
    
    for i = type_indices'
        final_error = errors{i}(end);
        if final_error > 0 && final_error < min_error
            min_error = final_error;
            best_idx = i;
        end
    end
    
    if best_idx > 0
        fprintf('%s methods: %s (error: %.2e)\n', ...
            type_name, integrator_names{best_idx}, min_error);
    else
        fprintf('%s methods: No valid results\n', type_name);
    end
end

fprintf('\nNote: Different step size ranges were used due to stability constraints.\n');
fprintf('Explicit methods required much smaller time steps for the stiff problem.\n');
fprintf('Splitting methods use operator splitting: [Linear part] + [Source part].\n');