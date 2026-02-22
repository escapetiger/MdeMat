%% ex06_custom_step_functions.m
% Demonstrate the new name-value interface and custom step function capability
%
% Problem: du/dt = λu + S(t), u(0) = 1
% where λ = -10 and S(t) = sin(t)
% Exact solution can be computed analytically
%
% This example shows:
% 1. Using the new name-value syntax for integrator step methods
% 2. Implementing custom step functions for specialized physics
% 3. Factor-specific customization in splitting methods

clear; close all; clc;

%% Problem Setup
lambda = -10;            % Linear coefficient
source_fcn = @(u, t) sin(t);  % Time-dependent source
u0 = 1;                  % Initial condition
final_time = 2.0;
dt = 0.01;

% For splitting: decompose into λu and source term
L_split = {lambda, []};
F_split = {[], []};
S_split = {[], source_fcn};
M_split = {[], []};

fprintf('=== Demonstrating New ODE Integrator Interface ===\n\n');

%% Example 1: Basic usage with new name-value syntax
fprintf('1. Basic usage with new name-value syntax:\n');

% Create a simple implicit integrator
integrator1 = approx.odeint.BeIntegrator(final_time);
integrator1.reset();
integrator1.add(u0);
integrator1.setTimeStep(dt=dt);

% Use new interface: explicit parameter names
U_new = integrator1.step(L=lambda, F=[], S=source_fcn, M=[]);
fprintf('   BE step result: %.6f\n', U_new);

%% Example 2: Custom step function for specialized physics
fprintf('\n2. Custom step function example:\n');

% Define a custom step function that logs its execution
customStepFcn = @(integrator, options) customPhysicsStep(integrator, options);

% Reset and use custom function
integrator1.reset();
integrator1.add(u0);

% Step with custom function
U_custom = integrator1.step(L=lambda, S=source_fcn, M=[], fn=customStepFcn, ...
                           context=struct('problemName', 'CustomPhysics'));
fprintf('   Custom step result: %.6f\n', U_custom);

%% Example 3: Splitting methods with factor-specific custom functions
fprintf('\n3. Splitting with custom factor functions:\n');

% Create Strang integrator with BE factors
factor1 = approx.odeint.BeIntegrator(final_time);
factor2 = approx.odeint.BeIntegrator(final_time);
strang = approx.odeint.StrangIntegrator(factors={factor1, factor2});

strang.reset();
strang.add(u0);
strang.setTimeStep(dt=dt);

% Custom functions for each factor
linearStepFcn = @(integrator, options) customLinearStep(integrator, options, 'Linear Factor');
sourceStepFcn = @(integrator, options) customSourceStep(integrator, options, 'Source Factor');

% Execute with factor-specific custom functions
U_strang = strang.step(L=L_split, F=F_split, S=S_split, M=M_split, ...
                      fn={linearStepFcn, sourceStepFcn});
fprintf('   Strang with custom factors: %.6f\n', U_strang);

%% Example 4: Advanced usage with time override
fprintf('\n4. Advanced usage with custom time step:\n');

% Create SDIRK integrator
integrator2 = approx.odeint.Sdirk2Integrator(final_time);
integrator2.reset();
integrator2.add(u0);

% Use different time step than integrator's timeline
custom_dt = 0.005;
U_sdirk = integrator2.step(L=lambda, S=source_fcn, M=[], dt=custom_dt);
fprintf('   SDIRK2 with dt=%.3f: %.6f\n', custom_dt, U_sdirk);

%% Example 5: Convergence study with multiple integrators
fprintf('\n5. Convergence study using new interface:\n');

integrators = {
    approx.odeint.BeIntegrator(final_time)
    approx.odeint.Sdirk2Integrator(final_time)
    approx.odeint.FeIntegrator(final_time)
};

integrator_names = {'BE', 'SDIRK2', 'FE'};
dt_values = [0.1, 0.05, 0.01, 0.005];

fprintf('   Method    dt=0.1    dt=0.05   dt=0.01   dt=0.005\n');
fprintf('   ------    ------    -------   -------   --------\n');

for i = 1:length(integrators)
    integrator = integrators{i};
    results = zeros(1, length(dt_values));

    for j = 1:length(dt_values)
        dt_test = dt_values(j);
        integrator.reset();
        integrator.add(u0);

        % Single step with custom dt
        if strcmp(integrator_names{i}, 'FE')
            % Explicit method
            U_result = integrator.step(L=[], F=lambda*u0 + source_fcn(u0, 0), S=[], M=[], dt=dt_test);
        else
            % Implicit method
            U_result = integrator.step(L=lambda, S=source_fcn, M=[], dt=dt_test);
        end
        results(j) = U_result;
    end

    fprintf('   %-8s', integrator_names{i});
    for j = 1:length(results)
        fprintf(' %8.5f', results(j));
    end
    fprintf('\n');
end

fprintf('\n=== All examples completed ===\n');

%% Helper Functions

function U = customPhysicsStep(integrator, options)
    % CUSTOMPHYSICSSTEP Example custom step function with logging.

    fprintf('   → Custom physics step called for: %s\n', options.context.problemName);
    fprintf('   → Integrator type: %s\n', class(integrator));

    % Use standard BE implementation with logging
    integrator.setLhs(L=options.L, M=options.M);
    integrator.setRhs(S=options.S, M=options.M);
    U = integrator.Solver.solve();

    fprintf('   → Custom step completed: U = %.6f\n', U);
end

function U = customLinearStep(integrator, options, factorName)
    % CUSTOMLINEARSTEP Custom step for linear factor with specialized handling.

    fprintf('   → %s: Custom linear step (dt=%.3f)\n', factorName, options.dt);

    % Example: Could implement optimized solver for linear systems
    integrator.setLhs(L=options.L, M=options.M, dt=options.dt);
    integrator.setRhs(S=options.S, M=options.M, dt=options.dt);
    U = integrator.Solver.solve();
end

function U = customSourceStep(integrator, options, factorName)
    % CUSTOMSOURCESTEP Custom step for source factor with analytical treatment.

    fprintf('   → %s: Custom source step (dt=%.3f)\n', factorName, options.dt);

    % Example: Analytical integration of source term
    % For S(t) = sin(t), could use exact integration
    integrator.setLhs(L=options.L, M=options.M, dt=options.dt);
    integrator.setRhs(S=options.S, M=options.M, dt=options.dt);
    U = integrator.Solver.solve();
end