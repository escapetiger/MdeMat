function dev06_convergence_test()
% DEV06_CONVERGENCE_TEST Test spatial convergence rates for orthogonal decomposition.
%
%   Runs simulations at multiple spatial resolutions and computes
%   convergence rates by comparing with exact solution.

    clc; close all;

    fprintf('=== Spatial Convergence Test ===\n\n');

    % Test configurations
    grid_sizes = [21, 31, 41, 51, 71];  % Different spatial resolutions
    fdOrders = [1, 2];  % Test both first and second order schemes

    % Run convergence tests for each order
    for orderIdx = 1:length(fdOrders)
        fdOrder = fdOrders(orderIdx);
        fprintf('\n========== Testing Order %d Scheme ==========\n', fdOrder);

        results = runConvergenceStudy(grid_sizes, fdOrder);

        % Compute and display convergence rates
        analyzeConvergence(results, fdOrder);

        % Plot convergence results
        plotConvergenceRates(results, fdOrder);
    end

    fprintf('\nConvergence test complete!\n');

end

function results = runConvergenceStudy(grid_sizes, fdOrder)
% RUNCONVERGENCESTUDY Run simulations at different resolutions.

    n_grids = length(grid_sizes);
    results = struct();
    results.grid_sizes = grid_sizes;
    results.h_values = zeros(n_grids, 1);
    results.l2_errors = zeros(n_grids, 1);
    results.linf_errors = zeros(n_grids, 1);
    results.mass_errors = zeros(n_grids, 1);

    fprintf('\nRunning simulations at %d different resolutions...\n', n_grids);

    for idx = 1:n_grids
        Nx = grid_sizes(idx);
        fprintf('\n--- Resolution %d/%d: %d×%d grid ---\n', idx, n_grids, Nx, Nx);

        % Setup parameters for this resolution
        params = setupConvergenceParameters(Nx, fdOrder);

        % Run simulation
        [l2_err, linf_err, mass_err, h] = runSingleSimulation(params);

        % Store results
        results.h_values(idx) = h;
        results.l2_errors(idx) = l2_err;
        results.linf_errors(idx) = linf_err;
        results.mass_errors(idx) = mass_err;

        fprintf('  Grid spacing h = %.4f\n', h);
        fprintf('  L2 error   = %.6e\n', l2_err);
        fprintf('  L∞ error   = %.6e\n', linf_err);
        fprintf('  Mass error = %.6e\n', mass_err);
    end

end

function params = setupConvergenceParameters(Nx, fdOrder)
% SETUPCONVERGENCEPARAMETERS Setup parameters for convergence test.

    % Spatial domain
    params.x_min = -1.0;
    params.x_max = 1.0;
    params.y_min = -1.0;
    params.y_max = 1.0;
    params.Nx = Nx;
    params.Ny = Nx;  % Square grid

    % Velocity space
    params.Nv = 16;  % Fixed velocity discretization

    % Decomposition parameters
    params.Nu = 3;  % {1, cos(θ), sin(θ)}
    params.Nc = 0;

    % Wendland parameter
    params.epsilon = 2.0;

    % Spatial discretization
    params.fdOrder = fdOrder;

    % Time integration - use small time step to minimize temporal error
    params.T_final = 0.2;  % Shorter time to focus on spatial error
    h = (params.x_max - params.x_min) / (Nx - 1);
    params.dt = 0.1 * h;  % CFL-based time step
    params.Nt = round(params.T_final / params.dt);

    % Initial condition
    params.sigma = 0.1;  % Gaussian width
    params.x0 = 0.0;
    params.y0 = 0.0;

    % Output
    params.n_snapshots = 1;  % Only save final state

end

function [l2_error, linf_error, mass_error, h] = runSingleSimulation(params)
% RUNSINGLESI MULATION Run one simulation and compute errors.

    % Setup grids
    x_grid = linspace(params.x_min, params.x_max, params.Nx);
    y_grid = linspace(params.y_min, params.y_max, params.Ny);
    [X, Y] = meshgrid(x_grid, y_grid);
    theta_grid = linspace(0, 2*pi, params.Nv + 1);
    theta_grid = theta_grid(1:end-1);

    h = x_grid(2) - x_grid(1);

    % Build basis functions
    basis = buildBasisFunctions(theta_grid, params);

    % Assemble flux matrices
    [A1, A2] = assembleFluxMatrices(basis, params);

    % Setup initial condition
    W0 = setupInitialCondition(X, Y, theta_grid, basis, params);

    % Time integration
    [W_all, times] = timeIntegration(W0, A1, A2, x_grid, y_grid, params, basis);

    % Compute exact solution at final time
    [f_exact, ~] = computeExactSolution(X, Y, theta_grid, times(end), params);

    % Reconstruct solution from W
    W_final = W_all(:, :, :, end);
    f_decomp = reconstructDistribution(W_final, basis, params);

    % Compute errors
    [l2_error, linf_error, mass_error] = computeErrors(f_decomp, f_exact, params);

end

function f = reconstructDistribution(W, basis, params)
% RECONSTRUCTDISTRIBUTION Reconstruct f from state vector W.

    Nx = params.Nx;
    Ny = params.Ny;
    Nv = params.Nv;
    Nu = basis.Nu;
    Nc = basis.Nc;

    f = zeros(Ny, Nx, Nv);

    for i = 1:Ny
        for j = 1:Nx
            W_local = squeeze(W(i, j, :));
            U = W_local(1:Nu);
            C = W_local(Nu+1:Nu+Nc);
            R = W_local(Nu+Nc+1:end);

            f(i, j, :) = basis.E * U + basis.B * C + R;
        end
    end

end

function [l2_error, linf_error, mass_error] = computeErrors(f_computed, f_exact, params)
% COMPUTEERRORS Compute L2, Linf, and mass conservation errors.

    Nx = params.Nx;
    Ny = params.Ny;
    Nv = params.Nv;

    dx = (params.x_max - params.x_min) / (Nx - 1);
    dy = (params.y_max - params.y_min) / (Ny - 1);
    dv = 2*pi / Nv;

    % Compute densities
    rho_computed = sum(f_computed, 3) * dv;
    rho_exact = sum(f_exact, 3) * dv;

    % L2 error in density
    diff = rho_computed - rho_exact;
    l2_error = sqrt(trapz2d(diff.^2, dx, dy));

    % L∞ error in density
    linf_error = max(abs(diff(:)));

    % Mass conservation error
    mass_computed = trapz2d(rho_computed, dx, dy);
    mass_exact = trapz2d(rho_exact, dx, dy);
    mass_error = abs(mass_computed - mass_exact) / mass_exact;

end

function analyzeConvergence(results, fdOrder)
% ANALYZECONVERGENCE Compute and display convergence rates.

    fprintf('\n========== Convergence Analysis (Order %d) ==========\n', fdOrder);
    fprintf('\n  h          L2 Error    L∞ Error    Mass Error\n');
    fprintf('  -----------------------------------------------\n');

    for i = 1:length(results.h_values)
        fprintf('  %.5f   %.6e  %.6e  %.6e\n', ...
                results.h_values(i), results.l2_errors(i), ...
                results.linf_errors(i), results.mass_errors(i));
    end

    % Compute convergence rates
    n = length(results.h_values);
    if n >= 2
        fprintf('\n  Convergence Rates (between consecutive grids):\n');
        fprintf('  -----------------------------------------------\n');

        for i = 2:n
            h_ratio = results.h_values(i-1) / results.h_values(i);

            l2_rate = log(results.l2_errors(i-1) / results.l2_errors(i)) / log(h_ratio);
            linf_rate = log(results.linf_errors(i-1) / results.linf_errors(i)) / log(h_ratio);

            fprintf('  h: %.5f -> %.5f:  L2 rate = %.2f,  L∞ rate = %.2f\n', ...
                    results.h_values(i-1), results.h_values(i), l2_rate, linf_rate);
        end

        % Average convergence rate
        if n >= 3
            avg_l2_rate = mean(log(results.l2_errors(1:end-1) ./ results.l2_errors(2:end)) ./ ...
                              log(results.h_values(1:end-1) ./ results.h_values(2:end)));
            avg_linf_rate = mean(log(results.linf_errors(1:end-1) ./ results.linf_errors(2:end)) ./ ...
                                log(results.h_values(1:end-1) ./ results.h_values(2:end)));

            fprintf('\n  Average convergence rates:\n');
            fprintf('    L2:  %.2f  (expected: %d)\n', avg_l2_rate, fdOrder);
            fprintf('    L∞:  %.2f  (expected: %d)\n', avg_linf_rate, fdOrder);
        end
    end

end

function plotConvergenceRates(results, fdOrder)
% PLOTCONVERGENCERATES Plot convergence results.

    figure;

    % L2 error plot
    subplot(1, 2, 1);
    loglog(results.h_values, results.l2_errors, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;

    % Add reference slopes
    h_ref = results.h_values;
    if fdOrder == 1
        loglog(h_ref, results.l2_errors(1) * (h_ref / h_ref(1)).^1, 'k--', 'LineWidth', 1.5, 'DisplayName', 'O(h)');
    else
        loglog(h_ref, results.l2_errors(1) * (h_ref / h_ref(1)).^2, 'k--', 'LineWidth', 1.5, 'DisplayName', 'O(h²)');
    end

    xlabel('Grid spacing h', 'FontSize', 12);
    ylabel('L2 Error', 'FontSize', 12);
    title(sprintf('L2 Convergence (Order %d FD)', fdOrder), 'FontSize', 14);
    legend('Computed', sprintf('O(h^%d)', fdOrder), 'Location', 'best');
    grid on;

    % L∞ error plot
    subplot(1, 2, 2);
    loglog(results.h_values, results.linf_errors, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;

    % Add reference slopes
    if fdOrder == 1
        loglog(h_ref, results.linf_errors(1) * (h_ref / h_ref(1)).^1, 'k--', 'LineWidth', 1.5, 'DisplayName', 'O(h)');
    else
        loglog(h_ref, results.linf_errors(1) * (h_ref / h_ref(1)).^2, 'k--', 'LineWidth', 1.5, 'DisplayName', 'O(h²)');
    end

    xlabel('Grid spacing h', 'FontSize', 12);
    ylabel('L∞ Error', 'FontSize', 12);
    title(sprintf('L∞ Convergence (Order %d FD)', fdOrder), 'FontSize', 14);
    legend('Computed', sprintf('O(h^%d)', fdOrder), 'Location', 'best');
    grid on;

    sgtitle(sprintf('Spatial Convergence Study (Order %d Finite Difference)', fdOrder), 'FontSize', 16);

end

% ========== Helper Functions (copied from main file) ==========

function basis = buildBasisFunctions(theta_grid, params)
% BUILDBBASISFUNCTIONS Build macro, meso, and micro basis functions.

    Nv = params.Nv;
    Nu = params.Nu;
    Nc = params.Nc;

    w = (2*pi / Nv) * ones(Nv, 1);
    basis.W = diag(w);

    basis.v_x = cos(theta_grid)';
    basis.v_y = sin(theta_grid)';
    basis.theta = theta_grid';

    if Nu > 0
        basis.E = zeros(Nv, Nu);
        basis.E(:, 1) = 1 / sqrt(2*pi);
        if Nu >= 2
            basis.E(:, 2) = sqrt(2/(2*pi)) * cos(theta_grid)';
        end
        if Nu >= 3
            basis.E(:, 3) = sqrt(2/(2*pi)) * sin(theta_grid)';
        end
    else
        basis.E = zeros(Nv, 0);
    end

    if Nc > 0
        basis.B = buildWendlandBasis(theta_grid, params);
    else
        basis.B = zeros(Nv, 0);
    end

    basis.Nu = Nu;
    basis.Nc = Nc;
    basis.Nv = Nv;
    basis.N = Nu + Nc + Nv;

end

function B = buildWendlandBasis(theta_grid, params)
% BUILDWENDLANDBASIS Build orthonormalized Wendland kernel basis.

    Nv = length(theta_grid);
    Nc = params.Nc;
    epsilon = params.epsilon;

    w = (2*pi / Nv) * ones(Nv, 1);
    W = diag(w);

    if params.Nu > 0
        E = zeros(Nv, params.Nu);
        E(:, 1) = 1 / sqrt(2*pi);
        if params.Nu >= 2
            E(:, 2) = sqrt(2/(2*pi)) * cos(theta_grid)';
        end
        if params.Nu >= 3
            E(:, 3) = sqrt(2/(2*pi)) * sin(theta_grid)';
        end
        P_E = E * (E' * W);
    else
        P_E = zeros(Nv, Nv);
    end

    theta_centers = linspace(0, 2*pi, Nc + 1);
    theta_centers = theta_centers(1:end-1);

    B_raw = zeros(Nv, Nc);
    for j = 1:Nc
        for i = 1:Nv
            d_theta = angularDistance(theta_grid(i), theta_centers(j));
            B_raw(i, j) = wendlandC2(d_theta, epsilon);
        end
    end

    B_proj = (eye(Nv) - P_E) * B_raw;

    B = zeros(Nv, Nc);
    for j = 1:Nc
        b_j = B_proj(:, j);
        for k = 1:j-1
            b_j = b_j - (b_j' * W * B(:, k)) * B(:, k);
        end
        norm_bj = sqrt(b_j' * W * b_j);
        if norm_bj > 1e-12
            B(:, j) = b_j / norm_bj;
        else
            B(:, j) = b_j;
        end
    end

end

function d = angularDistance(theta1, theta2)
% ANGULARDISTANCE Compute angular distance on circle.

    d = abs(theta1 - theta2);
    d = min(d, 2*pi - d);

end

function phi = wendlandC2(r, epsilon)
% WENDLANDC2 Wendland C² compactly supported RBF.

    s = epsilon * r;
    if s >= 1
        phi = 0;
    else
        phi = (1 - s)^4 * (4*s + 1);
    end

end

function [A1, A2] = assembleFluxMatrices(basis, params)
% ASSEMBLEFLUXMATRICES Assemble flux matrices.

    E = basis.E;
    B = basis.B;
    W = basis.W;
    Nu = basis.Nu;
    Nc = basis.Nc;
    Nv = basis.Nv;
    N = basis.N;

    V1 = diag(basis.v_x);
    V2 = diag(basis.v_y);

    if Nu > 0
        P_E = E * (E' * W);
    else
        P_E = zeros(Nv, Nv);
    end
    if Nc > 0
        P_B = B * (B' * W);
    else
        P_B = zeros(Nv, Nv);
    end

    if Nu > 0
        if Nu == 1
            L = 0;
        elseif Nu == 2
            L = [0, 1];
        elseif Nu == 3
            L = [0, 1, 1];
        else
            k = ceil((Nu - 1) / 2);
            l = [0, repelem(1:k, 2)];
            L = l(1:Nu);
        end

        max_level = max(L);
        S = ones(1, Nu);
        S(L < max_level) = 0;
        S = diag(S);

        S_P_E = E * S * (E' * W);
    else
        S_P_E = zeros(Nv, Nv);
    end

    A1 = zeros(N, N);
    if Nu > 0
        A1(1:Nu, 1:Nu) = E' * W * V1 * E;
    end
    if Nu > 0 && Nc > 0
        A1(1:Nu, Nu+1:Nu+Nc) = E' * W * V1 * B;
    end
    if Nu > 0
        A1(1:Nu, Nu+Nc+1:N) = S * E' * W * V1;
    end
    if Nc > 0 && Nu > 0
        A1(Nu+1:Nu+Nc, 1:Nu) = B' * W * V1 * E;
    end
    if Nc > 0
        A1(Nu+1:Nu+Nc, Nu+1:Nu+Nc) = B' * W * V1 * B;
    end
    if Nc > 0
        A1(Nu+1:Nu+Nc, Nu+Nc+1:N) = B' * W * V1;
    end
    if Nu > 0
        A1(Nu+Nc+1:N, 1:Nu) = (eye(Nv) - P_B) * (eye(Nv) - P_E) * V1 * E;
    end
    if Nc > 0
        A1(Nu+Nc+1:N, Nu+1:Nu+Nc) = (eye(Nv) - P_B) * (eye(Nv) - P_E) * V1 * B;
    end
    A1(Nu+Nc+1:N, Nu+Nc+1:N) = (eye(Nv) - P_B) * (eye(Nv) - S_P_E) * V1;

    A2 = zeros(N, N);
    if Nu > 0
        A2(1:Nu, 1:Nu) = E' * W * V2 * E;
    end
    if Nu > 0 && Nc > 0
        A2(1:Nu, Nu+1:Nu+Nc) = E' * W * V2 * B;
    end
    if Nu > 0
        A2(1:Nu, Nu+Nc+1:N) = S * E' * W * V2;
    end
    if Nc > 0 && Nu > 0
        A2(Nu+1:Nu+Nc, 1:Nu) = B' * W * V2 * E;
    end
    if Nc > 0
        A2(Nu+1:Nu+Nc, Nu+1:Nu+Nc) = B' * W * V2 * B;
    end
    if Nc > 0
        A2(Nu+1:Nu+Nc, Nu+Nc+1:N) = B' * W * V2;
    end
    if Nu > 0
        A2(Nu+Nc+1:N, 1:Nu) = (eye(Nv) - P_B) * (eye(Nv) - P_E) * V2 * E;
    end
    if Nc > 0
        A2(Nu+Nc+1:N, Nu+1:Nu+Nc) = (eye(Nv) - P_B) * (eye(Nv) - P_E) * V2 * B;
    end
    A2(Nu+Nc+1:N, Nu+Nc+1:N) = (eye(Nv) - P_B) * (eye(Nv) - S_P_E) * V2;

end

function W0 = setupInitialCondition(X, Y, theta_grid, basis, params)
% SETUPINITIALCONDITION Project Gaussian initial data onto W.

    Nx = params.Nx;
    Ny = params.Ny;
    Nv = params.Nv;
    Nu = basis.Nu;
    Nc = basis.Nc;
    N = basis.N;

    W0 = zeros(Ny, Nx, N);

    sigma = params.sigma;
    x0 = params.x0;
    y0 = params.y0;

    for i = 1:Ny
        for j = 1:Nx
            x = X(i, j);
            y = Y(i, j);

            r2 = (x - x0)^2 + (y - y0)^2;
            f_spatial = exp(-r2 / (2*sigma^2));

            f_v = f_spatial * ones(Nv, 1) / (2*pi);

            U = basis.E' * basis.W * f_v;
            C = basis.B' * basis.W * f_v;

            f_reconstructed = basis.E * U + basis.B * C;
            R = f_v - f_reconstructed;

            W0(i, j, :) = [U; C; R];
        end
    end

end

function [W_all, times] = timeIntegration(W0, A1, A2, x_grid, y_grid, params, basis)
% TIMEINTEGRATION Forward Euler time integration.

    Nx = params.Nx;
    Ny = params.Ny;
    N = size(W0, 3);
    Nt = params.Nt;
    dt = params.dt;
    dx = x_grid(2) - x_grid(1);
    dy = y_grid(2) - y_grid(1);

    n_snapshots = params.n_snapshots;
    snapshot_interval = max(1, floor(Nt / n_snapshots));
    W_all = zeros(Ny, Nx, N, n_snapshots + 1);
    times = zeros(n_snapshots + 1, 1);

    W_all(:, :, :, 1) = W0;
    times(1) = 0;
    snapshot_idx = 2;

    W_current = W0;

    for n = 1:Nt
        dWdx = computeSpatialDerivative(W_current, basis, dx, 'x', params);
        dWdy = computeSpatialDerivative(W_current, basis, dy, 'y', params);

        for i = 1:Ny
            for j = 1:Nx
                W_local = squeeze(W_current(i, j, :));
                dWdx_local = squeeze(dWdx(i, j, :));
                dWdy_local = squeeze(dWdy(i, j, :));

                W_new = W_local - dt * (A1 * dWdx_local + A2 * dWdy_local);
                W_current(i, j, :) = W_new;
            end
        end

        if mod(n, snapshot_interval) == 0 && snapshot_idx <= n_snapshots + 1
            W_all(:, :, :, snapshot_idx) = W_current;
            times(snapshot_idx) = n * dt;
            snapshot_idx = snapshot_idx + 1;
        end
    end

    if snapshot_idx <= n_snapshots + 1
        W_all(:, :, :, end) = W_current;
        times(end) = Nt * dt;
    end

end

function dWdx = computeSpatialDerivative(W, basis, h, direction, params)
% COMPUTESPATIALDERIVATIVE Compute spatial derivative.

    [Ny, Nx, N] = size(W);
    Nu = basis.Nu;
    Nc = basis.Nc;
    Nv = basis.Nv;
    fdOrder = params.fdOrder;

    dWdx = zeros(Ny, Nx, N);

    if strcmp(direction, 'x')
        v_comp = basis.v_x;
    else
        v_comp = basis.v_y;
    end

    if strcmp(direction, 'x')
        for i = 1:Ny
            for j = 1:Nx
                if fdOrder == 1
                    if j == 1
                        dWdx(i, j, 1:Nu+Nc) = (W(i, j+1, 1:Nu+Nc) - W(i, j, 1:Nu+Nc)) / h;
                    elseif j == Nx
                        dWdx(i, j, 1:Nu+Nc) = (W(i, j, 1:Nu+Nc) - W(i, j-1, 1:Nu+Nc)) / h;
                    else
                        dWdx(i, j, 1:Nu+Nc) = (W(i, j+1, 1:Nu+Nc) - W(i, j-1, 1:Nu+Nc)) / (2*h);
                    end
                else
                    if j == 1
                        dWdx(i, j, 1:Nu+Nc) = (-3*W(i, j, 1:Nu+Nc) + 4*W(i, j+1, 1:Nu+Nc) - W(i, j+2, 1:Nu+Nc)) / (2*h);
                    elseif j == Nx
                        dWdx(i, j, 1:Nu+Nc) = (3*W(i, j, 1:Nu+Nc) - 4*W(i, j-1, 1:Nu+Nc) + W(i, j-2, 1:Nu+Nc)) / (2*h);
                    else
                        dWdx(i, j, 1:Nu+Nc) = (W(i, j+1, 1:Nu+Nc) - W(i, j-1, 1:Nu+Nc)) / (2*h);
                    end
                end

                for k = 1:Nv
                    idx = Nu + Nc + k;
                    vx = v_comp(k);

                    if fdOrder == 1
                        if j == 1
                            dWdx(i, j, idx) = (W(i, j+1, idx) - W(i, j, idx)) / h;
                        elseif j == Nx
                            dWdx(i, j, idx) = (W(i, j, idx) - W(i, j-1, idx)) / h;
                        else
                            if vx > 0
                                dWdx(i, j, idx) = (W(i, j, idx) - W(i, j-1, idx)) / h;
                            else
                                dWdx(i, j, idx) = (W(i, j+1, idx) - W(i, j, idx)) / h;
                            end
                        end
                    else
                        if j <= 1
                            dWdx(i, j, idx) = (W(i, j+1, idx) - W(i, j, idx)) / h;
                        elseif j >= Nx
                            dWdx(i, j, idx) = (W(i, j, idx) - W(i, j-1, idx)) / h;
                        else
                            if vx > 0
                                if j == 2
                                    dWdx(i, j, idx) = (W(i, j, idx) - W(i, j-1, idx)) / h;
                                else
                                    dWdx(i, j, idx) = (3*W(i, j, idx) - 4*W(i, j-1, idx) + W(i, j-2, idx)) / (2*h);
                                end
                            else
                                if j == Nx-1
                                    dWdx(i, j, idx) = (W(i, j+1, idx) - W(i, j, idx)) / h;
                                else
                                    dWdx(i, j, idx) = (-3*W(i, j, idx) + 4*W(i, j+1, idx) - W(i, j+2, idx)) / (2*h);
                                end
                            end
                        end
                    end
                end
            end
        end
    else
        for i = 1:Ny
            for j = 1:Nx
                if fdOrder == 1
                    if i == 1
                        dWdx(i, j, 1:Nu+Nc) = (W(i+1, j, 1:Nu+Nc) - W(i, j, 1:Nu+Nc)) / h;
                    elseif i == Ny
                        dWdx(i, j, 1:Nu+Nc) = (W(i, j, 1:Nu+Nc) - W(i-1, j, 1:Nu+Nc)) / h;
                    else
                        dWdx(i, j, 1:Nu+Nc) = (W(i+1, j, 1:Nu+Nc) - W(i-1, j, 1:Nu+Nc)) / (2*h);
                    end
                else
                    if i == 1
                        dWdx(i, j, 1:Nu+Nc) = (-3*W(i, j, 1:Nu+Nc) + 4*W(i+1, j, 1:Nu+Nc) - W(i+2, j, 1:Nu+Nc)) / (2*h);
                    elseif i == Ny
                        dWdx(i, j, 1:Nu+Nc) = (3*W(i, j, 1:Nu+Nc) - 4*W(i-1, j, 1:Nu+Nc) + W(i-2, j, 1:Nu+Nc)) / (2*h);
                    else
                        dWdx(i, j, 1:Nu+Nc) = (W(i+1, j, 1:Nu+Nc) - W(i-1, j, 1:Nu+Nc)) / (2*h);
                    end
                end

                for k = 1:Nv
                    idx = Nu + Nc + k;
                    vy = v_comp(k);

                    if fdOrder == 1
                        if i == 1
                            dWdx(i, j, idx) = (W(i+1, j, idx) - W(i, j, idx)) / h;
                        elseif i == Ny
                            dWdx(i, j, idx) = (W(i, j, idx) - W(i-1, j, idx)) / h;
                        else
                            if vy > 0
                                dWdx(i, j, idx) = (W(i, j, idx) - W(i-1, j, idx)) / h;
                            else
                                dWdx(i, j, idx) = (W(i+1, j, idx) - W(i, j, idx)) / h;
                            end
                        end
                    else
                        if i <= 1
                            dWdx(i, j, idx) = (W(i+1, j, idx) - W(i, j, idx)) / h;
                        elseif i >= Ny
                            dWdx(i, j, idx) = (W(i, j, idx) - W(i-1, j, idx)) / h;
                        else
                            if vy > 0
                                if i == 2
                                    dWdx(i, j, idx) = (W(i, j, idx) - W(i-1, j, idx)) / h;
                                else
                                    dWdx(i, j, idx) = (3*W(i, j, idx) - 4*W(i-1, j, idx) + W(i-2, j, idx)) / (2*h);
                                end
                            else
                                if i == Ny-1
                                    dWdx(i, j, idx) = (W(i+1, j, idx) - W(i, j, idx)) / h;
                                else
                                    dWdx(i, j, idx) = (-3*W(i, j, idx) + 4*W(i+1, j, idx) - W(i+2, j, idx)) / (2*h);
                                end
                            end
                        end
                    end
                end
            end
        end
    end

end

function [f_exact, times] = computeExactSolution(X, Y, theta_grid, t_final, params)
% COMPUTEEXACTSOLUTION Compute exact solution at final time.

    Nx = params.Nx;
    Ny = params.Ny;
    Nv = params.Nv;

    f_exact = zeros(Ny, Nx, Nv);
    times = t_final;

    sigma = params.sigma;
    x0 = params.x0;
    y0 = params.y0;

    v_x = cos(theta_grid);
    v_y = sin(theta_grid);

    for k = 1:Nv
        vx = v_x(k);
        vy = v_y(k);

        for i = 1:Ny
            for j = 1:Nx
                x = X(i, j);
                y = Y(i, j);

                x_char = x - vx * t_final;
                y_char = y - vy * t_final;

                r2_char = (x_char - x0)^2 + (y_char - y0)^2;
                f_exact(i, j, k) = exp(-r2_char / (2*sigma^2)) / (2*pi);
            end
        end
    end

end

function integral = trapz2d(f, dx, dy)
% TRAPZ2D 2D trapezoidal integration.

    [Ny, Nx] = size(f);
    weights = ones(Ny, Nx);
    weights([1, end], :) = 0.5;
    weights(:, [1, end]) = 0.5;
    weights([1, end], [1, end]) = 0.25;

    integral = dx * dy * sum(sum(weights .* f));

end
