function dev01_rbf_vs_sn()
% DEV01_RBF_VS_SN Compare discrete ordinate vs orthogonal decomposition RBF methods.
%
%   Solves the kinetic transport equation ∂_t f + v·∇_x f = 0
%   using two methods:
%   1. Standard discrete ordinate (finite difference)
%   2. Orthogonal decomposition with RBF interpolation

    clc; close all; 

    % Setup parameters and grids
    params = setupParameters();
    [x_grid, y_grid, v_grids, X, Y] = setupGrids(params);

    % Update params with actual grid sizes (may have been adjusted)
    params.Nx = length(x_grid);
    params.Ny = length(y_grid);

    % Setup initial conditions for each method
    f0_do = setupInitialCondition(X, Y, params, params.Nv_do);
    f0_rbf = setupInitialCondition(X, Y, params, params.Nv_rbf);
    f0_exact = setupInitialCondition(X, Y, params, params.Nv_exact);

    % Time array
    times = linspace(0, params.T_final, params.Nt+1);

    % Solve using discrete ordinate method
    fprintf('Solving with discrete ordinate method (%d velocities)...\n', params.Nv_do);
    f_do = solveDiscreteOrdinate(f0_do, params, x_grid, y_grid, v_grids.do, times);

    % Solve using orthogonal decomposition RBF method
    fprintf('Solving with orthogonal RBF method (%d velocities)...\n', params.Nv_rbf);
    [f_rbf, rho_rbf_analytical] = solveOrthogonalRBF(f0_rbf, params, x_grid, y_grid, v_grids.rbf, times);

    % Compute exact solution
    fprintf('Computing exact solution (%d velocities)...\n', params.Nv_exact);
    f_exact = computeExactSolution(f0_exact, params, x_grid, y_grid, v_grids.exact, times);

    % Analyze results
    analyzeResults(f_exact, f_do, f_rbf, rho_rbf_analytical, params, x_grid, y_grid, times);

end

function params = setupParameters()
% SETUPPARAMETERS Initialize problem parameters.

    % Spatial domain
    params.x_min = -0.6; params.x_max = 0.6;
    params.y_min = -0.6; params.y_max = 0.6;
    params.Nx = 64; params.Ny = 64;

    % Velocity directions (discrete ordinates on unit circle)
    params.Nv_do = 8;      % Discrete ordinate method
    params.Nv_rbf = 8;      % Orthogonal RBF method
    params.Nv_exact = 32;   % Exact solution (high resolution)

    % Time stepping
    params.T_final = 0.25;
    params.Nt = 25;
    params.dt = params.T_final / params.Nt;

    % Initial condition: sharp Gaussian
    params.sigma = 0.01;  % Very sharp, delta-like
    params.x0 = 0; params.y0 = 0;  % Centered

    % Orthogonal decomposition parameters
    params.Nu = 3;  % Number of macroscopic modes: 1, vx, vy

    % RBF parameters
    params.Nc = 32;   % Number of RBF centers (more than Nv for good representation)
    params.epsilon = 1;  % RBF shape parameter

end

function [x_grid, y_grid, v_grids, X, Y] = setupGrids(params)
% SETUPGRIDS Create computational grids.

    % Spatial grids - ensure origin is on the grid
    % For symmetric domain [-L, L], use odd number of points to include origin
    if mod(params.Nx, 2) == 0
        % Even number of points: adjust to odd
        params.Nx = params.Nx + 1;
    end
    if mod(params.Ny, 2) == 0
        % Even number of points: adjust to odd
        params.Ny = params.Ny + 1;
    end

    x_grid = linspace(params.x_min, params.x_max, params.Nx);
    y_grid = linspace(params.y_min, params.y_max, params.Ny);
    [X, Y] = meshgrid(x_grid, y_grid);

    % Verify origin is on grid
    [~, ix_center] = min(abs(x_grid));
    [~, iy_center] = min(abs(y_grid));
    fprintf('Grid center: x = %.6f, y = %.6f\n', x_grid(ix_center), y_grid(iy_center));

    % Create velocity grids for each method
    v_grids = struct();

    % Discrete ordinate method
    theta_do = linspace(0, 2*pi, params.Nv_do+1);
    theta_do = theta_do(1:end-1);
    v_grids.do = [cos(theta_do)', sin(theta_do)'];

    % Orthogonal RBF method
    theta_rbf = linspace(0, 2*pi, params.Nv_rbf+1);
    theta_rbf = theta_rbf(1:end-1);
    v_grids.rbf = [cos(theta_rbf)', sin(theta_rbf)'];

    % Exact solution
    theta_exact = linspace(0, 2*pi, params.Nv_exact+1);
    theta_exact = theta_exact(1:end-1);
    v_grids.exact = [cos(theta_exact)', sin(theta_exact)'];

end

function f0 = setupInitialCondition(X, Y, params, Nv)
% SETUPINITIALCONDITION Create sharp Gaussian initial condition.

    r2 = (X - params.x0).^2 + (Y - params.y0).^2;
    f_spatial = exp(-r2 / (2 * params.sigma^2));

    % Same initial condition for all velocity directions
    f0 = repmat(f_spatial, [1, 1, Nv]);

end

function f_all = solveDiscreteOrdinate(f0, params, x_grid, y_grid, v_grid, times)
% SOLVEDISCRETEORDINATE Solve using upwind finite differences.

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

            % Apply operator splitting: first x-direction, then y-direction
            f_temp = f_current(:,:,k);
            f_x = f_temp;  % Store original values for x-sweep

            % x-direction transport with proper upwind (all at time level n)
            if abs(vx) > 1e-12  % Only apply if velocity is significant
                if vx > 0
                    % Positive velocity: use backward difference (upwind)
                    for i = 2:size(f_x,2)
                        f_temp(:,i) = f_x(:,i) - ...
                            params.dt * vx * (f_x(:,i) - f_x(:,i-1)) / dx;
                    end
                    % Left boundary: zero inflow
                    f_temp(:,1) = f_x(:,1) - ...
                        params.dt * vx * (f_x(:,1) - 0) / dx;
                else
                    % Negative velocity: use forward difference (upwind)
                    for i = size(f_x,2)-1:-1:1
                        f_temp(:,i) = f_x(:,i) - ...
                            params.dt * vx * (f_x(:,i+1) - f_x(:,i)) / dx;
                    end
                    % Right boundary: zero inflow
                    f_temp(:,end) = f_x(:,end) - ...
                        params.dt * vx * (0 - f_x(:,end)) / dx;
                end
            end

            % y-direction transport with proper upwind (use updated values from x-sweep)
            f_y = f_temp;  % Store x-updated values for y-sweep
            if abs(vy) > 1e-12  % Only apply if velocity is significant
                if vy > 0
                    % Positive velocity: use backward difference (upwind)
                    for j = 2:size(f_y,1)
                        f_temp(j,:) = f_y(j,:) - ...
                            params.dt * vy * (f_y(j,:) - f_y(j-1,:)) / dy;
                    end
                    % Bottom boundary: zero inflow
                    f_temp(1,:) = f_y(1,:) - ...
                        params.dt * vy * (f_y(1,:) - 0) / dy;
                else
                    % Negative velocity: use forward difference (upwind)
                    for j = size(f_y,1)-1:-1:1
                        f_temp(j,:) = f_y(j,:) - ...
                            params.dt * vy * (f_y(j+1,:) - f_y(j,:)) / dy;
                    end
                    % Top boundary: zero inflow
                    f_temp(end,:) = f_y(end,:) - ...
                        params.dt * vy * (0 - f_y(end,:)) / dy;
                end
            end

            f_new(:,:,k) = f_temp;
        end

        f_all(:,:,:,n+1) = f_new;
    end

end

function [f_all, rho_all] = solveOrthogonalRBF(f0, params, x_grid, y_grid, v_grid, times)
% SOLVEORTHOGONALRBF Solve using discrete ordinate method + RBF interpolation.
% Returns both reconstructed f and analytical mass density from unified interpolation.

    Nv = size(v_grid, 1);
    Nt = length(times);
    f_all = zeros(params.Ny, params.Nx, Nv, Nt);
    rho_all = zeros(params.Ny, params.Nx, Nt);
    f_all(:,:,:,1) = f0;

    % Setup RBF centers in velocity space
    rbf_centers_v = setupRBFCentersVelocity(params);

    % Compute initial analytical mass density
    [~, rho_all(:,:,1)] = applyRBFVelocityInterpolation(f0, rbf_centers_v, params, v_grid);

    % Current solution
    f_current = f0;

    for n = 1:Nt-1
        dt = times(n+1) - times(n);

        % Step 1: Apply one time step of discrete ordinate method
        times_single = [times(n), times(n+1)];
        f_evolved = solveDiscreteOrdinate(f_current, params, x_grid, y_grid, v_grid, times_single);
        f_after_transport = f_evolved(:,:,:,2);

        % Step 2: Apply RBF interpolation in velocity space
        [f_rbf, rho_rbf] = applyRBFVelocityInterpolation(f_after_transport, rbf_centers_v, params, v_grid);

        % Store result and update for next iteration
        f_all(:,:,:,n+1) = f_rbf;
        rho_all(:,:,n+1) = rho_rbf;
        f_current = f_rbf;
    end

end


function rbf_centers_v = setupRBFCentersVelocity(params)
% SETUPRBFCENTERSVELOCITY Setup RBF centers in velocity space (unit circle).

    % Distribute RBF centers uniformly on the unit circle in velocity space
    fprintf('Setting up %d RBF centers...\n', params.Nc);
    theta_rbf = linspace(0, 2*pi, params.Nc+1);
    theta_rbf = theta_rbf(1:end-1);  % Remove duplicate point

    % RBF centers on unit circle
    rbf_centers_v = [cos(theta_rbf)', sin(theta_rbf)'];

    fprintf('Created %d RBF centers on unit circle\n', size(rbf_centers_v, 1));

end

function [s_eval, rho_analytic] = unifiedInterpolation(eval_points, centers, values, epsilon, poly_degree)
% UNIFIEDINTERPOLATION Unified kernel and polynomial interpolation from paper equation (1).
%
%   [s_eval, rho_analytic] = unifiedInterpolation(eval_points, centers, values, epsilon, poly_degree)
%   implements the unified interpolant: s(x) = Σ c_k φ(ε||x-x_k||) + Σ d_j p_j(x)
%   subject to interpolation constraints s(x_k) = f(x_k) and
%   polynomial reproduction constraints Σ c_k p_j(x_k) = 0.
%   Also returns analytical mass density from coefficients.

    N_data = length(values);      % Number of data points
    N_centers = size(centers, 1); % Number of centers
    M = size(eval_points, 1);     % Number of evaluation points

    % For unified interpolation, we interpolate at data points (velocity grid points)
    % and use centers for the kernel expansion
    data_points = eval_points(1:N_data, :); % Use first N_data points as data locations

    % Build polynomial basis matrices
    P_data = buildPolynomialBasis(data_points, poly_degree);
    P_centers = buildPolynomialBasis(centers, poly_degree);
    P_eval = buildPolynomialBasis(eval_points, poly_degree);
    Q = size(P_data, 2);     % Number of polynomial basis functions

    % Build kernel matrix from data points to centers
    Phi_data_centers = zeros(N_data, N_centers);
    for i = 1:N_data
        for j = 1:N_centers
            r = norm(data_points(i,:) - centers(j,:));
            Phi_data_centers(i,j) = wendlandC2(r, epsilon);
        end
    end

    % Build evaluation kernel matrix from eval points to centers
    Phi_eval = zeros(M, N_centers);
    for i = 1:M
        for j = 1:N_centers
            r = norm(eval_points(i,:) - centers(j,:));
            Phi_eval(i,j) = wendlandC2(r, epsilon);
        end
    end

    % Solve modified unified interpolation system
    % [Phi_data_centers P_data; P_centers^T 0] [c; d] = [values; 0]
    [c, d] = solveUnifiedSystemEfficient(Phi_data_centers, P_data, P_centers, values);

    % Evaluate interpolant at evaluation points
    s_eval = Phi_eval * c + P_eval * d;

    % Compute analytical mass density from unified interpolation coefficients
    % ρ = ∫ s(v) dv = ∫ [Σ c_k φ(||v-v_k||) + Σ d_j p_j(v)] dv
    % = Σ c_k ∫ φ(||v-v_k||) dv + Σ d_j ∫ p_j(v) dv
    rho_analytic = computeAnalyticalMassDensity(c, d, centers, epsilon);

end

function P = buildPolynomialBasis(points, degree)
% BUILDPOLYNOMIALBASIS Build polynomial basis matrix for 2D points.
%
%   P = buildPolynomialBasis(points, degree) constructs polynomial basis
%   for degree 0: {1}, degree 1: {1, x, y}, degree 2: {1, x, y, x^2, xy, y^2}

    N = size(points, 1);
    x = points(:, 1);
    y = points(:, 2);

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

function [f_filtered, rho_reconstructed] = applyRBFVelocityInterpolation(f_input, rbf_centers_v, params, v_grid)
% APPLYRBFVELOCITYINTERPOLATION Apply unified RBF interpolation in velocity space.
% Uses the unified kernel and polynomial interpolation from paper equation (1).
% Returns both reconstructed f and analytically computed mass density.

    [Ny, Nx, Nv] = size(f_input);
    f_filtered = zeros(size(f_input));
    rho_reconstructed = zeros(Ny, Nx);

    % Polynomial degree for reproduction
    poly_degree = 1; % Linear polynomials: {1, vx, vy}

    % For each spatial point, apply unified interpolation
    for i = 1:Ny
        for j = 1:Nx
            % Get velocity distribution at this spatial point
            f_v = squeeze(f_input(i,j,:));

            % Apply unified interpolation from paper equation (1)
            [f_reconstructed, rho_analytic] = unifiedInterpolation(v_grid, rbf_centers_v, f_v, ...
                                                                  params.epsilon, poly_degree);

            f_filtered(i,j,:) = f_reconstructed;
            rho_reconstructed(i,j) = rho_analytic;
        end
    end

end

function rho = computeAnalyticalMassDensity(c, d, centers, epsilon)
% COMPUTEANALYTICALMASSDESITY Compute mass density analytically from RBF coefficients.
%
%   rho = computeAnalyticalMassDensity(c, d, centers, epsilon) computes:
%   ρ = ∫ s(v) dv = Σ c_k ∫ φ(||v-v_k||) dv + Σ d_j ∫ p_j(v) dv

    N_centers = length(c);
    rho = 0;

    % Contribution from RBF kernels: Σ c_k ∫ φ(||v-v_k||) dv
    for k = 1:N_centers
        % Analytical integral of Wendland C² kernel over entire velocity space
        % ∫ φ(||v-v_k||) dv = ∫ ψ(r) r dr dθ from 0 to support radius R = 1/ε
        kernel_integral = analyticalWendlandIntegral(epsilon);
        rho = rho + c(k) * kernel_integral;
    end

    % Contribution from polynomial basis: Σ d_j ∫ p_j(v) dv
    % For linear polynomials {1, vx, vy}:
    % ∫ 1 dv = ∞ (over infinite domain) - but we integrate over bounded support
    % ∫ vx dv = 0 (odd function over symmetric domain)
    % ∫ vy dv = 0 (odd function over symmetric domain)
    % For practical purposes, we assume the support is bounded by kernel supports
    if length(d) >= 1
        % d(1) corresponds to constant polynomial: ∫ 1 dv over bounded support
        support_area = N_centers * analyticalWendlandIntegral(epsilon);
        rho = rho + d(1) * support_area;
    end
    % d(2) and d(3) correspond to vx and vy terms which integrate to 0

end

function integral = analyticalWendlandIntegral(epsilon)
% ANALYTICALWENDLANDINTEGRAL Compute ∫ φ(||v||) dv for Wendland C² kernel.
%
%   For Wendland C²: φ(r) = (1 - εr)⁴(4εr + 1) for εr < 1, 0 otherwise
%   Support radius R = 1/ε, so integral is over disk of radius R.

    % ∫₀^(1/ε) ∫₀^(2π) φ(r) r dr dθ = 2π ∫₀^(1/ε) φ(r) r dr
    % Substitution: s = εr, r = s/ε, dr = ds/ε
    % = 2π (1/ε²) ∫₀¹ (1-s)⁴(4s+1) s ds

    % Analytical result for Wendland C² kernel
    % ∫₀¹ s(1-s)⁴(4s+1) ds = ∫₀¹ [4s²(1-s)⁴ + s(1-s)⁴] ds
    % Using beta function: ∫₀¹ sᵃ(1-s)ᵇ ds = B(a+1,b+1) = Γ(a+1)Γ(b+1)/Γ(a+b+2)

    % First term: ∫₀¹ 4s²(1-s)⁴ ds = 4 * B(3,5) = 4 * (2!*4!)/(6!) = 4 * 2*24/720 = 4/15
    % Second term: ∫₀¹ s(1-s)⁴ ds = B(2,5) = (1!*4!)/(5!) = 24/120 = 1/5
    inner_integral = 4/15 + 1/5;  % = 4/15 + 3/15 = 7/15

    integral = 2 * pi * inner_integral / (epsilon^2);

end

function [c, d] = solveUnifiedSystemEfficient(Phi, P_data, P_centers, f)
% SOLVEUNIFIEDSYSTEMEFFICIENT Solve unified interpolation system efficiently.
%
%   [c, d] = solveUnifiedSystemEfficient(Phi, P_data, P_centers, f) solves:
%   [Phi P_data; P_centers^T 0] [c; d] = [f; 0]
%   using block elimination as described in paper Section 2.

    N_data = size(Phi, 1);     % Number of data points
    N_centers = size(Phi, 2);  % Number of centers
    Q = size(P_data, 2);       % Number of polynomial basis functions

    if N_centers <= 1000 && Q <= 20
        % Small system: use direct method
        A = [Phi, P_data; P_centers', zeros(Q, Q)];
        b = [f; zeros(Q, 1)];
        coeffs = A \ b;
        c = coeffs(1:N_centers);
        d = coeffs(N_centers+1:N_centers+Q);
    else
        % Large system: use block elimination
        % From block structure: Phi*c + P_data*d = f, P_centers^T*c = 0

        % Step 1: Compute P_centers^T*Phi^(-1)*P_data (if Phi is invertible)
        % For rectangular Phi, use least squares approach
        if N_data >= N_centers
            % Overdetermined: use least squares
            c = pinv(Phi) * f;
            d = zeros(Q, 1); % No polynomial part needed
        else
            % Underdetermined: use minimal norm solution
            c = Phi' * ((Phi * Phi') \ f);
            d = zeros(Q, 1);
        end
    end

end

function phi = wendlandC2(r, epsilon)
% WENDLANDC2 Wendland C2 compactly supported RBF.

    s = epsilon * r;
    if s >= 1
        phi = 0;
    else
        phi = (1 - s)^4 * (4*s + 1);
    end

end

function f_exact = computeExactSolution(f0, params, x_grid, y_grid, v_grid, times)
% COMPUTEEXACTSOLUTION Compute exact solution using method of characteristics.

    Nv_exact = size(v_grid, 1);
    Nt = length(times);
    f_exact = zeros(params.Ny, params.Nx, Nv_exact, Nt);
    f_exact(:,:,:,1) = f0;

    [X, Y] = meshgrid(x_grid, y_grid);

    for n = 2:Nt
        t = times(n);

        for k = 1:Nv_exact
            vx = v_grid(k,1);
            vy = v_grid(k,2);

            % Characteristic foot
            X_char = X - vx * t;
            Y_char = Y - vy * t;

            % Evaluate initial condition at characteristic foot
            r2 = (X_char - params.x0).^2 + (Y_char - params.y0).^2;
            f_char = exp(-r2 / (2 * params.sigma^2));

            % For inflow boundaries: if characteristic foot is outside domain,
            % the value should be zero (no inflow)
            % For outflow: if current point has characteristic foot outside,
            % but the current point is inside, we still evaluate the initial condition

            % The correct condition: zero only if we're asking for data
            % that would have come from outside the initial domain
            % But since initial condition was zero outside anyway, this is automatic

            f_exact(:,:,k,n) = f_char;
        end
    end

end

function analyzeResults(f_exact, f_do, f_rbf, rho_rbf_analytical, params, x_grid, y_grid, times)
% ANALYZERESULTS Analyze and visualize results.

    % Compute masses
    mass_exact = computeMass(f_exact, x_grid, y_grid);
    mass_do = computeMass(f_do, x_grid, y_grid);
    mass_rbf = computeMass(f_rbf, x_grid, y_grid);

    % Print mass conservation
    fprintf('\nMass Conservation:\n');
    fprintf('Method\t\t\tInitial\t\tFinal\t\tChange\n');
    fprintf('Exact\t\t\t%.6f\t%.6f\t%.2e\n', mass_exact(1), mass_exact(end), ...
            abs(mass_exact(end) - mass_exact(1))/mass_exact(1));
    fprintf('Discrete Ordinate\t%.6f\t%.6f\t%.2e\n', mass_do(1), mass_do(end), ...
            abs(mass_do(end) - mass_do(1))/mass_do(1));
    fprintf('Unified RBF\t\t%.6f\t%.6f\t%.2e\n', mass_rbf(1), mass_rbf(end), ...
            abs(mass_rbf(end) - mass_rbf(1))/mass_rbf(1));

    % Additional diagnostics for exact solution
    fprintf('\nExact Solution Analysis:\n');
    fprintf('Max travel distance: %.3f (max |v| × T = 1.0 × %.3f)\n', params.T_final, params.T_final);
    fprintf('Distance to boundary: %.3f\n', min(abs([params.x_min, params.x_max, params.y_min, params.y_max])));

    % Unified interpolation diagnostics
    fprintf('\nUnified Interpolation Analysis:\n');
    fprintf('Number of velocity points (RBF): %d\n', params.Nv_rbf);
    fprintf('Number of RBF centers: %d\n', params.Nc);
    fprintf('RBF shape parameter ε: %.3f\n', params.epsilon);
    fprintf('Polynomial degree: 1 (linear: {1, vx, vy})\n');

    % Compute errors
    l2_error_do = computeL2Error(f_do, f_exact);
    l2_error_rbf = computeL2Error(f_rbf, f_exact);

    fprintf('\nL2 Errors vs Exact:\n');
    fprintf('Discrete Ordinate:\t%.6f\n', l2_error_do(end));
    fprintf('Unified RBF:\t\t%.6f\n', l2_error_rbf(end));

    % Plot results
    plotComparison(f_exact, f_do, f_rbf, rho_rbf_analytical, times, x_grid, y_grid);
    plotMassEvolution(mass_exact, mass_do, mass_rbf, times);
    plotErrorEvolution(l2_error_do, l2_error_rbf, times);

end

function mass = computeMass(f, x_grid, y_grid)
% COMPUTEMASS Compute total mass by numerical integration.

    [~, ~, Nv, Nt] = size(f);
    dx = x_grid(2) - x_grid(1);
    dy = y_grid(2) - y_grid(1);
    dv = 2*pi / Nv;  % Use actual velocity grid size

    mass = zeros(Nt, 1);
    for n = 1:Nt
        % Integrate over velocity
        rho = sum(f(:,:,:,n), 3) * dv;

        % More accurate spatial integration using trapezoidal rule with edge corrections
        [Ny, Nx] = size(rho);

        % Apply trapezoidal weights
        weights = ones(Ny, Nx);
        weights([1,end], :) = 0.5;    % Top and bottom edges
        weights(:, [1,end]) = 0.5;    % Left and right edges
        weights([1,end], [1,end]) = 0.25;  % Corners

        mass(n) = dx * dy * sum(sum(weights .* rho));
    end

end

function l2_error = computeL2Error(f_numerical, f_exact)
% COMPUTEL2ERROR Compute L2 error vs exact solution.

    [~, ~, Nv_num, Nt] = size(f_numerical);
    [~, ~, Nv_exact, ~] = size(f_exact);
    l2_error = zeros(Nt, 1);

    if Nv_num == Nv_exact
        % Same velocity grid - direct comparison
        for n = 1:Nt
            diff = f_numerical(:,:,:,n) - f_exact(:,:,:,n);
            l2_error(n) = sqrt(mean(diff(:).^2));
        end
    else
        % Different velocity grids - compare integrated densities
        for n = 1:Nt
            % Integrate over velocity to get spatial density
            rho_num = sum(f_numerical(:,:,:,n), 3) * (2*pi / Nv_num);
            rho_exact = sum(f_exact(:,:,:,n), 3) * (2*pi / Nv_exact);

            diff = rho_num - rho_exact;
            l2_error(n) = sqrt(mean(diff(:).^2));
        end
    end

end

function plotComparison(f_exact, f_do, f_rbf, rho_rbf_analytical, times, x_grid, y_grid)
% PLOTCOMPARISON Plot snapshots at different times.

    figure(1);

    % Select time snapshots
    time_indices = round(linspace(1, length(times), 4));

    % Set fixed color limits for consistent visualization
    clim = [0, 0.8];

    for i = 1:length(time_indices)
        ti = time_indices(i);
        t = times(ti);

        % Exact solution
        subplot(3, 4, i);
        rho = sum(f_exact(:,:,:,ti), 3);
        imagesc(x_grid, y_grid, rho, clim);
        colorbar; title(sprintf('Exact t=%.2f', t));
        axis equal tight;

        % Discrete ordinate
        subplot(3, 4, i + 4);
        rho = sum(f_do(:,:,:,ti), 3);
        imagesc(x_grid, y_grid, rho, clim);
        colorbar; title(sprintf('Discrete Ordinate t=%.2f', t));
        axis equal tight;

        % Unified RBF interpolation (analytical mass density)
        subplot(3, 4, i + 8);
        rho = rho_rbf_analytical(:,:,ti);
        imagesc(x_grid, y_grid, rho, clim);
        colorbar; title(sprintf('Unified RBF ρ_{analytical} t=%.2f', t));
        axis equal tight;
    end

end

function plotMassEvolution(mass_exact, mass_do, mass_rbf, times)
% PLOTMASSEVOLUTION Plot mass evolution over time.

    figure(2);

    subplot(1, 2, 1);
    plot(times, mass_exact, 'k-', 'LineWidth', 2); hold on;
    plot(times, mass_do, 'r--', 'LineWidth', 2);
    plot(times, mass_rbf, 'b:', 'LineWidth', 2);
    xlabel('Time'); ylabel('Total Mass');
    title('Mass Evolution');
    legend('Exact', 'Discrete Ordinate', 'Orthogonal RBF', 'Location', 'best');
    grid on;

    subplot(1, 2, 2);
    plot(times, mass_exact/mass_exact(1), 'k-', 'LineWidth', 2); hold on;
    plot(times, mass_do/mass_do(1), 'r--', 'LineWidth', 2);
    plot(times, mass_rbf/mass_rbf(1), 'b:', 'LineWidth', 2);
    xlabel('Time'); ylabel('Relative Mass');
    title('Mass Conservation');
    legend('Exact', 'Discrete Ordinate', 'Orthogonal RBF', 'Location', 'best');
    grid on;

end

function plotErrorEvolution(l2_error_do, l2_error_rbf, times)
% PLOTERROREVOLUTION Plot L2 error evolution.

    figure(3);
    semilogy(times, l2_error_do, 'r--', 'LineWidth', 2); hold on;
    semilogy(times, l2_error_rbf, 'b:', 'LineWidth', 2);
    xlabel('Time'); ylabel('L2 Error (log scale)');
    title('Error vs Exact Solution');
    legend('Discrete Ordinate', 'Orthogonal RBF', 'Location', 'best');
    grid on;

end