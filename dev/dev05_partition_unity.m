function dev05_partition_unity()
% DEV05_PARTITION_UNITY Implement partition of unity method for velocity coupling.
%
%   Compares standard discrete ordinate method with partition of unity
%   approach that naturally introduces angular coupling through overlapping
%   basis functions without adding artificial terms to transport equation.

    clc; close all;

    fprintf('=== Partition of Unity Method for Velocity Coupling ===\n\n');

    % Setup parameters
    params = setupParameters();

    % Compare methods
    fprintf('Running standard discrete ordinate method...\n');
    [rho_standard, f_standard, times] = runStandardMethod(params);

    fprintf('Running partition of unity method...\n');
    [rho_partition, f_partition, times, partition_basis] = runPartitionUnityMethod(params);

    fprintf('Running exact solution...\n');
    [rho_exact, f_exact] = runExactSolution(params, times);

    % Analyze and visualize results
    analyzeResults(rho_standard, rho_partition, rho_exact, f_standard, f_partition, f_exact, times, params, partition_basis);

    fprintf('\nPartition of unity analysis complete.\n');

end

function params = setupParameters()
% SETUPPARAMETERS Initialize simulation parameters.

    % Spatial domain
    params.x_min = -1; params.x_max = 1;
    params.y_min = -1; params.y_max = 1;
    params.Nx = 51; params.Ny = 51;

    % Time domain
    params.T_final = 0.8;
    params.Nt = 80;
    params.dt = params.T_final / params.Nt;

    % Velocity space
    params.Nv = 8; % Discrete ordinate directions (standard resolution)
    params.Nv_high = 32; % High resolution for initial representation

    % Partition of unity parameters
    params.angular_support = pi/2; % Angular width of basis functions
    params.basis_type = 'hat'; % 'hat', 'gaussian', 'bspline'
    params.overlap_factor = 2; % Controls basis function overlap

    % Initial condition parameters
    params.ic_type = 'isotropic_source';
    params.beam_direction = [1, 0.5]; % Initial beam direction
    params.beam_width = 0.2;
    params.beam_center = [-0.5, 0];

end

function [rho_all, f_all, times] = runStandardMethod(params)
% RUNSTANDARDMETHOD Run standard discrete ordinate method.

    % Setup grids
    [X, Y, x_grid, y_grid] = setupSpatialGrids(params);
    [v_grid, theta] = setupVelocityGrid(params.Nv);

    % Initial condition
    f0 = setupInitialCondition(X, Y, v_grid, params);

    % Time integration
    times = linspace(0, params.T_final, params.Nt+1);
    f_all = zeros(params.Ny, params.Nx, params.Nv, length(times));
    rho_all = zeros(params.Ny, params.Nx, length(times));

    f_current = f0;
    f_all(:,:,:,1) = f_current;
    rho_all(:,:,1) = sum(f_current, 3) * (2*pi/params.Nv);

    dx = (params.x_max - params.x_min) / (params.Nx - 1);
    dy = (params.y_max - params.y_min) / (params.Ny - 1);

    for t_idx = 2:length(times)
        f_current = discreteOrdinateStep(f_current, params.dt, dx, dy, v_grid);
        f_all(:,:,:,t_idx) = f_current;
        rho_all(:,:,t_idx) = sum(f_current, 3) * (2*pi/params.Nv);
    end

end

function [rho_all, f_all, times, partition_basis] = runPartitionUnityMethod(params)
% RUNPARTITIONUNITYMETHOD Run high-resolution partition of unity method.

    % Setup grids
    [X, Y, x_grid, y_grid] = setupSpatialGrids(params);
    [v_grid_8, theta_8] = setupVelocityGrid(params.Nv);           % 8 standard directions
    [v_grid_32, theta_32] = setupVelocityGrid(params.Nv_high);    % 32 high-res directions

    % Setup partition of unity basis systems
    partition_basis = setupHighResPartitionBasis(v_grid_8, v_grid_32, params);

    % Analyze basis properties
    analyzeHighResBasisProperties(partition_basis, params);

    % Diagnose potential rank deficiency issues
    diagnoseRankDeficiency(partition_basis, params);

    % Initial condition at HIGH resolution (32 directions)
    f0_high = setupInitialCondition(X, Y, v_grid_32, params);

    % Project high-res initial condition to 8 coefficients
    c0 = projectHighResToCoefficients(f0_high, partition_basis);

    % Time integration
    times = linspace(0, params.T_final, params.Nt+1);
    f_all = zeros(params.Ny, params.Nx, params.Nv_high, length(times));  % 32 directions output
    rho_all = zeros(params.Ny, params.Nx, length(times));

    % Initial reconstruction at high resolution (32 directions)
    c_current = c0;
    f_current = reconstructHighResFromCoefficients(c_current, partition_basis);
    f_all(:,:,:,1) = f_current;
    rho_all(:,:,1) = sum(f_current, 3) * (2*pi/params.Nv_high);

    dx = (params.x_max - params.x_min) / (params.Nx - 1);
    dy = (params.y_max - params.y_min) / (params.Ny - 1);

    fprintf('    Weak form partition unity: solve %dx%d least squares per spatial point\n', ...
            params.Nv_high, params.Nv);

    for t_idx = 2:length(times)
        % Solve weak form: Φ(c^{n+1}-c^n)/dt + v·Φ∇c^n = 0 via least squares
        c_current = partitionUnityStep(c_current, params.dt, dx, dy, v_grid_32, partition_basis);

        % Reconstruct at 32 high-resolution directions
        f_current = reconstructHighResFromCoefficients(c_current, partition_basis);
        f_all(:,:,:,t_idx) = f_current;
        rho_all(:,:,t_idx) = sum(f_current, 3) * (2*pi/params.Nv_high);
    end

end

function [rho_all, f_all] = runExactSolution(params, times)
% RUNEXACTSOLUTION Compute exact solution using method of characteristics.

    [X, Y, x_grid, y_grid] = setupSpatialGrids(params);
    [v_grid, theta] = setupVelocityGrid(64); % High resolution for exact

    f_all = zeros(params.Ny, params.Nx, length(v_grid), length(times));
    rho_all = zeros(params.Ny, params.Nx, length(times));

    for t_idx = 1:length(times)
        t = times(t_idx);

        for i = 1:params.Ny
            for j = 1:params.Nx
                x = x_grid(j);
                y = y_grid(i);

                for k = 1:length(v_grid)
                    % Characteristic: trace back from (x,y) at time t
                    x0 = x - v_grid(k,1) * t;
                    y0 = y - v_grid(k,2) * t;

                    % Evaluate initial condition at (x0, y0)
                    f_all(i,j,k,t_idx) = evaluateInitialCondition(x0, y0, v_grid(k,:), params);
                end
            end
        end

        rho_all(:,:,t_idx) = sum(f_all(:,:,:,t_idx), 3) * (2*pi/length(v_grid));
    end

end

function [X, Y, x_grid, y_grid] = setupSpatialGrids(params)
% SETUPSPATIALGRIDS Create spatial grids.

    x_grid = linspace(params.x_min, params.x_max, params.Nx);
    y_grid = linspace(params.y_min, params.y_max, params.Ny);
    [X, Y] = meshgrid(x_grid, y_grid);

end

function [v_grid, theta] = setupVelocityGrid(Nv)
% SETUPVELOCITYGRID Create velocity grid.

    theta = linspace(0, 2*pi, Nv+1);
    theta = theta(1:end-1);
    v_grid = [cos(theta)', sin(theta)'];

end

function partition_basis = setupHighResPartitionBasis(v_grid_8, v_grid_32, params)
% SETUPHIGHRESPARTITIONBASIS Create high-resolution partition of unity basis.

    partition_basis = struct();
    partition_basis.v_centers_8 = v_grid_8;      % 8 coefficient centers
    partition_basis.v_high_32 = v_grid_32;       % 32 high-res directions
    partition_basis.support_width = params.angular_support;
    partition_basis.basis_type = params.basis_type;

    % Basis matrix 1: [32×8] - project 32 directions to 8 coefficients
    partition_basis.phi_32to8 = zeros(32, 8);
    for i = 1:32
        phi_values = createPartitionBasisFunction(v_grid_32(i,:), v_grid_8, params);
        partition_basis.phi_32to8(i,:) = phi_values;
    end

    % Basis matrix 2: [32×8] - reconstruct at 32 high-resolution directions
    partition_basis.phi_32to8_recon = zeros(32, 8);
    for i = 1:32
        phi_values = createPartitionBasisFunction(v_grid_32(i,:), v_grid_8, params);
        partition_basis.phi_32to8_recon(i,:) = phi_values;
    end

    % Basis matrix 3: [8×8] - reconstruct at 8 standard directions (for compatibility)
    partition_basis.phi_8to8 = zeros(8, 8);
    for i = 1:8
        phi_values = createPartitionBasisFunction(v_grid_8(i,:), v_grid_8, params);
        partition_basis.phi_8to8(i,:) = phi_values;
    end

    % Store legacy matrix for compatibility
    partition_basis.phi_matrix = partition_basis.phi_8to8;

end

function partition_basis = setupPartitionBasis(v_grid, params)
% SETUPPARTITIONBASIS Create partition of unity basis functions (legacy).

    Nv = size(v_grid, 1);
    partition_basis = struct();
    partition_basis.v_centers = v_grid;
    partition_basis.support_width = params.angular_support;
    partition_basis.basis_type = params.basis_type;

    % Precompute basis matrix: phi_matrix(i,j) = φ_j(v_i)
    partition_basis.phi_matrix = zeros(Nv, Nv);

    for i = 1:Nv
        phi_values = createPartitionBasisFunction(v_grid(i,:), v_grid, params);
        partition_basis.phi_matrix(i,:) = phi_values;
    end

end

function phi = createPartitionBasisFunction(v_eval, v_centers, params)
% CREATEPARTITIONBASISFUNCTION Create partition of unity basis function values.

    Nv = size(v_centers, 1);
    phi = zeros(1, Nv);

    % Get angle of evaluation point
    theta_eval = atan2(v_eval(2), v_eval(1));

    for k = 1:Nv
        % Get angle of center k
        theta_k = atan2(v_centers(k,2), v_centers(k,1));

        % Angular distance (handle periodicity)
        d_theta = min(abs(theta_eval - theta_k), 2*pi - abs(theta_eval - theta_k));

        % Basis function value
        switch params.basis_type
            case 'hat'
                if d_theta <= params.angular_support
                    phi(k) = max(0, 1 - d_theta/params.angular_support);
                end
            case 'gaussian'
                sigma = params.angular_support / 3; % 3-sigma support
                phi(k) = exp(-d_theta^2 / (2*sigma^2));
            case 'bspline'
                % Quadratic B-spline
                s = d_theta / params.angular_support;
                if s <= 0.5
                    phi(k) = 0.75 - s^2;
                elseif s <= 1.5
                    phi(k) = 0.5 * (1.5 - s)^2;
                end
        end
    end

    % Normalize to ensure partition of unity
    if sum(phi) > 1e-12
        phi = phi / sum(phi);
    end

end

function f0 = setupInitialCondition(X, Y, v_grid, params)
% SETUPINITIALCONDITION Create initial condition.

    [Ny, Nx] = size(X);
    Nv = size(v_grid, 1);
    f0 = zeros(Ny, Nx, Nv);

    for k = 1:Nv
        for i = 1:Ny
            for j = 1:Nx
                f0(i,j,k) = evaluateInitialCondition(X(i,j), Y(i,j), v_grid(k,:), params);
            end
        end
    end

end

function f_val = evaluateInitialCondition(x, y, v, params)
% EVALUATEINITIALCONDITION Evaluate initial condition at given point.

    switch params.ic_type
        case 'gaussian_beam'
            % Spatial Gaussian
            x0 = params.beam_center(1);
            y0 = params.beam_center(2);
            spatial_factor = exp(-((x-x0)^2 + (y-y0)^2) / (2*params.beam_width^2));

            % Directional delta function (approximate with narrow Gaussian)
            vx_target = params.beam_direction(1);
            vy_target = params.beam_direction(2);
            norm_target = sqrt(vx_target^2 + vy_target^2);
            vx_target = vx_target / norm_target;
            vy_target = vy_target / norm_target;

            % Angular difference
            dot_product = v(1)*vx_target + v(2)*vy_target;
            angular_factor = exp(50*(dot_product - 1)); % Sharp directional beam

            f_val = spatial_factor * angular_factor;

        case 'isotropic_source'
            % Isotropic source at center
            r = sqrt((x-0)^2 + (y-0)^2);
            f_val = exp(-r^2 / (2*0.1^2)) / (2*pi); % Isotropic

        otherwise
            f_val = 0;
    end

end

function c = projectHighResToCoefficients(f_high, partition_basis)
% PROJECTHIGHRESTOCOEFFICENTS Project high-res (32) distribution to coefficients (8).

    [Ny, Nx, Nv_high] = size(f_high);
    c = zeros(Ny, Nx, 8);  % 8 coefficients

    % Precompute pseudoinverse: [32×8] → [8×32]
    phi_32to8_pinv = pinv(partition_basis.phi_32to8);

    % Project at each spatial point
    for i = 1:Ny
        for j = 1:Nx
            f_local = squeeze(f_high(i,j,:));  % [32 x 1]

            % Project: min ||Φ₃₂₋₈ * c_local - f_local||²
            c_local = phi_32to8_pinv * f_local;  % [8×32] * [32×1] = [8×1]

            c(i,j,:) = c_local;
        end
    end

end

function c = projectToCoefficients(f, partition_basis)
% PROJECTTOCOEFFICIENTS Project distribution to coefficient representation using least squares.

    [Ny, Nx, Nv] = size(f);
    c = zeros(Ny, Nx, Nv);

    % Precompute pseudoinverse for efficiency
    phi_pinv = pinv(partition_basis.phi_matrix);

    % Solve least squares at each spatial point
    for i = 1:Ny
        for j = 1:Nx
            f_local = squeeze(f(i,j,:));  % [Nv x 1]

            % Solve: min ||Φ * c_local - f_local||²
            c_local = phi_pinv * f_local;

            c(i,j,:) = c_local;
        end
    end

end

function f = reconstructFromCoefficients(c, partition_basis)
% RECONSTRUCTFROMCOEFFICIENTS Reconstruct distribution from coefficients at 8 standard directions.

    [Ny, Nx, Nc] = size(c);  % Nc = 8 coefficients
    f = zeros(Ny, Nx, 8);    % Reconstruct at 8 standard directions

    % Use 8×8 reconstruction matrix
    phi_8to8 = partition_basis.phi_8to8;

    % Reconstruct at each spatial point
    for i = 1:Ny
        for j = 1:Nx
            c_local = squeeze(c(i,j,:));  % [8×1]
            f_local = phi_8to8 * c_local; % [8×8] * [8×1] = [8×1]
            f(i,j,:) = f_local;
        end
    end

end

function f = reconstructHighResFromCoefficients(c, partition_basis)
% RECONSTRUCTHIGHRESFROMCOEFFICIENTS Reconstruct at 32 high-resolution directions from 8 coefficients.

    [Ny, Nx, Nc] = size(c);  % Nc = 8 coefficients
    f = zeros(Ny, Nx, 32);   % Reconstruct at 32 high-resolution directions

    % Use 32×8 reconstruction matrix
    phi_32to8_recon = partition_basis.phi_32to8_recon;

    % Reconstruct at each spatial point
    for i = 1:Ny
        for j = 1:Nx
            c_local = squeeze(c(i,j,:));  % [8×1]
            f_local = phi_32to8_recon * c_local; % [32×8] * [8×1] = [32×1]
            f(i,j,:) = f_local;
        end
    end

end

function f_new = discreteOrdinateStep(f, dt, dx, dy, v_grid)
% DISCRETEORDINATESTEP Standard discrete ordinate time step.

    [Ny, Nx, Nv] = size(f);
    f_new = zeros(size(f));

    for k = 1:Nv
        vx = v_grid(k,1);
        vy = v_grid(k,2);

        for i = 1:Ny
            for j = 1:Nx
                % Upwind finite differences
                if vx > 0
                    dfx = (f(i,j,k) - f(i,max(1,j-1),k)) / dx;
                else
                    dfx = (f(i,min(Nx,j+1),k) - f(i,j,k)) / dx;
                end

                if vy > 0
                    dfy = (f(i,j,k) - f(max(1,i-1),j,k)) / dy;
                else
                    dfy = (f(min(Ny,i+1),j,k) - f(i,j,k)) / dy;
                end

                f_new(i,j,k) = f(i,j,k) - dt * (vx*dfx + vy*dfy);
            end
        end
    end

end

function c_new = partitionUnityStep(c, dt, dx, dy, v_grid_32, partition_basis)
% PARTITIONUNITYSTEP Weak form transport: solve Φ(c^{n+1}-c^n)/dt + v·Φ∇c^n = 0.

    [Ny, Nx, Nc] = size(c);  % Nc = 8 coefficients
    c_new = zeros(size(c));

    % Extract reconstruction matrix Φ [32×8]
    Phi = partition_basis.phi_32to8_recon;

    % Solve weak form at each spatial point
    for i = 1:Ny
        for j = 1:Nx
            c_new(i,j,:) = solveWeakFormAtPoint(c, i, j, dt, dx, dy, v_grid_32, Phi, Ny, Nx, i, j);
        end
    end

end

function c_new_local = solveWeakFormAtPoint(c, i, j, dt, dx, dy, v_grid_32, Phi, Ny, Nx, i_idx, j_idx)
% SOLVEWEAKFORMATPOINT Solve weak form least squares at a spatial point.

    Nc = size(c, 3);  % 8 coefficients
    Nv = size(v_grid_32, 1);  % 32 velocity directions

    % Current coefficients at point (i,j)
    c_local = squeeze(c(i, j, :));  % [8×1]

    % Build least squares system: A * c^{n+1} = b
    A = zeros(Nv, Nc);  % [32×8]
    b = zeros(Nv, 1);   % [32×1]

    for k = 1:Nv  % For each velocity direction
        vx = v_grid_32(k, 1);
        vy = v_grid_32(k, 2);

        % Spatial derivatives of coefficients using upwind differences
        dcx = zeros(Nc, 1);
        dcy = zeros(Nc, 1);

        for m = 1:Nc  % For each coefficient
            % x-direction spatial derivative
            if vx > 0
                c_left = c(i, max(1, j-1), m);
                dcx(m) = (c_local(m) - c_left) / dx;
            else
                c_right = c(i, min(Nx, j+1), m);
                dcx(m) = (c_right - c_local(m)) / dx;
            end

            % y-direction spatial derivative
            if vy > 0
                c_down = c(max(1, i-1), j, m);
                dcy(m) = (c_local(m) - c_down) / dy;
            else
                c_up = c(min(Ny, i+1), j, m);
                dcy(m) = (c_up - c_local(m)) / dy;
            end
        end

        % Weak form equation for direction k:
        % Φ_k·(c^{n+1} - c^n)/dt + v_k·Φ_k·∇c^n = 0
        % Rearranged: Φ_k·c^{n+1}/dt = Φ_k·c^n/dt - v_k·Φ_k·∇c^n

        phi_k = Phi(k, :)';  % [8×1] - basis functions for direction k

        % Left-hand side: Φ_k * c^{n+1} / dt
        A(k, :) = phi_k' / dt;

        % Right-hand side: Φ_k * c^n / dt - v_k * Φ_k * ∇c^n
        rhs_time = sum(phi_k .* c_local) / dt;
        rhs_spatial = vx * sum(phi_k .* dcx) + vy * sum(phi_k .* dcy);
        b(k) = rhs_time - rhs_spatial;
    end

    % Check condition number and add regularization if needed
    condition_num = cond(A);

    if condition_num > 1e12 || rank(A) < Nc
        % Add Tikhonov regularization: min ||A*c - b||² + λ||c - c_old||²
        lambda = 1e-6 * max(diag(A'*A));
        A_reg = [A; sqrt(lambda) * eye(Nc)];
        b_reg = [b; sqrt(lambda) * c_local];
        c_new_local = A_reg \ b_reg;

        if condition_num > 1e12
            fprintf('  Warning: Ill-conditioned system at (%d,%d), cond=%.2e, using regularization\n', i_idx, j_idx, condition_num);
        end
    else
        % Standard least squares
        c_new_local = A \ b;
    end

end

function analyzeHighResBasisProperties(partition_basis, params)
% ANALYZEHIGHRESBASISPROPERTIES Analyze high-resolution partition basis properties.

    fprintf('  High-resolution basis analysis:\n');
    fprintf('    Basis type: %s\n', params.basis_type);
    fprintf('    Angular support: %.3f rad (%.1f degrees)\n', params.angular_support, params.angular_support*180/pi);

    % Analyze 32→8 projection matrix
    phi_32to8 = partition_basis.phi_32to8;
    fprintf('    Projection matrix [32×8]:\n');
    fprintf('      Condition number: %.2e\n', cond(phi_32to8));

    % Check if each high-res direction contributes to coefficients
    for k = 1:8
        contrib_count = sum(abs(phi_32to8(:,k)) > 1e-6);
        fprintf('      Coefficient %d receives from %d/32 directions\n', k, contrib_count);
    end

    % Analyze 32×8 high-resolution reconstruction matrix
    phi_32to8_recon = partition_basis.phi_32to8_recon;
    fprintf('    High-res reconstruction matrix [32×8]:\n');
    fprintf('      Condition number: %.2e\n', cond(phi_32to8_recon));

    % Check if each high-res output gets contributions from multiple coefficients
    for k = 1:32
        contrib_count = sum(abs(phi_32to8_recon(k,:)) > 1e-6);
        if mod(k-1, 8) == 0  % Print every 8th direction
            fprintf('      Direction %d influenced by %d/8 coefficients\n', k, contrib_count);
        end
    end

    % Analyze 8×8 compatibility reconstruction matrix
    phi_8to8 = partition_basis.phi_8to8;
    fprintf('    Compatibility reconstruction matrix [8×8]:\n');
    fprintf('      Condition number: %.2e\n', cond(phi_8to8));

    % Check partition of unity for 8×8 reconstruction
    row_sums = sum(phi_8to8, 2);
    unity_error = max(abs(row_sums - 1));
    fprintf('      Partition of unity error: %.2e\n', unity_error);

    % Check coupling strength in high-res reconstruction
    non_zero_entries = sum(abs(phi_32to8_recon) > 1e-6, 'all');
    total_entries = 32 * 8;
    sparsity = 1 - non_zero_entries / total_entries;
    fprintf('      High-res reconstruction sparsity: %.1f%%\n', sparsity * 100);

    % Check maximum coupling strength
    coupling_strength = max(abs(phi_32to8_recon(:)));
    fprintf('      Max reconstruction weight: %.3f\n', coupling_strength);

    if coupling_strength > 0.1
        fprintf('      → Strong angular coupling detected!\n');
    elseif coupling_strength > 0.01
        fprintf('      → Moderate angular coupling detected.\n');
    else
        fprintf('      → Weak angular coupling.\n');
    end

end

function diagnoseRankDeficiency(partition_basis, params)
% DIAGNOSErankdeficiency Diagnose potential rank deficiency in weak form system.

    fprintf('  Rank deficiency diagnosis:\n');

    Phi = partition_basis.phi_32to8_recon;  % [32×8]
    dt = params.dt;

    % Analyze the matrix A = Phi/dt that will be used in weak form
    A_template = Phi / dt;

    % Check rank and condition number
    rank_A = rank(A_template);
    cond_A = cond(A_template);

    fprintf('    Matrix A = Φ/dt properties:\n');
    fprintf('      Size: [%d×%d]\n', size(A_template, 1), size(A_template, 2));
    fprintf('      Rank: %d (should be %d)\n', rank_A, size(A_template, 2));
    fprintf('      Condition number: %.2e\n', cond_A);

    if rank_A < size(A_template, 2)
        fprintf('      → Rank deficient! %d missing dimensions\n', size(A_template, 2) - rank_A);

        % Find which directions contribute least
        [U, S, V] = svd(A_template);
        singular_values = diag(S);

        fprintf('      Singular values (largest to smallest):\n');
        for i = 1:min(8, length(singular_values))
            fprintf('        σ_%d = %.2e\n', i, singular_values(i));
        end

        % Check null space
        null_space = null(A_template');
        if ~isempty(null_space)
            fprintf('      Null space dimensions: %d\n', size(null_space, 2));
            fprintf('      → These coefficient combinations are unobservable\n');
        end
    else
        fprintf('      → Full rank system\n');
    end

    % Check basis function properties
    fprintf('    Basis function analysis:\n');
    row_sums = sum(Phi, 2);
    fprintf('      Partition of unity error: %.2e\n', max(abs(row_sums - 1)));

    % Check linear independence of basis functions
    rank_phi = rank(Phi);
    fprintf('      Rank of Φ: %d (should be ≤ %d)\n', rank_phi, min(size(Phi)));

    if rank_phi < size(Phi, 2)
        fprintf('      → Basis functions are linearly dependent!\n');
        fprintf('      → Consider reducing number of coefficients or changing support width\n');
    end

end

function analyzeBasisProperties(partition_basis, params)
% ANALYZEBASISPROPERTIES Analyze partition of unity basis properties.

    phi_matrix = partition_basis.phi_matrix;
    fprintf('  Basis analysis:\n');
    fprintf('    Basis type: %s\n', params.basis_type);
    fprintf('    Angular support: %.3f rad (%.1f degrees)\n', params.angular_support, params.angular_support*180/pi);
    fprintf('    Condition number: %.2e\n', cond(phi_matrix));

    % Check partition of unity property
    row_sums = sum(phi_matrix, 2);
    unity_error = max(abs(row_sums - 1));
    fprintf('    Partition of unity error: %.2e\n', unity_error);

    % Check symmetry
    symmetry_error = max(max(abs(phi_matrix - phi_matrix')));
    fprintf('    Matrix symmetry error: %.2e\n', symmetry_error);

    % Display basis matrix
    fprintf('    Basis matrix Φ:\n');
    for i = 1:size(phi_matrix,1)
        fprintf('      [');
        for j = 1:size(phi_matrix,2)
            fprintf(' %6.3f', phi_matrix(i,j));
        end
        fprintf(' ]\n');
    end

end

function analyzeResults(rho_standard, rho_partition, rho_exact, f_standard, f_partition, f_exact, times, params, partition_basis)
% ANALYZERESULTS Compare and visualize results.

    % Plot density evolution comparison
    plotDensityComparison(rho_standard, rho_partition, rho_exact, times, params);

    % Plot angular distribution comparison
    plotAngularComparison(f_standard, f_partition, f_exact, times, params);

    % Quantitative analysis
    analyzeAngularCoupling(f_standard, f_partition, times, params);

    % Error analysis
    analyzeErrors(rho_standard, rho_partition, rho_exact, times);

    % Mass conservation analysis
    analyzeMassConservation(rho_standard, rho_partition, rho_exact, times, params);

    % Note: Coefficient correlation analysis disabled for high-res output
    % (would need to store coefficient evolution separately)
    fprintf('    Note: Coefficient correlation analysis requires coefficient time series\n');

end

function plotDensityComparison(rho_standard, rho_partition, rho_exact, times, params)
% PLOTDENSITYCOMPARISON Plot density evolution for all methods.

    [~, ~, x_grid, y_grid] = setupSpatialGrids(params);

    % Select time points for comparison
    time_indices = [1, round(length(times)/3), round(2*length(times)/3), length(times)];

    figure(1);
    for i = 1:length(time_indices)
        ti = time_indices(i);

        % Standard DO method
        subplot(4, 3, (i-1)*3 + 1);
        imagesc(x_grid, y_grid, rho_standard(:,:,ti));
        colorbar; clim([0, 0.3]);
        title(sprintf('Standard DO, t=%.3f', times(ti)));
        xlabel('x'); ylabel('y'); axis equal tight;

        % Partition of unity method
        subplot(4, 3, (i-1)*3 + 2);
        imagesc(x_grid, y_grid, rho_partition(:,:,ti));
        colorbar; clim([0, 0.3]);
        title(sprintf('Partition Unity, t=%.3f', times(ti)));
        xlabel('x'); ylabel('y'); axis equal tight;

        % Exact solution
        subplot(4, 3, (i-1)*3 + 3);
        imagesc(x_grid, y_grid, rho_exact(:,:,ti));
        colorbar; clim([0, 0.3]);
        title(sprintf('Exact, t=%.3f', times(ti)));
        xlabel('x'); ylabel('y'); axis equal tight;
    end

    sgtitle('Density Evolution Comparison');

end

function plotAngularComparison(f_standard, f_partition, f_exact, times, params)
% PLOTANGULARCOMPARISON Compare angular distributions.

    % Extract angular distribution at a fixed spatial point
    center_i = round(params.Ny/2);
    center_j = round(params.Nx/2);

    % Angular grids for different methods
    theta_standard = linspace(0, 2*pi, params.Nv+1);         % 8 directions
    theta_standard = theta_standard(1:end-1);

    theta_partition = linspace(0, 2*pi, params.Nv_high+1);   % 32 directions
    theta_partition = theta_partition(1:end-1);

    theta_exact = linspace(0, 2*pi, size(f_exact,3)+1);      % Exact resolution
    theta_exact = theta_exact(1:end-1);

    figure(2);

    time_indices = [1, round(length(times)/2), length(times)];

    for i = 1:length(time_indices)
        ti = time_indices(i);

        subplot(1, 3, i);

        % Standard method (8 directions)
        f_std = squeeze(f_standard(center_i, center_j, :, ti));
        plot(theta_standard, f_std, 'ro-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Standard DO');
        hold on;

        % Partition method (32 directions)
        f_part = squeeze(f_partition(center_i, center_j, :, ti));
        plot(theta_partition, f_part, 'bs-', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Partition Unity');

        % Exact solution
        f_ex = squeeze(f_exact(center_i, center_j, :, ti));
        plot(theta_exact, f_ex, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Exact');

        xlabel('Angle (rad)'); ylabel('f');
        title(sprintf('Angular Distribution, t=%.3f', times(ti)));
        legend('Location', 'best');
        grid on;
    end

    sgtitle('Angular Distribution Comparison at Center');

end

function analyzeAngularCoupling(f_standard, f_partition, times, params)
% ANALYZEANGULARCOUPLING Analyze the effectiveness of angular coupling.

    center_i = round(params.Ny/2);
    center_j = round(params.Nx/2);

    angular_variance_std = zeros(size(times));
    angular_variance_part = zeros(size(times));

    % Angular grids for different methods
    theta_std = linspace(0, 2*pi, params.Nv+1);
    theta_std = theta_std(1:end-1);

    theta_part = linspace(0, 2*pi, params.Nv_high+1);
    theta_part = theta_part(1:end-1);

    for ti = 1:length(times)
        % Standard method angular variance
        f_ang_std = squeeze(f_standard(center_i, center_j, :, ti));
        if sum(f_ang_std) > 1e-10
            f_norm = f_ang_std / sum(f_ang_std);
            z = sum(f_norm(:) .* exp(1i * theta_std(:)));
            angular_variance_std(ti) = 1 - abs(z);
        end

        % Partition method angular variance
        f_ang_part = squeeze(f_partition(center_i, center_j, :, ti));
        if sum(f_ang_part) > 1e-10
            f_norm = f_ang_part / sum(f_ang_part);
            z = sum(f_norm(:) .* exp(1i * theta_part(:)));
            angular_variance_part(ti) = 1 - abs(z);
        end
    end

    figure(3);
    plot(times, angular_variance_std, 'r-', 'LineWidth', 2, 'DisplayName', 'Standard DO');
    hold on;
    plot(times, angular_variance_part, 'b-', 'LineWidth', 2, 'DisplayName', 'Partition Unity');
    xlabel('Time'); ylabel('Angular Variance');
    title('Angular Diffusion Analysis');
    legend('Location', 'best');
    grid on;

    fprintf('  Angular coupling analysis:\n');
    fprintf('    Standard DO final variance: %.6f\n', angular_variance_std(end));
    fprintf('    Partition Unity final variance: %.6f\n', angular_variance_part(end));
    fprintf('    Improvement ratio: %.2f\n', angular_variance_part(end) / angular_variance_std(end));

end

function analyzeErrors(rho_standard, rho_partition, rho_exact, times)
% ANALYZEERRORS Compute and display error metrics.

    l2_error_std = zeros(size(times));
    l2_error_part = zeros(size(times));

    for ti = 1:length(times)
        diff_std = rho_standard(:,:,ti) - rho_exact(:,:,ti);
        diff_part = rho_partition(:,:,ti) - rho_exact(:,:,ti);

        l2_error_std(ti) = sqrt(mean(diff_std(:).^2));
        l2_error_part(ti) = sqrt(mean(diff_part(:).^2));
    end

    figure(4);
    semilogy(times, l2_error_std, 'r-', 'LineWidth', 2, 'DisplayName', 'Standard DO');
    hold on;
    semilogy(times, l2_error_part, 'b-', 'LineWidth', 2, 'DisplayName', 'Partition Unity');
    xlabel('Time'); ylabel('L2 Error (log scale)');
    title('Error Evolution');
    legend('Location', 'best');
    grid on;

    fprintf('  Error analysis:\n');
    fprintf('    Standard DO final L2 error: %.6e\n', l2_error_std(end));
    fprintf('    Partition Unity final L2 error: %.6e\n', l2_error_part(end));
    fprintf('    Error reduction ratio: %.2f\n', l2_error_std(end) / l2_error_part(end));

end

function analyzeMassConservation(rho_standard, rho_partition, rho_exact, times, params)
% ANALYZEMASSCONSERVATION Analyze mass conservation properties.

    dx = (params.x_max - params.x_min) / (params.Nx - 1);
    dy = (params.y_max - params.y_min) / (params.Ny - 1);

    mass_std = zeros(size(times));
    mass_part = zeros(size(times));
    mass_exact = zeros(size(times));

    for ti = 1:length(times)
        mass_std(ti) = sum(sum(rho_standard(:,:,ti))) * dx * dy;
        mass_part(ti) = sum(sum(rho_partition(:,:,ti))) * dx * dy;
        mass_exact(ti) = sum(sum(rho_exact(:,:,ti))) * dx * dy;
    end

    figure(5);
    plot(times, mass_std, 'r-', 'LineWidth', 2, 'DisplayName', 'Standard DO');
    hold on;
    plot(times, mass_part, 'b-', 'LineWidth', 2, 'DisplayName', 'Partition Unity');
    plot(times, mass_exact, 'k--', 'LineWidth', 2, 'DisplayName', 'Exact');
    xlabel('Time'); ylabel('Total Mass');
    title('Mass Conservation');
    legend('Location', 'best');
    grid on;

    fprintf('  Mass conservation analysis:\n');
    fprintf('    Initial masses: Std=%.6f, Part=%.6f, Exact=%.6f\n', ...
            mass_std(1), mass_part(1), mass_exact(1));
    fprintf('    Final masses: Std=%.6f, Part=%.6f, Exact=%.6f\n', ...
            mass_std(end), mass_part(end), mass_exact(end));
    fprintf('    Mass loss: Std=%.2e%%, Part=%.2e%%\n', ...
            100*(mass_std(1)-mass_std(end))/mass_std(1), ...
            100*(mass_part(1)-mass_part(end))/mass_part(1));

end

function analyzeCoefficientCorrelation(f_standard, f_partition, partition_basis, times, params)
% ANALYZECOEFFICIENTCORRELATION Analyze correlation between coefficients.

    fprintf('\n  Coefficient Correlation Analysis:\n');

    % Convert f values to coefficients for analysis
    [Ny, Nx, Nv_partition] = size(f_partition);  % f_partition now has 32 directions

    % Sample spatial points for analysis
    sample_points = [
        round(Ny*0.3), round(Nx*0.3);   % Point 1
        round(Ny*0.5), round(Nx*0.5);   % Center
        round(Ny*0.7), round(Nx*0.7)    % Point 3
    ];

    figure(6);

    % Analyze correlations at different times
    time_indices = [1, round(length(times)/2), length(times)];

    for t_idx = 1:length(time_indices)
        ti = time_indices(t_idx);

        subplot(2, 3, t_idx);

        % Extract coefficients at sample points
        coeffs_at_points = [];
        labels = {};

        for sp = 1:size(sample_points, 1)
            i = sample_points(sp, 1);
            j = sample_points(sp, 2);

            % Get function values and project to coefficients
            f_local = squeeze(f_partition(i, j, :, ti));
            c_local = projectToCoefficients(reshape(f_local, 1, 1, length(f_local)), partition_basis);
            c_local = squeeze(c_local);

            coeffs_at_points = [coeffs_at_points, c_local];
            labels{sp} = sprintf('Point %d', sp);
        end

        % Plot coefficients
        Nc = size(coeffs_at_points, 1);
        for k = 1:Nc
            plot(1:size(coeffs_at_points, 2), coeffs_at_points(k, :), 'o-', ...
                 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', sprintf('c_%d', k));
            hold on;
        end

        xlabel('Spatial Points'); ylabel('Coefficient Values');
        title(sprintf('Coefficients at t=%.3f', times(ti)));
        legend('Location', 'best'); grid on;

        % Compute correlation matrix
        subplot(2, 3, t_idx + 3);

        % Correlation between coefficients across all spatial points
        Nc = 8;  % Number of coefficients
        all_coeffs = zeros(Nc, Ny*Nx);
        idx = 1;
        for i = 1:Ny
            for j = 1:Nx
                f_local = squeeze(f_partition(i, j, :, ti));
                if sum(abs(f_local)) > 1e-12  % Only non-zero points
                    c_local = projectToCoefficients(reshape(f_local, 1, 1, length(f_local)), partition_basis);
                    all_coeffs(:, idx) = squeeze(c_local);
                    idx = idx + 1;
                end
            end
        end

        % Trim unused columns
        all_coeffs = all_coeffs(:, 1:idx-1);

        if size(all_coeffs, 2) > 1
            % Compute correlation matrix
            corr_matrix = corrcoef(all_coeffs');

            imagesc(corr_matrix);
            colorbar; clim([-1, 1]);
            title(sprintf('Coefficient Correlation t=%.3f', times(ti)));
            xlabel('Coefficient Index'); ylabel('Coefficient Index');

            % Add correlation values as text
            for i = 1:Nc
                for j = 1:Nc
                    text(j, i, sprintf('%.2f', corr_matrix(i,j)), ...
                         'HorizontalAlignment', 'center', 'Color', 'white', 'FontSize', 8);
                end
            end
        else
            text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center');
        end
    end

    sgtitle('Coefficient Analysis: Values and Correlations');

    % Compute temporal correlation analysis
    analyzeTemporalCorrelation(f_partition, partition_basis, times, params);

    % Compare with standard discrete ordinate
    compareCorrelationPatterns(f_standard, f_partition, partition_basis, times, params);

end

function analyzeTemporalCorrelation(f_partition, partition_basis, times, params)
% ANALYZETEMPORALCORRELATION Analyze how coefficient correlations evolve in time.

    fprintf('    Temporal correlation analysis:\n');

    % Sample a fixed spatial point
    center_i = round(size(f_partition, 1)/2);
    center_j = round(size(f_partition, 2)/2);

    % Extract coefficient time series
    Nv_output = size(f_partition, 3);  % This should be 8
    Nc = 8;  % Number of coefficients
    coeff_time_series = zeros(Nc, length(times));

    for ti = 1:length(times)
        f_local = squeeze(f_partition(center_i, center_j, :, ti));
        if sum(abs(f_local)) > 1e-12
            c_local = projectToCoefficients(reshape(f_local, 1, 1, length(f_local)), partition_basis);
            coeff_time_series(:, ti) = squeeze(c_local);
        end
    end

    % Compute correlations between adjacent coefficients
    adjacent_correlations = zeros(Nc, 1);
    for k = 1:Nc
        k_next = mod(k, Nc) + 1;  % Periodic boundary

        if std(coeff_time_series(k, :)) > 1e-12 && std(coeff_time_series(k_next, :)) > 1e-12
            correlation = corrcoef(coeff_time_series(k, :), coeff_time_series(k_next, :));
            adjacent_correlations(k) = correlation(1, 2);
        end
    end

    fprintf('      Adjacent coefficient correlations:\n');
    for k = 1:Nc
        k_next = mod(k, Nc) + 1;
        fprintf('        corr(c_%d, c_%d) = %.3f\n', k, k_next, adjacent_correlations(k));
    end

    mean_correlation = mean(adjacent_correlations);
    fprintf('      Mean adjacent correlation: %.3f\n', mean_correlation);

    if mean_correlation > 0.5
        fprintf('      → Strong coupling detected!\n');
    elseif mean_correlation > 0.1
        fprintf('      → Moderate coupling detected.\n');
    else
        fprintf('      → Weak/no coupling detected.\n');
    end

end

function compareCorrelationPatterns(f_standard, f_partition, partition_basis, times, params)
% COMPARECORRELATIONPATTERNS Compare correlation patterns between methods.

    fprintf('    Comparison with standard discrete ordinate:\n');

    % Sample spatial point
    center_i = round(size(f_standard, 1)/2);
    center_j = round(size(f_standard, 2)/2);

    % Final time analysis
    ti = length(times);

    % Standard method: f values directly
    f_std = squeeze(f_standard(center_i, center_j, :, ti));

    % Partition method: convert to coefficients
    f_part = squeeze(f_partition(center_i, center_j, :, ti));
    c_part = projectToCoefficients(reshape(f_part, 1, 1, length(f_part)), partition_basis);
    c_part = squeeze(c_part);

    % Compute correlations between adjacent directions
    corr_std = zeros(length(f_std)-1, 1);
    corr_part = zeros(length(c_part)-1, 1);

    for k = 1:length(f_std)-1
        if std([f_std(k), f_std(k+1)]) > 1e-12
            corr_std(k) = f_std(k) * f_std(k+1) / (norm(f_std(k:k+1))^2);
        end

        if std([c_part(k), c_part(k+1)]) > 1e-12
            corr_part(k) = c_part(k) * c_part(k+1) / (norm(c_part(k:k+1))^2);
        end
    end

    fprintf('      Standard DO adjacent correlations: %.3f ± %.3f\n', ...
            mean(corr_std), std(corr_std));
    fprintf('      Partition Unity adjacent correlations: %.3f ± %.3f\n', ...
            mean(corr_part), std(corr_part));

    % Angular smoothness measure
    smoothness_std = -sum(abs(diff(f_std)));  % Negative total variation
    smoothness_part = -sum(abs(diff(c_part)));

    fprintf('      Angular smoothness: Standard=%.3f, Partition=%.3f\n', ...
            smoothness_std, smoothness_part);

    if smoothness_part > smoothness_std
        fprintf('      → Partition method produces smoother distributions!\n');
    else
        fprintf('      → No significant smoothness improvement.\n');
    end

end