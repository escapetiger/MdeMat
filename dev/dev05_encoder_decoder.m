function dev05_encoder_decoder()
% DEV05_ENCODER_DECODER Encoder-decoder framework for transport equations.
%
%   Three-step latent space approach:
%   1. Encoder: f(x,v) → c(x)  [32 directions → 8 features]
%   2. Transport: evolve c in latent space
%   3. Decoder: c(x) → f(x,v)  [8 features → 32 directions]

    clc; close all;

    fprintf('=== Encoder-Decoder Transport Framework ===\n\n');

    % Setup parameters
    params = setupParameters();

    % Compare methods
    fprintf('Running standard discrete ordinate method...\n');
    [rho_standard, f_standard, times] = runStandardMethod(params);

    fprintf('Running encoder-decoder method...\n');
    [rho_encoded, f_encoded, times] = runEncoderDecoderMethod(params);

    fprintf('Running exact solution...\n');
    [rho_exact, f_exact] = runExactSolution(params, times);

    % Analyze results
    analyzeResults(rho_standard, rho_encoded, rho_exact, f_standard, f_encoded, f_exact, times, params);

    fprintf('\nEncoder-decoder analysis complete.\n');

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
    params.Nv_standard = 8;    % Standard discrete ordinate
    params.Nv_high = 32;       % High-resolution physical space
    params.Nv_latent = 8;      % Latent space dimension

    % Encoder-decoder parameters
    params.rbf_type = 'wendland_c2';        % RBF type: 'wendland_c2', 'gaussian', 'multiquadric'
    params.rbf_epsilon = 0.5;               % RBF shape parameter (higher = narrower support)
    params.angular_support = pi/2;          % Base angular support width
    params.regularization_lambda = 1e-8;    % Tikhonov regularization parameter
    params.normalize_partition_unity = false; % Enforce partition of unity
    params.use_regularized_encoder = false;  % Use regularized vs standard pseudoinverse

    % Initial condition
    params.ic_type = 'isotropic_source';

end

function [rho_all, f_all, times] = runStandardMethod(params)
% RUNSTANDARDMETHOD Standard discrete ordinate method.

    % Setup grids
    [X, Y, x_grid, y_grid] = setupSpatialGrids(params);
    [v_grid, ~] = setupVelocityGrid(params.Nv_standard);

    % Initial condition
    f0 = setupInitialCondition(X, Y, v_grid, params);

    % Time integration
    times = linspace(0, params.T_final, params.Nt+1);
    f_all = zeros(params.Ny, params.Nx, params.Nv_standard, length(times));
    rho_all = zeros(params.Ny, params.Nx, length(times));

    f_current = f0;
    f_all(:,:,:,1) = f_current;
    rho_all(:,:,1) = sum(f_current, 3) * (2*pi/params.Nv_standard);

    dx = (params.x_max - params.x_min) / (params.Nx - 1);
    dy = (params.y_max - params.y_min) / (params.Ny - 1);

    for t_idx = 2:length(times)
        f_current = discreteOrdinateStep(f_current, params.dt, dx, dy, v_grid);
        f_all(:,:,:,t_idx) = f_current;
        rho_all(:,:,t_idx) = sum(f_current, 3) * (2*pi/params.Nv_standard);
    end

end

function [rho_all, f_all, times] = runEncoderDecoderMethod(params)
% RUNENCODEDECODERMETHOD Encoder-decoder transport method.

    % Setup grids
    [X, Y, x_grid, y_grid] = setupSpatialGrids(params);
    [v_grid_high, ~] = setupVelocityGrid(params.Nv_high);
    [v_grid_latent, ~] = setupVelocityGrid(params.Nv_latent);

    % Setup encoder-decoder system
    [encoder, decoder] = setupEncoderDecoder(v_grid_latent, v_grid_high, params);

    % Analyze encoder-decoder properties
    analyzeEncoderDecoder(encoder, decoder, params);

    % Initial condition at high resolution
    f0_high = setupInitialCondition(X, Y, v_grid_high, params);

    % Step 1: Encode to latent space
    c0 = encodeToLatent(f0_high, encoder);

    % Time integration in latent space
    times = linspace(0, params.T_final, params.Nt+1);
    f_all = zeros(params.Ny, params.Nx, params.Nv_high, length(times));
    rho_all = zeros(params.Ny, params.Nx, length(times));

    % Initial decode for output
    f_current = decodeFromLatent(c0, decoder);
    f_all(:,:,:,1) = f_current;
    rho_all(:,:,1) = sum(f_current, 3) * (2*pi/params.Nv_high);

    dx = (params.x_max - params.x_min) / (params.Nx - 1);
    dy = (params.y_max - params.y_min) / (params.Ny - 1);

    fprintf('    Encoder-decoder: %d→%d→%d dimensions\n', ...
            params.Nv_high, params.Nv_latent, params.Nv_high);

    c_current = c0;
    for t_idx = 2:length(times)
        % Step 2: Transport in latent space
        c_current = latentSpaceTransport(c_current, params.dt, dx, dy, v_grid_high, encoder, decoder);

        % Step 3: Decode for output
        f_current = decodeFromLatent(c_current, decoder);
        f_all(:,:,:,t_idx) = f_current;
        rho_all(:,:,t_idx) = sum(f_current, 3) * (2*pi/params.Nv_high);
    end

end

function [rho_all, f_all] = runExactSolution(params, times)
% RUNEXACTSOLUTION Exact solution using method of characteristics.

    [X, Y, x_grid, y_grid] = setupSpatialGrids(params);
    [v_grid, ~] = setupVelocityGrid(64);  % High resolution

    f_all = zeros(params.Ny, params.Nx, 64, length(times));
    rho_all = zeros(params.Ny, params.Nx, length(times));

    for t_idx = 1:length(times)
        t = times(t_idx);
        for i = 1:params.Ny
            for j = 1:params.Nx
                x = x_grid(j); y = y_grid(i);
                for k = 1:64
                    x0 = x - v_grid(k,1) * t;
                    y0 = y - v_grid(k,2) * t;
                    f_all(i,j,k,t_idx) = evaluateInitialCondition(x0, y0, v_grid(k,:), params);
                end
            end
        end
        rho_all(:,:,t_idx) = sum(f_all(:,:,:,t_idx), 3) * (2*pi/64);
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

function [encoder, decoder] = setupEncoderDecoder(v_grid_latent, v_grid_high, params)
% SETUPENCODERDECODER Create encoder and decoder matrices with better reconstruction.

    Nv_latent = size(v_grid_latent, 1);   % 8
    Nv_high = size(v_grid_high, 1);       % 32

    fprintf('  Setting up encoder-decoder with controllable parameters...\n');
    fprintf('    RBF type: %s\n', params.rbf_type);
    fprintf('    Shape parameter: %.3f\n', params.rbf_epsilon);
    fprintf('    Angular support: %.3f rad\n', params.angular_support);
    fprintf('    Regularization: %.2e\n', params.regularization_lambda);

    % Create RBF basis functions
    Phi = zeros(Nv_high, Nv_latent);     % [32×8] decoder matrix

    % Shape parameter (controls support width)
    epsilon = params.rbf_epsilon / params.angular_support;  % Scale by angular support

    for i = 1:Nv_high
        for j = 1:Nv_latent
            % Compute angular distance on unit circle
            theta_high = atan2(v_grid_high(i,2), v_grid_high(i,1));
            theta_latent = atan2(v_grid_latent(j,2), v_grid_latent(j,1));

            % Angular distance with periodic boundary
            d_theta = min(abs(theta_high - theta_latent), 2*pi - abs(theta_high - theta_latent));

            % Evaluate RBF based on type
            Phi(i,j) = evaluateRBF(d_theta, epsilon, params.rbf_type);
        end
    end

    % Ensure partition of unity if requested
    if params.normalize_partition_unity
        fprintf('    Enforcing partition of unity...\n');
        for i = 1:Nv_high
            row_sum = sum(Phi(i,:));
            if row_sum > 1e-12
                % Normalize to achieve partition of unity
                Phi(i,:) = Phi(i,:) / row_sum;
            else
                % If no RBFs cover this direction (support too narrow), use nearest neighbor
                theta_high = atan2(v_grid_high(i,2), v_grid_high(i,1));
                distances = zeros(Nv_latent, 1);
                for j = 1:Nv_latent
                    theta_latent = atan2(v_grid_latent(j,2), v_grid_latent(j,1));
                    distances(j) = min(abs(theta_high - theta_latent), 2*pi - abs(theta_high - theta_latent));
                end
                [~, nearest] = min(distances);
                Phi(i, nearest) = 1.0;
            end
        end

        % Verify partition of unity
        row_sums = sum(Phi, 2);
        unity_error = max(abs(row_sums - 1));
        fprintf('    Partition of unity error: %.2e\n', unity_error);
    else
        fprintf('    Using raw RBF values (no normalization)\n');
    end

    % Store decoder
    decoder.Phi = Phi;

    % Construct encoder based on parameters
    if params.use_regularized_encoder
        fprintf('    Using regularized encoder...\n');
        % Use Tikhonov regularization for better conditioning
        lambda = params.regularization_lambda;
        Phi_reg = Phi' * Phi + lambda * eye(Nv_latent);
        encoder.Phi_pinv = Phi_reg \ Phi';

        fprintf('    Regularization parameter: %.2e\n', lambda);
        fprintf('    Original condition number: %.2e\n', cond(Phi));
        fprintf('    Regularized condition number: %.2e\n', cond(Phi_reg));
    else
        fprintf('    Using standard pseudoinverse encoder...\n');
        encoder.Phi_pinv = pinv(Phi);
        fprintf('    Condition number: %.2e\n', cond(Phi));
    end

    % Store both for comparison if needed
    encoder.Phi_pinv_std = pinv(Phi);
    encoder.params = params;  % Store parameters for reference

    % Test reconstruction quality
    testReconstructionQuality(encoder, decoder, params);

end

function testReconstructionQuality(encoder, decoder, params)
% TESTRECONSTRUCTIONQUALITY Test encoder-decoder reconstruction quality.

    fprintf('  Testing reconstruction quality...\n');

    Phi = decoder.Phi;
    Phi_pinv = encoder.Phi_pinv;

    % Test 1: Identity reconstruction
    I_approx = Phi * Phi_pinv;
    reconstruction_error = norm(I_approx - eye(size(I_approx)), 'fro');
    fprintf('    ||ΦΦ⁺ - I||_F: %.2e\n', reconstruction_error);

    % Test 2: Reconstruction of smooth test functions
    Nv_high = size(Phi, 1);
    theta_high = linspace(0, 2*pi, Nv_high+1);
    theta_high = theta_high(1:end-1);

    % Test function 1: Smooth isotropic
    f_test1 = ones(Nv_high, 1);
    c_test1 = Phi_pinv * f_test1;
    f_recon1 = Phi * c_test1;
    error1 = norm(f_test1 - f_recon1) / norm(f_test1);
    fprintf('    Isotropic reconstruction error: %.2e\n', error1);

    % Test function 2: Smooth directional
    f_test2 = 1 + 0.5 * cos(2*theta_high');
    c_test2 = Phi_pinv * f_test2;
    f_recon2 = Phi * c_test2;
    error2 = norm(f_test2 - f_recon2) / norm(f_test2);
    fprintf('    Smooth directional reconstruction error: %.2e\n', error2);

    % Test function 3: Single peak
    [~, peak_idx] = max(cos(theta_high'));
    f_test3 = zeros(Nv_high, 1);
    f_test3(peak_idx) = 1;
    c_test3 = Phi_pinv * f_test3;
    f_recon3 = Phi * c_test3;
    error3 = norm(f_test3 - f_recon3) / norm(f_test3);
    fprintf('    Single peak reconstruction error: %.2e\n', error3);

    % Check condition number
    cond_num = cond(Phi);
    fprintf('    Decoder condition number: %.2e\n', cond_num);

    if reconstruction_error > 0.1
        fprintf('    WARNING: Large reconstruction error - consider adjusting parameters\n');
    end

end

function c = encodeToLatent(f, encoder)
% ENCODETOLATENT Encode physical space to latent space: f → c.

    [Ny, Nx, Nv] = size(f);
    c = zeros(Ny, Nx, size(encoder.Phi_pinv, 1));

    for i = 1:Ny
        for j = 1:Nx
            f_local = squeeze(f(i,j,:));
            c(i,j,:) = encoder.Phi_pinv * f_local;
        end
    end

end

function f = decodeFromLatent(c, decoder)
% DECODEFROMLATENT Decode latent space to physical space: c → f.

    [Ny, Nx, Nc] = size(c);
    f = zeros(Ny, Nx, size(decoder.Phi, 1));

    for i = 1:Ny
        for j = 1:Nx
            c_local = squeeze(c(i,j,:));
            f(i,j,:) = decoder.Phi * c_local;
        end
    end

end

function c_new = latentSpaceTransport(c, dt, dx, dy, v_grid, encoder, decoder)
% LATENTSPACETRANSPORT Transport in latent space using transformed equation.

    [Ny, Nx, Nc] = size(c);

    % Compute latent transport operator: T_latent = (Φ⁺Φ)⁻¹ Φ⁺(v·Φ)
    Phi = decoder.Phi;                    % [32×8]
    Phi_pinv = encoder.Phi_pinv;          % [8×32]

    % For each latent dimension, compute transport
    dcx = zeros(size(c));
    dcy = zeros(size(c));
    dcx(:,1,:) = (c(:,2,:)-c(:,1,:))/dx;
    dcx(:,2:Nx-1,:) = (c(:,3:Nx,:)- c(:,1:Nx-2,:))/(2*dx);
    dcx(:,Nx,:) = (c(:,Nx,:)-c(:,Nx-1,:))/dx;
    dcy(1,:,:) = (c(2,:,:)-c(1,:,:))/dy;
    dcy(2:Ny-1,:,:) = (c(3:Ny,:,:)- c(1:Ny-2,:,:))/(2*dy);
    dcy(Ny,:,:) = (c(Ny,:,:)-c(Ny-1,:,:))/dy;
    vx = v_grid(:, 1);
    vy = v_grid(:, 2);
    dcx = reshape(dcx, [], 1, Nc);
    dcy = reshape(dcy, [], 1, Nc);
    Px = reshape(vx .* Phi, 1, [], Nc);
    Py = reshape(vy .* Phi, 1, [], Nc);
    dc = sum(Px .* dcx, 3) + sum(Py .* dcy, 3);
    dc = reshape(dc, [], 1, size(dc, 2));
    Q = reshape(Phi_pinv, 1, Nc, []);
    dc = reshape(sum(Q .* dc, 3), Ny, Nx, []);
    c_new = c - dt * dc;
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
        case 'isotropic_source'
            r = sqrt(x^2 + y^2);
            f_val = exp(-r^2 / (2*0.1^2)) / (2*pi);
        otherwise
            f_val = 0;
    end

end

function f_new = discreteOrdinateStep(f, dt, dx, dy, v_grid)
% DISCRETEORDINATESTEP Standard discrete ordinate time step.

    [Ny, Nx, Nv] = size(f);
    f_new = zeros(size(f));

    for k = 1:Nv
        vx = v_grid(k,1); vy = v_grid(k,2);
        for i = 1:Ny
            for j = 1:Nx
                % Upwind differences
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

function analyzeEncoderDecoder(encoder, decoder, params)
% ANALYZEENCODERDECODER Analyze encoder-decoder properties.

    fprintf('\n  Final encoder-decoder analysis:\n');

    Phi = decoder.Phi;
    Phi_pinv = encoder.Phi_pinv;

    fprintf('    Decoder matrix Φ: [%d×%d]\n', size(Phi, 1), size(Phi, 2));
    fprintf('    Encoder matrix Φ⁺: [%d×%d]\n', size(Phi_pinv, 1), size(Phi_pinv, 2));

    % Check reconstruction quality
    I_approx = Phi * Phi_pinv;
    reconstruction_error = norm(I_approx - eye(size(I_approx)), 'fro');
    fprintf('    Reconstruction error ||ΦΦ⁺ - I||_F: %.2e\n', reconstruction_error);

    % Check individual reconstruction errors
    max_diag_error = max(abs(diag(I_approx) - 1));

    % Off-diagonal error: create matrix with zeros on diagonal
    I_offdiag = I_approx - diag(diag(I_approx));
    max_offdiag_error = max(abs(I_offdiag(:)));

    fprintf('    Max diagonal error: %.2e\n', max_diag_error);
    fprintf('    Max off-diagonal error: %.2e\n', max_offdiag_error);

    % Check rank and conditioning
    fprintf('    Rank of Φ: %d (should be ≤ %d)\n', rank(Phi), min(size(Phi)));
    fprintf('    Condition number: %.2e\n', cond(Phi));

    % Check compression properties
    compression_ratio = size(Phi, 1) / size(Phi, 2);
    fprintf('    Compression ratio: %.1f:1\n', compression_ratio);

    % Estimate information loss
    [U, S, V] = svd(Phi);
    singular_values = diag(S);
    info_retained = sum(singular_values.^2) / sum(singular_values(1)^2 * ones(size(singular_values)));
    fprintf('    Information retention: %.1f%%\n', info_retained * 100);

    if reconstruction_error < 0.01
        fprintf('    ✓ Excellent reconstruction quality\n');
    elseif reconstruction_error < 0.1
        fprintf('    ✓ Good reconstruction quality\n');
    else
        fprintf('    ⚠ Poor reconstruction quality - consider adjusting parameters\n');
    end

end

function phi = evaluateRBF(r, epsilon, rbf_type)
% EVALUATERBF Evaluate radial basis function of specified type.
%
%   phi = evaluateRBF(r, epsilon, rbf_type) evaluates the RBF at distance r
%   with shape parameter epsilon.
%
%   Supported RBF types:
%   - 'wendland_c2': Wendland C2 compactly supported RBF
%   - 'wendland_c4': Wendland C4 compactly supported RBF
%   - 'gaussian': Gaussian RBF (global support)
%   - 'multiquadric': Multiquadric RBF (global support)
%   - 'inverse_multiquadric': Inverse multiquadric RBF
%   - 'thin_plate_spline': Thin plate spline RBF

    switch lower(rbf_type)
        case 'wendland_c2'
            phi = wendlandC2RBF(r, epsilon);

        case 'wendland_c4'
            phi = wendlandC4RBF(r, epsilon);

        case 'gaussian'
            phi = exp(-(epsilon * r)^2);

        case 'multiquadric'
            phi = sqrt(1 + (epsilon * r)^2);

        case 'inverse_multiquadric'
            phi = 1 / sqrt(1 + (epsilon * r)^2);

        case 'thin_plate_spline'
            if r == 0
                phi = 0;
            else
                phi = (epsilon * r)^2 * log(epsilon * r);
            end

        otherwise
            error('Unknown RBF type: %s', rbf_type);
    end

end

function phi = wendlandC2RBF(r, epsilon)
% WENDLANDC2RBF Wendland C2 compactly supported radial basis function.
%
%   Wendland C2 RBF: φ(r) = (1-εr)₊⁴(4εr + 1)
%   where (·)₊ = max(·, 0) and ε controls the support width.

    s = epsilon * r;

    if s >= 1
        phi = 0;  % Compact support
    else
        phi = (1 - s)^4 * (4*s + 1);  % Wendland C2 formula
    end

end

function phi = wendlandC4RBF(r, epsilon)
% WENDLANDC4RBF Wendland C4 compactly supported radial basis function.
%
%   Wendland C4 RBF: φ(r) = (1-εr)₊⁶(35(εr)² + 18εr + 3)
%   Higher order smoothness than C2.

    s = epsilon * r;

    if s >= 1
        phi = 0;  % Compact support
    else
        phi = (1 - s)^6 * (35*s^2 + 18*s + 3);  % Wendland C4 formula
    end

end

function analyzeResults(rho_standard, rho_encoded, rho_exact, f_standard, f_encoded, f_exact, times, params)
% ANALYZERESULTS Compare and visualize results.

    fprintf('\n  Results analysis:\n');

    % Plot density comparison
    plotDensityComparison(rho_standard, rho_encoded, rho_exact, times, params);

    % Plot angular comparison
    plotAngularComparison(f_standard, f_encoded, f_exact, times, params);

    % Compute errors
    computeErrors(rho_standard, rho_encoded, rho_exact, times);

end

function plotDensityComparison(rho_standard, rho_encoded, rho_exact, times, params)
% PLOTDENSITYCOMPARISON Plot density evolution comparison.

    [~, ~, x_grid, y_grid] = setupSpatialGrids(params);
    time_indices = [1, round(length(times)/2), length(times)];

    figure(1);
    for i = 1:length(time_indices)
        ti = time_indices(i);

        subplot(3, 3, (i-1)*3 + 1);
        imagesc(x_grid, y_grid, rho_standard(:,:,ti));
        colorbar; clim([0, 0.3]);
        title(sprintf('Standard DO, t=%.3f', times(ti)));
        axis equal tight;

        subplot(3, 3, (i-1)*3 + 2);
        imagesc(x_grid, y_grid, rho_encoded(:,:,ti));
        colorbar; clim([0, 0.3]);
        title(sprintf('Encoder-Decoder, t=%.3f', times(ti)));
        axis equal tight;

        subplot(3, 3, (i-1)*3 + 3);
        imagesc(x_grid, y_grid, rho_exact(:,:,ti));
        colorbar; clim([0, 0.3]);
        title(sprintf('Exact, t=%.3f', times(ti)));
        axis equal tight;
    end

    sgtitle('Density Evolution Comparison');

end

function plotAngularComparison(f_standard, f_encoded, f_exact, times, params)
% PLOTANGULARCOMPARISON Compare angular distributions.

    center_i = round(params.Ny/2);
    center_j = round(params.Nx/2);

    theta_standard = linspace(0, 2*pi, params.Nv_standard+1);
    theta_standard = theta_standard(1:end-1);

    theta_encoded = linspace(0, 2*pi, params.Nv_high+1);
    theta_encoded = theta_encoded(1:end-1);

    theta_exact = linspace(0, 2*pi, size(f_exact,3)+1);
    theta_exact = theta_exact(1:end-1);

    figure(2);
    time_indices = [1, round(length(times)/2), length(times)];

    for i = 1:length(time_indices)
        ti = time_indices(i);

        subplot(1, 3, i);

        f_std = squeeze(f_standard(center_i, center_j, :, ti));
        plot(theta_standard, f_std, 'ro-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Standard DO');
        hold on;

        f_enc = squeeze(f_encoded(center_i, center_j, :, ti));
        plot(theta_encoded, f_enc, 'bs-', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'Encoder-Decoder');

        f_ex = squeeze(f_exact(center_i, center_j, :, ti));
        plot(theta_exact, f_ex, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Exact');

        xlabel('Angle (rad)'); ylabel('f');
        title(sprintf('Angular Distribution, t=%.3f', times(ti)));
        legend('Location', 'best'); grid on;
    end

    sgtitle('Angular Distribution Comparison');

end

function computeErrors(rho_standard, rho_encoded, rho_exact, times)
% COMPUTEERRORS Compute and display error metrics.

    l2_error_std = zeros(size(times));
    l2_error_enc = zeros(size(times));

    for ti = 1:length(times)
        diff_std = rho_standard(:,:,ti) - rho_exact(:,:,ti);
        diff_enc = rho_encoded(:,:,ti) - rho_exact(:,:,ti);

        l2_error_std(ti) = sqrt(mean(diff_std(:).^2));
        l2_error_enc(ti) = sqrt(mean(diff_enc(:).^2));
    end

    figure(3);
    semilogy(times, l2_error_std, 'r-', 'LineWidth', 2, 'DisplayName', 'Standard DO');
    hold on;
    semilogy(times, l2_error_enc, 'b-', 'LineWidth', 2, 'DisplayName', 'Encoder-Decoder');
    xlabel('Time'); ylabel('L2 Error (log scale)');
    title('Error Evolution');
    legend('Location', 'best'); grid on;

    fprintf('    Final L2 errors:\n');
    fprintf('      Standard DO: %.6e\n', l2_error_std(end));
    fprintf('      Encoder-Decoder: %.6e\n', l2_error_enc(end));
    fprintf('      Improvement ratio: %.2f\n', l2_error_std(end) / l2_error_enc(end));

end