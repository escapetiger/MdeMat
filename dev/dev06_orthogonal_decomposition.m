function dev06_orthogonal_decomposition()
% DEV06_ORTHOGONAL_DECOMPOSITION Compare three solvers for kinetic transport.
%
%   Solves the kinetic transport equation:
%       ∂_t f + v·∇_x f = 0
%   on 2D spatial domain with velocity v on unit circle S^1.
%
%   Three solvers compared:
%   1. Method of Characteristics (MoC): Exact solution via backward tracking
%   2. Discrete Ordinate Method (DoM): Direct discretization with upwind FD
%   3. URM Method: U+R decomposition
%
% See also:
%   configureParameters, configureGrids, encode, decode

    clc; close all;

    fprintf('[M] Simulation started.\n\n');

    fprintf('\n[M] configure parameters.\n');
    params = configureParameters();

    fprintf('\n[M] configure grids.\n');
    grid = configureGrids(params);

    fprintf('\n[M] Running MOC Solver...\n');
    sol_moc = mocRun(grid, params);

    fprintf('\n[M] Running DOM Solver...\n');
    sol_dom = domRun(grid, params);

    fprintf('\n[M] Running URM Solver...\n');
    sol_urm = urmRun(grid, params);

    fprintf('\n[M] Comparison and Analysis\n');
    compareResults(sol_moc, sol_dom, sol_urm, grid, params);

    fprintf('[M] Simulation completed.\n');
end

function params = configureParameters()
% configurePARAMETERS Initialize simulation parameters.
%
%   params = configureParameters() creates a structure containing all
%   simulation parameters including spatial/temporal discretization,
%   decomposition settings, initial/boundary conditions, and output
%   control.

    % Spatial domain
    params.Nd = 2;
    params.bbox = repmat([-1, 1], 1, params.Nd);
    params.Nxp = [51, 51];
    params.Nx = prod(params.Nxp);

    % Velocity space discretization (unit circle S^1)
    params.Nv_moc = 32;  % Quadrature points for MoC and visualization
    params.Nv_dom = 32;   % Quadrature points for DoM and visualization
    params.Nv_u = 32;  % Quadrature points for macro projection

    % URM decomposition parameters
    params.Nu = 3;  % Macro modes: {1, cos(θ), sin(θ), ...}
    params.Nr = 0;  % Residual modes (discrete velocity points)

    % RBF interpolation for residual
    params.rbfType = 'wendland_c2';  % RBF type: 'wendland_c2', 'gaussian', 'multiquadric'
    params.rbfEpsilon = 2.0;  % RBF shape parameter
    % NOTE: For Gaussian RBF, use smaller epsilon (0.5-1.5) to avoid ill-conditioning
    %       For sharp initial conditions, prefer 'wendland_c2' for better stability

    % Temporal discretization
    params.tDisc = 'ssprk2';  % Time integrator: 'fe', 'heun', 'ssprk2'
    params.tFinal = 0.5;      % Final time
    params.cfl = 0.1;         % CFL number for time step control

    % Spatial discretization
    params.xDisc = 'fd';     % Spatial method: 'fd' (finite difference)
    params.xOrder = 2;       % Finite difference order
    params.xFlux = 'alt_bf'; % Flux scheme for macro components:
    % 'ctr': central difference with artificial dissipation (Laiu et al.)
    % 'alt_fb': alternating (k<l forward, k>=l backward)
    % 'alt_bf': alternating (k<l backward, k>=l forward)
    params.thetaMinmod = 1.5; % Minmod parameter (1 < theta < 2) for artificial dissipation

    % Initial condition
    params.ic = 'beam';       % Type: 'gaussian', 'ring', 'two_peaks', 'beam'
    params.sigma = sqrt(1e-2);    % Gaussian width parameter
    params.x0 = [0, 0];           % Initial center position
    params.x1 = [-0.3, 0.0];      % First peak center for 'two_peaks'
    params.x2 = [0.3, 0.0];       % Second peak center for 'two_peaks'
    params.r0 = 0.3;              % Ring radius for 'ring' IC
    params.v0 = [0.5, 0.5];       % Beam direction for 'beam' IC
    params.sigma_v = 0.1;         % Angular width for 'beam' IC (in radians)

    % Boundary condition
    params.bc = 'zero';  % Boundary type: 'zero', 'periodic', 'reflect'

    % Filtering (for URM method)
    params.filtering = true;         % Enable spectral filtering
    params.filterType = 'rot';       % Filter type: 'rot', 'exp'
    params.filterStrength = 1;       % Filter strength
    params.filterOrder = 36;         % Filter order
    params.filterCutoff = 2/3;       % Cutoff ratio (fraction of modes)

    % Limiter (for oscillation control)
    params.limiting = false;           % Enable slope/positivity limiter
    params.limiterType = 'cutoff';    % Limiter type: 'cutoff', 'ls-r', 'opt-r'
    params.limiterEps = 1e-14;        % Minimum allowed density

    % Output control
    params.Ns = 5;  % Number of snapshots to save
end

function grid = configureGrids(params)
% configureGRIDS Create computational grids for spatial and velocity domains.
%
%   grid = configureGrids(params) creates all computational grids including
%   spatial grids, temporal grid, and velocity quadrature grids for
%   different basis projections.

    Nd = params.Nd;
    bbox = params.bbox;
    Nxp = params.Nxp;
    Nv_moc = params.Nv_moc;
    Nv_dom = params.Nv_dom;
    Nv_u = params.Nv_u;
    Nr = params.Nr;
    tFinal = params.tFinal;
    cfl = params.cfl;
    Ns = params.Ns;

    % Spatial grids
    grid.x = cell(1, Nd);
    grid.dx = zeros(1, Nd);
    for i = 1:Nd
        grid.x{i} = linspace(bbox(2*i-1), bbox(2*i), Nxp(i));
        grid.dx(i) = grid.x{i}(2) - grid.x{i}(1);
    end
    grid.h = max(grid.dx);

    % Temporal grid with CFL condition
    % CFL condition: dt <= cfl * h / |v_max|
    % For unit circle, |v_max| = 1
    grid.dt = cfl * grid.h;
    grid.Nt = ceil(tFinal / grid.dt);
    grid.dt = tFinal / grid.Nt;  % Adjust dt to hit tFinal exactly
    grid.t = linspace(0, tFinal, grid.Nt + 1);

    % Snapshot indices: select Ns + 1 uniformly distributed times from grid.t
    si = round(linspace(1, grid.Nt + 1, Ns + 1));
    si(end) = grid.Nt + 1;  % Ensure final time is included
    grid.si = si;
    grid.s = grid.t(si);

    % Velocity quadrature
    grid.angle_moc = linspace(0, 2*pi, Nv_moc + 1).';
    grid.angle_moc = grid.angle_moc(1:end-1);
    grid.v_moc = [cos(grid.angle_moc), sin(grid.angle_moc)];
    grid.w_moc = (2*pi / Nv_moc) * ones(Nv_moc, 1);

    % Velocity quadrature
    grid.angle_dom = linspace(0, 2*pi, Nv_dom + 1).';
    grid.angle_dom = grid.angle_dom(1:end-1);
    grid.v_dom = [cos(grid.angle_dom), sin(grid.angle_dom)];
    grid.w_dom = (2*pi / Nv_dom) * ones(Nv_dom, 1);

    % Velocity quadrature for macro projection
    grid.angle_u = linspace(0, 2*pi, Nv_u + 1).';
    grid.angle_u = grid.angle_u(1:end-1);
    grid.v_u = [cos(grid.angle_u), sin(grid.angle_u)];
    grid.w_u = (2*pi / Nv_u) * ones(Nv_u, 1);

    % Velocity quadrature for residual projection
    grid.angle_r = linspace(0, 2*pi, Nr + 1).';
    grid.angle_r = grid.angle_r(1:end-1);
    grid.v_r = [cos(grid.angle_r), sin(grid.angle_r)];
    grid.w_r = (2*pi / Nr) * ones(Nr, 1);

    % Report grid information
    fprintf('[M] Grid configure:\n');
    fprintf('[M] Spatial: [%s] points, dx = [%s]\n', ...
        num2str(Nxp), num2str(grid.dx));
    fprintf('[M] CFL number: %.4f\n', cfl);
    fprintf('[M] Time step: dt = %.4f (Nt = %d)\n', grid.dt, grid.Nt);
end

% ==========================================================================
% METHOD OF CHARACTERISTICS (MOC)
% ==========================================================================

function sol = mocRun(grid, params)
% RUNMOCSOLVER Run method of characteristics solver.
%
%   sol = mocRun(grid, params) solves the kinetic transport equation
%   using the method of characteristics. Returns solution structure with
%   phase space distribution @a sol.F at snapshot times @a sol.t.

    Nx = params.Nx;
    Nv = params.Nv_moc;
    Ns = params.Ns;

    fprintf('[M] Computing exact solution at %d time snapshots...\n', Ns + 1);

    % Compute solution at each snapshot time
    sol.F = zeros([Nx, Nv, Ns + 1]);

    [X, V] = computePhaseGrid(grid.x, grid.v_moc, params);
    for js = 1:(Ns + 1)
        t = grid.s(js);
        X0 = X - V * t;
        F0 = computeInitialDistribution(X0, V, params);
        F0 = reshape(F0, [Nx, Nv]);
        sol.F(:, :, js) = F0;
        fprintf('[M] Computed snapshot %d/%d (t=%.4f)\n', js, Ns+1, t);
    end
end

% ==========================================================================
% DISCRETE ORDINATE METHOD (DOM)
% ==========================================================================

function sol = domRun(grid, params)
% RUNDOMSOLVER Run discrete ordinate method solver.
%
%   sol = domRun(grid, params) solves the kinetic transport equation
%   using the discrete ordinate method. Returns solution structure with
%   phase space distribution @a sol.F at snapshot times @a sol.t.

    fprintf('[M] Using %d-order FD with %s time integrator\n', ...
        params.xOrder, upper(params.tDisc));

    Nx = params.Nx;
    Nv_moc = params.Nv_moc;
    Nv_dom = params.Nv_dom;
    Nt = grid.Nt;
    Ns = params.Ns;
    si = grid.si;

    % configure initial condition on residual grid (Nr points)
    fprintf('Setting up initial condition...\n');
    [X, V] = computePhaseGrid(grid.x, grid.v_dom, params);
    F0 = computeInitialDistribution(X, V, params);
    F0 = reshape(F0, [Nx, Nv_dom]);

    sol.F_dom = zeros([Nx, Nv_dom, Ns + 1]);
    sol.F_dom(:, :, 1) = F0;
    js = 2;

    % Time stepping loop
    for n = 1:Nt
        t = grid.t(n);

        switch lower(params.tDisc)
            case 'fe'
                F = domStepFe(F0, grid, params);
            case 'heun'
                F = domStepHeun(F0, grid, params);
            case 'ssprk2'
                F = domStepSsprk2(F0, grid, params);
            otherwise
                error('Unknown time integrator: %s', params.tDisc);
        end

        % Save snapshot if current time step is a snapshot time
        if js <= Ns + 1 && n + 1 == si(js)
            sol.F_dom(:, :, js) = F;
            fprintf('[M] t = %.4f (step %d/%d, snapshot %d/%d)\n', ...
                    grid.t(n+1), n, Nt, js, Ns+1);
            js = js + 1;
        end

        F0 = F;
    end

    % Interpolate from Nr points to Nv points for visualization
    sol.F = zeros([Nx, Nv_moc, Ns + 1]);
    for js = 1:(Ns + 1)
        F = sol.F_dom(:, :, js);
        F = rbfInterpolate(F, grid.angle_dom, grid.angle_moc, params);
        sol.F(:, :, js) = F;
    end
end

function F = domStepFe(F0, grid, params)
% DOMSTEPFE Single Forward Euler step for discrete ordinate method.
%
%   F = domStepFe(F0, grid, params) advances the solution by one time
%   step using forward Euler with upwind finite differences.

    Nd = params.Nd;
    Nxp = params.Nxp;
    Nx = params.Nx;
    Nv = params.Nv_dom;
    n = params.xOrder;
    bc = params.bc;
    dt = grid.dt;
    v = grid.v_r;

    dFdx = zeros([Nxp, Nv, Nd]);
    i = [repmat({':'}, 1, Nd), {1, 1}];
    j = [repmat({':'}, 1, Nd), {1}];
    F0 = reshape(F0, [Nxp, Nv]);
    for d = 1:Nd
        h = grid.dx(d);
        i{Nd + 2} = d;
        for k = 1:Nv
            i{Nd + 1} = k;
            j{Nd + 1} = k;
            dFdx(i{:}) = uwfd(F0(j{:}), n, d, h, v(k, d), bc);
        end
    end
    F0 = reshape(F0, [Nx, Nv]);
    dFdx = reshape(dFdx, [Nx, Nv, Nd]);
    dF = sum(reshape(v, [1, Nv, Nd]) .* dFdx, 3);
    F = F0 - dt * dF;
end

function F = domStepHeun(F0, grid, params)
% DOMSTEPHEUN Single Heun (RK2) step for discrete ordinate method.
%
%   F = domStepHeun(F0, grid, params) advances the solution by one time
%   step using Heun's method (2nd order Runge-Kutta).

    Nd = params.Nd;
    Nxp = params.Nxp;
    Nx = params.Nx;
    Nv = params.Nv_dom;
    n = params.xOrder;
    bc = params.bc;
    dt = grid.dt;
    v = grid.v_dom;

    % Stage 1: Compute dF/dt at F0
    dFdx1 = zeros([Nxp, Nv, Nd]);
    i = [repmat({':'}, 1, Nd), {1, 1}];
    j = [repmat({':'}, 1, Nd), {1}];
    F0 = reshape(F0, [Nxp, Nv]);
    for d = 1:Nd
        h = grid.dx(d);
        i{Nd + 2} = d;
        for k = 1:Nv
            i{Nd + 1} = k;
            j{Nd + 1} = k;
            dFdx1(i{:}) = uwfd(F0(j{:}), n, d, h, v(k, d), bc);
        end
    end
    F0 = reshape(F0, [Nx, Nv]);
    dFdx1 = reshape(dFdx1, [Nx, Nv, Nd]);
    dF1 = sum(reshape(v, [1, Nv, Nd]) .* dFdx1, 3);
    F1 = F0 - dt * dF1;

    % Stage 2: Compute dF/dt at F_star
    dFdx2 = zeros([Nxp, Nv, Nd]);
    i = [repmat({':'}, 1, Nd), {1, 1}];
    j = [repmat({':'}, 1, Nd), {1}];
    F1 = reshape(F1, [Nxp, Nv]);
    for d = 1:Nd
        h = grid.dx(d);
        i{Nd + 2} = d;
        for k = 1:Nv
            i{Nd + 1} = k;
            j{Nd + 1} = k;
            dFdx2(i{:}) = uwfd(F1(j{:}), n, d, h, v(k, d), bc);
        end
    end
    dFdx2 = reshape(dFdx2, [Nx, Nv, Nd]);
    dF2 = sum(reshape(v, [1, Nv, Nd]) .* dFdx2, 3);
    F = F0 - dt/2 * (dF1 + dF2);
end

function F = domStepSsprk2(F0, grid, params)
% DOMSTEPSSPRK2 Single SSPRK2 step for discrete ordinate method.
%
%   F = domStepSsprk2(F0, grid, params) advances the solution by one time
%   step using the strong stability preserving Runge-Kutta 2 method.

    Nd = params.Nd;
    Nxp = params.Nxp;
    Nx = params.Nx;
    Nv = params.Nv_dom;
    n = params.xOrder;
    bc = params.bc;
    dt = grid.dt;
    v = grid.v_dom;

    % Stage 1
    dFdx0 = zeros([Nxp, Nv, Nd]);
    i = [repmat({':'}, 1, Nd), {1, 1}];
    j = [repmat({':'}, 1, Nd), {1}];
    F0_reshaped = reshape(F0, [Nxp, Nv]);
    for d = 1:Nd
        h = grid.dx(d);
        i{Nd + 2} = d;
        for k = 1:Nv
            i{Nd + 1} = k;
            j{Nd + 1} = k;
            dFdx0(i{:}) = uwfd(F0_reshaped(j{:}), n, d, h, v(k, d), bc);
        end
    end
    F0 = reshape(F0, [Nx, Nv]);
    dFdx0 = reshape(dFdx0, [Nx, Nv, Nd]);
    dF0 = sum(reshape(v, [1, Nv, Nd]) .* dFdx0, 3);
    F1 = F0 - dt * dF0;

    % Stage 2
    dFdx1 = zeros([Nxp, Nv, Nd]);
    i = [repmat({':'}, 1, Nd), {1, 1}];
    j = [repmat({':'}, 1, Nd), {1}];
    F1_reshaped = reshape(F1, [Nxp, Nv]);
    for d = 1:Nd
        h = grid.dx(d);
        i{Nd + 2} = d;
        for k = 1:Nv
            i{Nd + 1} = k;
            j{Nd + 1} = k;
            dFdx1(i{:}) = uwfd(F1_reshaped(j{:}), n, d, h, v(k, d), bc);
        end
    end
    dFdx1 = reshape(dFdx1, [Nx, Nv, Nd]);
    dF1 = sum(reshape(v, [1, Nv, Nd]) .* dFdx1, 3);
    F = 0.5 * F0 + 0.5 * (F1 - dt * dF1);
end

% ==========================================================================
% URM METHOD
% ==========================================================================

function basis = configureBasis(grid, params)
% configureBASIS Build macro and residual basis functions.
%
%   Constructs:
%   - E: Fourier macro basis (orthonormalized)
%   - Quadrature weights W

    Nu = params.Nu;
    Nr = params.Nr;

    % Macro part: Fourier basis
    if Nu > 0
        basis.Eu = computeFourierBasis(grid.angle_u, params);
    else
        % No macro basis
        basis.Eu = [];
        fprintf('[M] Macro basis: DISABLED (Nu = 0)\n');
    end

    % Residual part
    if Nr > 0
        basis.Er = computeFourierBasis(grid.angle_r, params);
    else
        % No residual basis
        basis.Er = [];
        fprintf('[M] Residual basis: DISABLED (Nr = 0)\n');
    end

    % Store dimensions
    params.Nu = Nu;
    params.Nr = Nr;
    basis.Nz = Nu + Nr;  % Total dimension
end

function A = configureFlux(basis, grid, params)
% CONFIGUREFLUX Assemble flux matrices for URM method.
%
%   For URM method with Z = [U; R], the flux matrix has 2x2 block structure:
%   A^α(1,1) = E' * W * V^α * E
%   A^α(1,2) = E' * W * V^α (coupling U with R)
%   A^α(2,1) = (I - E * E' * W) * V^α * E (coupling R with U)
%   A^α(2,2) = (I - E * E' * W) * V^α (R self-coupling)

    Nd = params.Nd;
    Nu = params.Nu;
    Nr = params.Nr;
    Nz = Nu + Nr;

    % Selector matrix
    k = ceil((Nu - 1)/2);
    l = [0, repelem(1:k, 2)];
    L = l(1:Nu);
    S = ones(1, Nu);
    S(L < max(L)) = 0;
    S = diag(S);

    % Flux matrix
    A = zeros(Nz, Nz, Nd);

    for i = 1:Nd
        % Block (1,1): Nu × Nu (macro self-coupling)
        if Nu > 0
            V = diag(grid.v_u(:, i));
            W = diag(grid.w_u);
            E = basis.Eu;
            A(1:Nu, 1:Nu, i) = E' * W * V * E;
        end

        % Block (1,2): Nu × Nr (macro-residual coupling)
        if Nu > 0 && Nr > 0
            V = diag(grid.v_r(:, i));
            W = diag(grid.w_r);
            E = basis.Er;
            A(1:Nu, Nu+(1:Nr), i) = S * E' * W * V;
        end

        % Block (2,1): Nr × Nu (residual-macro coupling)
        if Nr > 0 && Nu > 0
            V = diag(grid.v_r(:, i));
            W = diag(grid.w_r);
            E = basis.Er;
            I = eye(Nr);
            P = E * (E' * W);
            A(Nu+(1:Nr), 1:Nu, i) = (I - P) * V * E;
        end

        % Block (2,2): Nr × Nr (residual self-coupling)
        if Nr > 0
            V = diag(grid.v_r(:, i));
            if Nu > 0
                W = diag(grid.w_r);
                E = basis.Er;
                I = eye(Nr);
                P = E * S * (E' * W);
                A(Nu+(1:Nr), Nu+(1:Nr), i) = (I - P) * V;
            else
                % When Nu = 0, use pure upwind (same as DOM)
                A(Nu+(1:Nr), Nu+(1:Nr), i) = V;
            end
        end
    end

    fprintf('[M] Flux matrices assembled (URM): %d x %d x %d\n', Nz, Nz, Nd);
end

function sol = urmRun(grid, params)
% RUNURMSOLVER Run URM solver.
%
%   sol = urmRun(grid, params) solves the kinetic transport equation
%   using the URM method (U+R decomposition).
%   Returns solution structure with phase space distribution @a sol.F at
%   snapshot times @a sol.t.

    fprintf('[M] Using %d-order FD with %s time integrator\n', ...
        params.xOrder, upper(params.tDisc));

    Nx = params.Nx;
    Nv_moc = params.Nv_moc;
    Nv_u = params.Nv_u;
    Nu = params.Nu;
    Nr = params.Nr;
    Ns = params.Ns;
    Nt = grid.Nt;
    si = grid.si;

    % configure basis functions
    fprintf('Setting up basis...\n');
    basis = configureBasis(grid, params);
    Nz = basis.Nz;

    % configure flux matrices
    fprintf('Setting up flux matrices...\n');
    A = configureFlux(basis, grid, params);

    % configure initial condition
    fprintf('Setting up initial condition...\n');
    if Nu > 0
        [X, V] = computePhaseGrid(grid.x, grid.v_u, params);
        F0 = computeInitialDistribution(X, V, params);
        F0 = reshape(F0, [Nx, Nv_u]);
        U0 = encodeU(F0, basis, grid, params);
    else
        U0 = [];
    end
    if Nr > 0
        [X, V] = computePhaseGrid(grid.x, grid.v_r, params);
        F0 = computeInitialDistribution(X, V, params);
        F0 = reshape(F0, [Nx, Nr]);
        R0 = encodeROnR(F0, U0, basis, grid, params);
    else
        R0 = [];
    end
    Z0 = [U0, R0];
    Z0 = reshape(Z0, [Nx, Nz]);
    if params.filtering
        Z0 = filtering(Z0, grid, params);
    end

    sol.Z = zeros([Nx, Nz, Ns + 1]);
    sol.Z(:, :, 1) = Z0;
    js = 2;

    % Time stepping loop
    for n = 1:Nt
        t = grid.t(n);

        switch lower(params.tDisc)
            case 'fe'
                Z = urmStepFe(Z0, A, basis, grid, params);
            case 'heun'
                Z = urmStepHeun(Z0, A, basis, grid, params);
            case 'ssprk2'
                Z = urmStepSsprk2(Z0, A, basis, grid, params);
            otherwise
                error('[E] Unknown time integrator: %s', params.tDisc);
        end

        % Save snapshot if current time step is a snapshot time
        if js <= Ns + 1 && n + 1 == si(js)
            sol.Z(:, :, js) = Z;
            fprintf('[M] t = %.4f (step %d/%d, snapshot %d/%d)\n', ...
                    grid.t(n+1), n, Nt, js, Ns+1);
            js = js + 1;
        end

        Z0 = Z;
    end

    % Decode Z to F for all snapshots
    sol.F = zeros([Nx, Nv_moc, Ns + 1]);
    for js = 1:(Ns + 1)
        if Nu > 0
            U = sol.Z(:, 1:Nu, js);
        else
            U = [];
        end
        if Nr > 0
            R = sol.Z(:, Nu+(1:Nr), js);
            R = rbfInterpolate(R, grid.angle_r, grid.angle_moc, params);
        else
            R = [];
        end
        Z = [U, R];
        F = decodeF(Z, basis, grid, params);

        if Nu > 0
            U = encodeUOnMoc(F, basis, grid, params);
        else
            U = [];
        end
        if Nr > 0
            R = encodeR(F, U, basis, grid, params);
        else
            R = [];
        end
        Z = [U, R];
        if params.limiting
            Z = limiting(Z, basis, grid, params);
        end
        F = decodeF(Z, basis, grid, params);

        sol.F(:, :, js) = F;
    end
end

function Z = urmStepFe(Z0, A, basis, grid, params)
% URMSTEPFE Single Forward Euler step for URM method.
%
%   Z = urmStepFe(Z0, A, basis, grid, params) advances the decomposed
%   solution Z = [U; R] by one time step using forward Euler.

    dt = grid.dt;
    dZ = computeSpatialDerivative(Z0, A, basis, grid, params);
    Z = Z0 - dt * dZ;
    if params.filtering
        Z = filtering(Z, grid, params);
    end

    if params.limiting
        Z = limiting(Z, basis, grid, params);
    end
end

function Z = urmStepHeun(Z0, A, basis, grid, params)
% URMSTEPHEUN Single Heun (RK2) step for URM method.
%
%   Z = urmStepHeun(Z0, A, basis, grid, params) advances the decomposed
%   solution Z = [U; R] by one time step using Heun's method.

    dt = grid.dt;

    % Stage 1: Compute dZ/dt at Z0
    dZ1 = computeSpatialDerivative(Z0, A, basis, grid, params);
    Z1 = Z0 - dt * dZ1;
    if params.filtering
        Z1 = filtering(Z1, grid, params);
    end
    if params.limiting
        Z1 = limiting(Z1, basis, grid, params);
    end

    % Stage 2: Compute dZ/dt at Z1
    dZ2 = computeSpatialDerivative(Z1, A, basis, grid, params);
    Z = Z0 - dt/2 * (dZ1 + dZ2);
    if params.filtering
        Z = filtering(Z, grid, params);
    end
    if params.limiting
        Z = limiting(Z, basis, grid, params);
    end
end

function Z = urmStepSsprk2(Z0, A, basis, grid, params)
% URMSTEPSSPRK2 Single SSPRK2 step for URM method.
%
%   Z = urmStepSsprk2(Z0, A, basis, grid, params) advances the decomposed
%   solution Z = [U; R] by one time step using the strong stability
%   preserving Runge-Kutta 2 method (Shu-Osher form).
%
%   SSPRK2 scheme:
%     Z1 = Z0 + dt * L(Z0)
%     Z = (1/2) * Z0 + (1/2) * (Z1 + dt * L(Z1))

    dt = grid.dt;

    % Stage 1
    dZ0 = computeSpatialDerivative(Z0, A, basis, grid, params);
    Z1 = Z0 - dt * dZ0;
    if params.filtering
        Z1 = filtering(Z1, grid, params);
    end
    if params.limiting
        Z1 = limiting(Z1, basis, grid, params);
    end

    % Stage 2
    dZ1 = computeSpatialDerivative(Z1, A, basis, grid, params);
    Z = 0.5 * Z0 + 0.5 * (Z1 - dt * dZ1);
    if params.filtering
        Z = filtering(Z, grid, params);
    end
    if params.limiting
        Z = limiting(Z, basis, grid, params);
    end
end

function dZ = computeSpatialDerivative(Z, A, basis, grid, params)
% COMPUTESPATIALDERIVATIVE Compute spatial derivative for URM method (vectorized).
%
%   dZ = computeSpatialDerivative(Z, A, basis, grid, params)
%
%   For URM with Z = [U; R], uses:
%   - Alternating forward/backward difference for U (macro components)
%     to suppress oscillations. Choice depends on sign of flux matrix A:
%     if A(i,j,d) > 0 use backward, if A(i,j,d) < 0 use forward.
%     This ensures that if A(i,j) uses forward, then A(j,i) uses backward.
%   - Velocity-based upwind for R (residual components at discrete velocities)

    Nd = params.Nd;
    Nxp = params.Nxp;
    Nx = params.Nx;
    Nu = params.Nu;
    Nr = params.Nr;
    n = params.xOrder;
    bc = params.bc;
    v = grid.v_r;  % Use residual velocity grid for upwind
    flux = params.xFlux;

    A_UU = A(1:Nu, 1:Nu, :);
    A_UR = A(1:Nu, Nu+(1:Nr), :);
    A_RU = A(Nu+(1:Nr), 1:Nu, :);
    A_RR = A(Nu+(1:Nr), Nu+(1:Nr), :);
    U = Z(:, 1:Nu);
    R = Z(:, Nu+(1:Nr));
    U = reshape(U, [Nxp, Nu]);
    R = reshape(R, [Nxp, Nr]);

    % Compute all spatial derivatives at once for each dimension
    % dZdx_all will be [Nxp, Nz, Nd, 3] where 3 = {fwd, bwd, ctr}
    if strcmp(flux, 'ctr')
        % Only need central difference
        dUdx_ctr = zeros([Nxp, Nu, Nd]);
        dRdx_ctr = zeros([Nxp, Nr, Nd]);
    else
        % Need forward and backward differences
        dUdx_fwd = zeros([Nxp, Nu, Nd]);
        dUdx_bwd = zeros([Nxp, Nu, Nd]);
        dRdx_fwd = zeros([Nxp, Nr, Nd]);
        dRdx_bwd = zeros([Nxp, Nr, Nd]);
    end
    dRdx_uw = zeros([Nxp, Nr, Nd]);  

    % Pre-compute all derivatives
    i = [repmat({':'}, 1, Nd), {1, 1}];
    j = [repmat({':'}, 1, Nd), {1}];
    for d = 1:Nd
        h = grid.dx(d);
        i{Nd + 2} = d;

        for k = 1:Nu
            i{Nd + 1} = k;
            j{Nd + 1} = k;
            if strcmp(flux, 'ctr')
                dUdx_ctr(i{:}) = ctrfd(U(j{:}), n, d, h, bc);
            else
                dUdx_fwd(i{:}) = fwfd(U(j{:}), n, d, h, bc);
                dUdx_bwd(i{:}) = bwfd(U(j{:}), n, d, h, bc);
            end
        end

        for k = 1:Nr
            i{Nd + 1} = k;
            j{Nd + 1} = k;
            if strcmp(flux, 'ctr')
                dRdx_ctr(i{:}) = ctrfd(R(j{:}), n, d, h, bc);
            else
                dRdx_fwd(i{:}) = fwfd(R(j{:}), n, d, h, bc);
                dRdx_bwd(i{:}) = bwfd(R(j{:}), n, d, h, bc);
            end
        end

        for k = 1:Nr
            i{Nd + 1} = k;
            j{Nd + 1} = k;
            dRdx_uw(i{:}) = uwfd(R(j{:}), n, d, h, v(k, d), bc);
        end
    end

    if strcmp(flux, 'ctr')
        dUdx_ctr = reshape(dUdx_ctr, [Nx, Nu, Nd]);
        dRdx_ctr = reshape(dRdx_ctr, [Nx, Nr, Nd]);
    else
        dUdx_fwd = reshape(dUdx_fwd, [Nx, Nu, Nd]);
        dUdx_bwd = reshape(dUdx_bwd, [Nx, Nu, Nd]);
        dRdx_fwd = reshape(dRdx_fwd, [Nx, Nr, Nd]);
        dRdx_bwd = reshape(dRdx_bwd, [Nx, Nr, Nd]);
    end
    dRdx_uw = reshape(dRdx_uw, [Nx, Nr, Nd]);

    tol = 1e-10;
    dZdx_RR = dRdx_uw;
    switch flux
        case 'ctr'
            dZdx_UR = dRdx_ctr;
            dZdx_RU = dUdx_ctr;
            dZdx_UU = dUdx_ctr;
        case 'alt_fb'
            [k, l] = ndgrid(1:Nu, 1:Nu);
            S = 1 * (k < l) + 2 * (k >= l); 
            mask_UU = abs(A_UU) > tol;
            useF = reshape(mask_UU & (S==1), [Nu, Nu, Nd]);
            useB = reshape(mask_UU & (S==2), [Nu, Nu, Nd]);
            dZdx_UU = zeros([Nx, Nu, Nd]);
            for d = 1:Nd
                [~, l] = find(useF(:, :, d));
                dZdx_UU(:, l, d) = dUdx_fwd(:, l, d);
                [~, l] = find(useB(:, :, d));
                dZdx_UU(:, l, d) = dUdx_bwd(:, l, d);
            end

            mask_UR = abs(A_UR) > tol;
            useF = reshape(mask_UR, [Nu, Nr, Nd]);
            dZdx_UR = zeros([Nx, Nr, Nd]);
            for d = 1:Nd
                [~, l] = find(useF(:, :, d));
                dZdx_UR(:, l, d) = dRdx_fwd(:, l, d);
            end

            mask_RU = abs(A_RU) > tol;
            useB = reshape(mask_RU, [Nr, Nu, Nd]);
            dZdx_RU = zeros([Nx, Nu, Nd]);
            for d = 1:Nd
                [~, l] = find(useB(:, :, d));
                dZdx_RU(:, l, d) = dUdx_bwd(:, l, d);
            end
        case 'alt_bf'
            [k, l] = ndgrid(1:Nu, 1:Nu);
            S = 2 * (k < l) + 1 * (k >= l); 
            mask_UU = abs(A_UU) > tol;
            useF = reshape(mask_UU & (S==1), [Nu, Nu, Nd]);
            useB = reshape(mask_UU & (S==2), [Nu, Nu, Nd]);
            dZdx_UU = zeros([Nx, Nu, Nd]);
            for d = 1:Nd
                [~, l] = find(useF(:, :, d));
                dZdx_UU(:, l, d) = dUdx_fwd(:, l, d);
                [~, l] = find(useB(:, :, d));
                dZdx_UU(:, l, d) = dUdx_bwd(:, l, d);
            end

            mask_UR = abs(A_UR) > tol;
            useB = reshape(mask_UR, [Nu, Nr, Nd]);
            dZdx_UR = zeros([Nx, Nr, Nd]);
            for d = 1:Nd
                [~, l] = find(useB(:, :, d));
                dZdx_UR(:, l, d) = dRdx_bwd(:, l, d);
            end

            mask_RU = abs(A_RU) > tol;
            useF = reshape(mask_RU, [Nr, Nu, Nd]);
            dZdx_RU = zeros([Nx, Nu, Nd]);
            for d = 1:Nd
                [~, l] = find(useF(:, :, d));
                dZdx_RU(:, l, d) = dUdx_fwd(:, l, d);
            end
    end

    A_UU = reshape(A_UU, [1, Nu, Nu, Nd]);
    dZdx_UU = reshape(dZdx_UU, [Nx, 1, Nu, Nd]);
    dZ_UU = sum(A_UU .* dZdx_UU, [3, 4]);

    A_UR = reshape(A_UR, [1, Nu, Nr, Nd]);
    dZdx_UR = reshape(dZdx_UR, [Nx, 1, Nr, Nd]);
    dZ_UR = sum(A_UR .* dZdx_UR, [3, 4]);

    A_RU = reshape(A_RU, [1, Nr, Nu, Nd]);
    dZdx_RU = reshape(dZdx_RU, [Nx, 1, Nu, Nd]);
    dZ_RU = sum(A_RU .* dZdx_RU, [3, 4]);

    A_RR = reshape(A_RR, [1, Nr, Nr, Nd]);
    dZdx_RR = reshape(dZdx_RR, [Nx, 1, Nr, Nd]);
    dZ_RR = sum(A_RR .* dZdx_RR, [3, 4]);

    dZ_U = dZ_UU + dZ_UR;
    dZ_R = dZ_RU + dZ_RR;

    if strcmp(flux, 'ctr')
        if Nu > 0
            U = reshape(U, [Nx, Nu]);
            U0 = reshape(U(:, 1), Nxp);
            D_U0 = computeArtificialDissipation(U0, grid.dx, params);
            D_U0 = reshape(D_U0, [Nx, 1]);
            dZ_U(:, 1) = dZ_U(:, 1) + D_U0;
        end
    end

    dZ = [dZ_U, dZ_R];
end

% ========================================================================
% COMPARISON AND ANALYSIS
% ========================================================================
function compareResults(sol_moc, sol_dom, sol_urm, grid, params)
% COMPARERESULTS Compare MoC vs DoM vs URM methods.
%
%   compareResults(sol_moc, sol_dom, sol_urm, grid, params) performs
%   comprehensive comparison of the three solvers including: (1) mass
%   conservation, (2) L2 errors, (3) 2D density snapshots, (4) 1D slices
%   along x=0, y=0, and y=x, and (5) angular distribution at selected
%   spatial locations.

    Nxp = params.Nxp;
    Ns = params.Ns + 1;
    w = grid.w_moc;

    % Extract distributions
    f_moc = sol_moc.F;
    f_dom = sol_dom.F;
    f_urm = sol_urm.F;

    % Compute densities by integrating over velocity
    rho_moc = zeros([Nxp, Ns]);
    rho_dom = zeros([Nxp, Ns]);
    rho_urm = zeros([Nxp, Ns]);

    for n = 1:Ns
        rho_moc(:, :, n) = reshape(sum(f_moc(:, :, n) .* w(:).', 2), Nxp);
        rho_dom(:, :, n) = reshape(sum(f_dom(:, :, n) .* w(:).', 2), Nxp);
        rho_urm(:, :, n) = reshape(sum(f_urm(:, :, n) .* w(:).', 2), Nxp);
    end

    % Compute masses
    mass_moc = zeros(Ns, 1);
    mass_dom = zeros(Ns, 1);
    mass_urm = zeros(Ns, 1);

    for n = 1:Ns
        mass_moc(n) = trapzNd(rho_moc(:, :, n), grid.dx);
        mass_dom(n) = trapzNd(rho_dom(:, :, n), grid.dx);
        mass_urm(n) = trapzNd(rho_urm(:, :, n), grid.dx);
    end

    % Print statistics
    fprintf('\n[M] Mass Conservation Comparison:\n');
    fprintf('    MoC (Exact):  Initial=%.6f, Final=%.6f, Change=%.2e\n', ...
            mass_moc(1), mass_moc(end), abs(mass_moc(end)-mass_moc(1))/mass_moc(1));
    fprintf('    DoM:          Initial=%.6f, Final=%.6f, Change=%.2e\n', ...
            mass_dom(1), mass_dom(end), abs(mass_dom(end)-mass_dom(1))/mass_dom(1));
    fprintf('    URM:          Initial=%.6f, Final=%.6f, Change=%.2e\n', ...
            mass_urm(1), mass_urm(end), abs(mass_urm(end)-mass_urm(1))/mass_urm(1));

    % Compute L2 errors relative to MoC solution
    l2_errors_dom = zeros(Ns, 1);
    l2_errors_urm = zeros(Ns, 1);

    for n = 1:Ns
        diff_dom = rho_dom(:, :, n) - rho_moc(:, :, n);
        diff_urm = rho_urm(:, :, n) - rho_moc(:, :, n);
        l2_errors_dom(n) = sqrt(trapzNd(diff_dom.^2, grid.dx));
        l2_errors_urm(n) = sqrt(trapzNd(diff_urm.^2, grid.dx));
    end

    fprintf('\n[M] L2 Errors vs MoC:\n');
    fprintf('    DoM: Initial=%.2e, Final=%.2e\n', ...
            l2_errors_dom(1), l2_errors_dom(end));
    fprintf('    URM: Initial=%.2e, Final=%.2e\n', ...
            l2_errors_urm(1), l2_errors_urm(end));

    % Visualize 2D density snapshots
    plot2dSnapshots(rho_moc, rho_dom, rho_urm, grid.s, grid);

    % Visualize 1D slices
    plot1dSlices(rho_moc, rho_dom, rho_urm, grid.s, grid);

    % Visualize angular distribution at center and initial time
    plotAngularDistribution(f_moc, f_dom, f_urm, grid, params);
end

function plot2dSnapshots(rho_moc, rho_dom, rho_urm, times, grid)
% PLOT2DSNAPSHOTS Visualize density snapshots for all three methods.
%
%   plot2dSnapshots(rho_moc, rho_dom, rho_urm, times, grid) creates a 4-row
%   subplot visualization comparing MoC (exact), DoM, URM methods, and
%   displays the DoM and URM errors relative to MoC at selected time
%   snapshots.

    Ns = length(times);
    time_indices = round(linspace(1, Ns, min(4, Ns)));

    figure(1);
    for idx = 1:length(time_indices)
        n = time_indices(idx);

        % Color limits based on MoC solution
        cmin = min(rho_moc(:, :, n), [], 'all');
        cmax = max(rho_moc(:, :, n), [], 'all');

        % MoC solution (exact)
        subplot(4, length(time_indices), idx);
        imagesc(grid.x{1}, grid.x{2}, rho_moc(:, :, n).');
        colorbar; axis equal tight; clim([cmin,cmax]);
        title(sprintf('MoC t=%.3f', times(n)));
        xlabel('x'); ylabel('y');

        % DoM solution
        subplot(4, length(time_indices), idx + length(time_indices));
        imagesc(grid.x{1}, grid.x{2}, rho_dom(:, :, n).');
        colorbar; axis equal tight; clim([cmin,cmax]);
        title(sprintf('DoM t=%.3f', times(n)));
        xlabel('x'); ylabel('y');

        % URM solution
        subplot(4, length(time_indices), idx + 2*length(time_indices));
        imagesc(grid.x{1}, grid.x{2}, rho_urm(:, :, n).');
        colorbar; axis equal tight; clim([cmin,cmax]);
        title(sprintf('URM t=%.3f', times(n)));
        xlabel('x'); ylabel('y');

        % Error: URM vs MoC
        subplot(4, length(time_indices), idx + 3*length(time_indices));
        diff = rho_urm(:, :, n) - rho_moc(:, :, n);
        imagesc(grid.x{1}, grid.x{2}, diff.');
        colorbar; axis equal tight;
        title(sprintf('URM Error t=%.3f', times(n)));
        xlabel('x'); ylabel('y');
    end
    sgtitle('Density Comparison: MoC vs DoM vs URM');

end

function plot1dSlices(rho_moc, rho_dom, rho_urm, times, grid)
% PLOT1DSLICES Visualize 1D density slices along x=0, y=0, and y=x.
%
%   plot1dSlices(rho_moc, rho_dom, rho_urm, times, grid) creates a 3x4
%   subplot visualization showing 1D slices of the density at multiple
%   time snapshots along three lines: (1) y=0 (horizontal), (2) x=0
%   (vertical), and (3) y=x (diagonal).

    Ns = length(times);
    time_indices = round(linspace(1, Ns, min(4, Ns)));

    % Get grid coordinates
    x = grid.x{1};
    y = grid.x{2};
    Nx = length(x);
    Ny = length(y);

    % Find center indices
    [~, idx_x0] = min(abs(x));
    [~, idx_y0] = min(abs(y));

    figure(2);

    for idx = 1:length(time_indices)
        n = time_indices(idx);

        % Slice along y=0 (horizontal)
        subplot(3, length(time_indices), idx);
        hold on;
        plot(x, rho_moc(:, idx_y0, n), 'k-', LineWidth=2);
        plot(x, rho_dom(:, idx_y0, n), 'b--', LineWidth=1.5);
        plot(x, rho_urm(:, idx_y0, n), 'r:', LineWidth=1.5);
        hold off;
        xlabel('x'); ylabel('\rho(x, 0)');
        title(sprintf('y=0, t=%.3f', times(n)));
        if idx == 1
            legend('MoC', 'DoM', 'URM', Location='best');
        end

        % Slice along x=0 (vertical)
        subplot(3, length(time_indices), idx + length(time_indices));
        hold on;
        plot(y, rho_moc(idx_x0, :, n), 'k-', LineWidth=2);
        plot(y, rho_dom(idx_x0, :, n), 'b--', LineWidth=1.5);
        plot(y, rho_urm(idx_x0, :, n), 'r:', LineWidth=1.5);
        hold off;
        xlabel('y'); ylabel('\rho(0, y)');
        title(sprintf('x=0, t=%.3f', times(n)));

        % Slice along y=x (diagonal)
        subplot(3, length(time_indices), idx + 2*length(time_indices));
        diag_moc = zeros(1, min(Nx, Ny));
        diag_dom = zeros(1, min(Nx, Ny));
        diag_urm = zeros(1, min(Nx, Ny));
        diag_coord = zeros(1, min(Nx, Ny));
        for i = 1:min(Nx, Ny)
            diag_moc(i) = rho_moc(i, i, n);
            diag_dom(i) = rho_dom(i, i, n);
            diag_urm(i) = rho_urm(i, i, n);
            diag_coord(i) = x(i);
        end
        hold on;
        plot(diag_coord, diag_moc, 'k-', LineWidth=2);
        plot(diag_coord, diag_dom, 'b--', LineWidth=1.5);
        plot(diag_coord, diag_urm, 'r:', LineWidth=1.5);
        hold off;
        xlabel('s'); ylabel('\rho(s, s)');
        title(sprintf('y=x, t=%.3f', times(n)));
    end

    sgtitle('1D Density Slices: MoC vs DoM vs URM');
end

function plotAngularDistribution(f_moc, f_dom, f_urm, grid, params)
% PLOTANGULARDISTRIBUTION Visualize angular distribution at center point.
%
%   plotAngularDistribution(f_moc, f_dom, f_urm, grid, params) creates a
%   single plot showing the angular distribution f(0, 0, θ) at the initial
%   time (t=0) for the center spatial location.

    Nxp = params.Nxp;
    n_initial = 1;

    % Get grid coordinates
    x = grid.x{1};
    y = grid.x{2};
    angle = grid.angle_moc;

    % Find center index
    [~, idx_x0] = min(abs(x));
    [~, idx_y0] = min(abs(y));
    idx_center = sub2ind(Nxp, idx_x0, idx_y0);

    figure(3);

    % Angular distribution at center and initial time
    hold on;
    plot(angle, f_moc(idx_center, :, n_initial), 'k-', LineWidth=2, ...
        DisplayName='MoC');
    plot(angle, f_dom(idx_center, :, n_initial), 'b--', LineWidth=1.5, ...
        DisplayName='DoM');
    plot(angle, f_urm(idx_center, :, n_initial), 'r:', LineWidth=1.5, ...
        DisplayName='URM');
    hold off;
    xlabel('\theta'); ylabel('f(0, 0, \theta)');
    title(sprintf('Angular distribution at center (0,0), t=%.3f', ...
        grid.s(n_initial)));
    legend(Location='best');
end

% ========================================================================
% HELPER FUNCTIONS
% ========================================================================

function [X, V] = computePhaseGrid(x, v, params)
% COMPUTEPHASEGRID Build phase space grid from spatial and velocity grids.
%
%   [X, V] = computePhaseGrid(x, v, params) constructs phase space coordinate
%   matrices @a X and @a V by taking the tensor product of spatial grid @a x
%   and velocity grid @a v.

    Nd = params.Nd;
    Nx = params.Nx;
    Nv = size(v, 1);
    [x{:}] = ndgrid(x{:});
    x = reshape(cat(Nd+1, x{:}), [Nx, Nd]);
    X = kron(ones(Nv, 1), x);
    V = kron(v, ones(Nx, 1));
end

function F0 = computeInitialDistribution(X, V, params)
% COMPUTEINITIALDISTRIBUTION Generate initial condition in phase space.
%
%   F0 = computeInitialDistribution(X, V, params) generates the initial
%   phase space distribution based on @a params.ic.
%
%   Parameters:
%   - X: Spatial coordinates matrix [Nx*Nv x Nd]
%   - V: Velocity grid matrix [Nx*Nv x Nd]
%   - params: Structure containing problem parameters
%
%   Supported initial condition types:
%   - 'gaussian': Isotropic Gaussian f = exp(-r²/(2σ²))/(2πσ²)
%   - 'ring': Ring-shaped distribution
%   - 'two_peaks': Two separated Gaussians
%   - 'beam': Gaussian in space and velocity (collimated beam)

    switch lower(params.ic)
        case 'gaussian'
            F0 = computeIsotropicGaussian(X, V, params);

        case 'ring'
            F0 = computeIsotropicRing(X, V, params);

        case 'two_peaks'
            F0 = computeIsotropicTwoPeaks(X, V, params);

        case 'beam'
            F0 = computeBeam(X, V, params);

        otherwise
            error('Unknown initial condition type: %s', params.ic);
    end
end

function F0 = computeIsotropicGaussian(X, V, params)
% COMPUTEISOTROPICGAUSSIAN Generate isotropic Gaussian initial distribution.
%
%   F0 = computeIsotropicGaussian(X, V, params) creates a spatially isotropic
%   Gaussian distribution centered at origin with width @a params.sigma.

    Nd = params.Nd;

    x0 = zeros(1, Nd);
    sigma = params.sigma;

    R = sqrt(sum((X - x0(:).').^2, 2));
    rho = exp(-R.^2 / (2*sigma^2));
    rho = rho / (sqrt(2*pi)*sigma)^Nd;
    F0 = rho / computeSphereArea(Nd);
end

function F0 = computeIsotropicRing(X, V, params)
% COMPUTEISOTROPICRING Generate ring-shaped initial distribution.
%
%   F0 = computeIsotropicRing(X, V, params) creates a spatially isotropic
%   ring distribution centered at @a params.x0 with radius @a params.r0.

    Nd = params.Nd;
    x0 = params.x0;
    r0 = params.r0;
    sigma = params.sigma;

    R = sqrt(sum((X - x0(:).').^2, 2));
    rho = exp(-(R - r0).^2 / (2*sigma^2));
    rho = rho / (sqrt(2*pi)*sigma)^Nd;
    F0 = rho / computeSphereArea(Nd);
end

function F0 = computeIsotropicTwoPeaks(X, V, params)
% COMPUTEISOTROPICTWOPEAKS Generate two-peaked initial distribution.
%
%   F0 = computeIsotropicTwoPeaks(X, V, params) creates a spatially isotropic
%   distribution with two Gaussian peaks centered at @a params.x1 and @a params.x2.

    Nd = params.Nd;

    x1 = params.x1;
    x2 = params.x2;
    sigma = params.sigma;

    R1 = sqrt(sum((X - x1(:).').^2, 2));
    R2 = sqrt(sum((X - x2(:).').^2, 2));

    rho = exp(-R1.^2 / (2*sigma^2)) + exp(-R2.^2 / (2*sigma^2));
    rho = rho / (sqrt(2*pi)*sigma)^Nd;
    F0 = rho / computeSphereArea(Nd);
end

function F0 = computeBeam(X, V, params)
% COMPUTEBEAM Generate collimated beam initial distribution.
%
%   F0 = computeBeam(X, V, params) creates a beam distribution that is
%   localized in both space and velocity. The spatial part is a Gaussian
%   centered at @a params.x0 with width @a params.sigma. The angular part
%   is a Gaussian in angle centered around direction @a params.v0 with
%   angular width @a params.sigma_v.

    Nd = params.Nd;
    x0 = params.x0;
    sigma = params.sigma;
    v0 = params.v0;
    sigma_v = params.sigma_v;

    % Spatial Gaussian density
    R = sqrt(sum((X - x0(:).').^2, 2));
    rho = exp(-R.^2 / (2*sigma^2));
    rho = rho / (sqrt(2*pi)*sigma)^Nd;

    % Compute target angle from v0
    theta0 = atan2(v0(2), v0(1));

    % Compute angle for each velocity vector
    theta = atan2(V(:, 2), V(:, 1));

    % Angular distance on circle (accounting for periodicity)
    dtheta = theta - theta0;
    dtheta = mod(dtheta + pi, 2*pi) - pi;

    % Angular Gaussian distribution (normalized on circle)
    g = exp(-dtheta.^2 / (2*sigma_v^2));
    g = g / (sqrt(2*pi) * sigma_v);

    % Combined distribution
    F0 = rho .* g;
end

function A = computeSphereArea(Nd)
% COMPUTESPHEREAREA Compute surface area of unit sphere in Nd dimensions.
%
%   A = computeSphereArea(Nd) returns the surface area of the unit sphere
%   S^{Nd-1} embedded in Nd-dimensional space using the formula:
%   A = 2π^{n/2} / Γ(n/2) where n = Nd - 1.

    n = Nd - 1;

    switch n
        case 0
            A = 2;
        case 1
            A = 2 * pi;
        case 2
            A = 4 * pi;
        otherwise
            A = 2 * pi^((n + 1) / 2) / gamma((n + 1)/2);
    end
end

function E = computeFourierBasis(angle, params)
% COMPUTEFOURIERBASIS Build orthonormalized Fourier basis on S^1.
%
%   E = computeFourierBasis(angle, params) constructs an orthonormal Fourier
%   basis {1, cos(θ), sin(θ), cos(2θ), sin(2θ), ...} on the unit circle,
%   normalized with respect to the L2 inner product with uniform measure.

    Nv = length(angle);
    Nu = params.Nu;

    E = zeros(Nv, Nu);
    E(:, 1) = 1 / sqrt(2*pi);  % Normalized constant

    for k = 1:ceil((Nu-1)/2)
        if 2*k <= Nu
            E(:, 2*k) = sqrt(1/pi) * cos(k*angle);
        end
        if 2*k+1 <= Nu
            E(:, 2*k+1) = sqrt(1/pi) * sin(k*angle);
        end
    end
end

function U = encodeU(F, basis, grid, params)
    Nx = params.Nx;
    Nu = params.Nu;

    U = zeros([Nx, Nu]);
    E = computeFourierBasis(grid.angle_u, params);
    W = diag(grid.w_u);
    U(:, 1:Nu) = F * W * E;
end


function U = encodeUOnMoc(F, basis, grid, params)
    Nx = params.Nx;
    Nu = params.Nu;

    U = zeros([Nx, Nu]);
    E = computeFourierBasis(grid.angle_moc, params);
    W = diag(grid.w_moc);
    U(:, 1:Nu) = F * W * E;
end

function R = encodeR(F, U, basis, grid, params)
    if isempty(U)
        R = F;
    else
        E = computeFourierBasis(grid.angle_moc, params);
        R = F - U * E.';
    end
end

function R = encodeROnR(F, U, basis, grid, params)
    if isempty(U)
        R = F;
    else
        E = computeFourierBasis(grid.angle_r, params);
        R = F - U * E.';
    end
end

function F = decodeF(Z, basis, grid, params)
    Nu = params.Nu;
    Nv = params.Nv_moc;
    Nr = params.Nr;

    if Nu > 0
        U = Z(:, 1:Nu);
        E = computeFourierBasis(grid.angle_moc, params);
        U = U * E';
    else
        U = 0;
    end

    if Nr > 0
        R = Z(:, Nu+(1:Nv));
    else
        R = 0;
    end

    F = U + R;
end

function F = decodeFOnR(Z, basis, grid, params)
    Nu = params.Nu;
    Nr = params.Nr;

    if Nu > 0
        U = Z(:, 1:Nu);
        E = computeFourierBasis(grid.angle_r, params);
        U = U * E';
    else
        U = 0;
    end

    if Nr > 0
        R = Z(:, Nu+(1:Nr));
    else
        R = 0;
    end

    F = U + R;
end

function Y = fdpad(X, n, options)
%FDPAD Pad array with NumPy-like modes and per-side widths.
%
%   Y = fdpad(X, n, options)
%     X : array (numeric or logical), any dimensions
%     n : vector of length 2*ndims(X):
%         [nl(1), nr(1), nl(2), nr(2), ..., nl(nd), nr(nd)]
%     options.mode : 'constant' (default) | 'edge' | 'reflect' | 'symmetric' | 'wrap'
%     options.val : scalar pad value for 'constant' (default 0)

    arguments
        X {mustBeNumericOrLogical}
        n {mustBeVector, mustBeNumeric}
        options.mode {mustBeTextScalar} = "constant"
        options.val {mustBeScalarOrEmpty} = 0
    end

    nd = ndims(X);
    if ~isvector(n) || ~isnumeric(n) || ~all(isfinite(n))
        error('n must be a finite numeric vector.');
    end
    if numel(n) ~= 2*nd
        error('n must have length 2*ndims(X) (= %d).', 2*nd);
    end
    if any(n < 0) || any(n ~= floor(n))
        error('n must contain nonnegative integers.');
    end

    mode = lower(string(options.mode));
    if ~ismember(mode, ["constant","edge","reflect","symmetric","wrap"])
        error('unknown mode. Use constant|edge|reflect|symmetric|wrap.');
    end

    nl = n(1:2:end);
    nr = n(2:2:end);
    sz = size(X);

    % -------- fast path: constant --------
    if mode == "constant"
        if ~(isscalar(options.val) && isnumeric(options.val) && isfinite(options.val))
            error('options.val must be a finite numeric scalar for mode ''constant''.');
        end
        szp = sz;
        for d = 1:nd
            szp(d) = sz(d) + nl(d) + nr(d);
        end

        Y = zeros(szp, 'like', X);
        if islogical(X)
            cval = logical(options.val ~= 0);
        else
            cval = cast(options.val, 'like', X);
        end
        if ~isequal(cval, 0)
            Y(:) = cval;
        end

        subs = cell(1, nd);
        for d = 1:nd
            subs{d} = (nl(d)+1) : (nl(d) + sz(d));
        end
        Y(subs{:}) = X;
        return
    end

    idx = cell(1, nd);
    for d = 1:nd
        szd = sz(d);
        if szd < 1
            error('fdpad: size(X,%d) must be >= 1.', d);
        end
        mL = nl(d); mR = nr(d);

        switch mode
            case "edge"
                left  = ones(1, mL);
                mid   = 1:szd;
                right = repmat(szd, 1, mR);
                idx{d} = [left, mid, right];

            case "wrap"
                % inline modular wrap (handles mL/mR = 0 too)
                left  = mod((-mL+1:0) - 1, szd) + 1;
                mid   = 1:szd;
                right = mod((1:mR) - 1, szd) + 1;
                idx{d} = [left, mid, right];

            case "reflect"  % no edge repetition
                if (mL > szd-1) || (mR > szd-1)
                    error(['fdpad: for ''reflect'', nl(%d) and nr(%d) must be <= size(X,%d)-1 ' ...
                                   '(size=%d).'], d, d, d, szd);
                end
                % left: [2,3,...,mL+1] reversed; right: [szd-1, ..., szd-mR]
                left  = (2:(mL+1));  left  = left(end:-1:1);
                mid   = 1:szd;
                right = (szd-1):-1:(szd-mR);
                idx{d} = [left, mid, right];

            case "symmetric" % includes edge
                if (mL > szd) || (mR > szd)
                    error(['fdpad: for ''symmetric'', nl(%d) and nr(%d) must be <= size(X,%d) ' ...
                                   '(size=%d).'], d, d, d, szd);
                end
                % left: [1,2,...,mL] reversed; right: [szd, ..., szd-mR+1]
                left  = (1:mL);      left  = left(end:-1:1);
                mid   = 1:szd;
                right =  szd:-1:(szd-mR+1);
                idx{d} = [left, mid, right];
        end
    end

    Y = X(idx{:});
end

function w = fdcoeffs(nodes)
% === helper: Fornberg-style weights for the first derivative at x0 on 'nodes' ===
% Return weights w such that f'(0) ≈ sum_j w(j) * f(nodes(j))
% nodes: vector of distinct grid positions (e.g., 0:-1:-n or 0:1:n)
    x = nodes(:).';    % 1×(n+1)
    x0 = 0;
    m  = 1;            % derivative order
    N  = numel(x);
    c  = zeros(N, m+1);
    c1 = 1;
    c4 = x(1) - x0;
    c(1,1) = 1;
    for i = 2:N
        mn = min(i, m+1);
        c2 = 1; c5 = c4;
        c4 = x(i) - x0;
        for j = 1:i-1
            c3 = x(i) - x(j);
            c2 = c2 * c3;
            if j == i-1
                for k = mn:-1:2
                    c(i,k) = c1 * ((k-1)*c(i-1,k-1) - c5*c(i-1,k)) / c2;
                end
                c(i,1) = -c1 * c5 * c(i-1,1) / c2;
            end
            for k = mn:-1:2
                c(j,k) = (c4*c(j,k) - (k-1)*c(j,k-1)) / c3;
            end
            c(j,1) = c4 * c(j,1) / c3;
        end
        c1 = c2;
    end
    w = c(:, m+1).';  % row vector of length N
end

function Y = uwfd(X, n, d, h, v, bc)
%UDFD Upwind finite-difference first derivative along dimension d.
%
%   Y = uwfd(X, n, d, h, v, bc)
%     X : input array
%     n : accuracy order (positive integer, e.g., 1,2,3,...)
%     d : dimension to differentiate along (1-based)
%     h : grid spacing (>0)
%     v : velocity component (scalar). Sign selects upwind direction:
%           v > 0  => backward stencil
%           v < 0  => forward  stencil
%           v = 0  => returns zeros
%     bc: boundary condition: 'zero' | 'periodic' | 'reflect'
%
% Notes
% - Uses fdpad(X, padvec, options) where padvec = [nl(1),nr(1),...,nl(D),nr(D)].
% - For arbitrary order n, coefficients are generated by a Fornberg-style routine
%   on a unit-spaced stencil, then scaled by 1/h.

    arguments
        X {mustBeNumeric}
        n (1,1) {mustBeInteger, mustBePositive}
        d (1,1) {mustBeInteger, mustBePositive}
        h {mustBeNumeric, mustBeFinite, mustBePositive}
        v (1,1) {mustBeNumeric, mustBeFinite}
        bc {mustBeTextScalar}
    end

    D = ndims(X);
    if d < 1 || d > D
        error('uwfd: d must be between 1 and ndims(X) (= %d).', D);
    end

    % Trivial case
    if v == 0
        Y = zeros(size(X), 'like', X);
        return
    end

    % Map boundary condition to fdpad options
    switch lower(bc)
        case 'zero'
            fd_mode = "constant"; pad_val = 0;
        case 'periodic'
            fd_mode = "wrap";     pad_val = []; %#ok<NASGU>
        case 'reflect'
            fd_mode = "edge";     pad_val = [];
        otherwise
            error('uwfd: bc must be ''zero'', ''periodic'', or ''reflect''.');
    end

    % Choose upwind direction and stencil nodes
    if v > 0
        % backward stencil uses nodes: 0, -1, -2, ..., -n (length n+1)
        nodes = 0:-1:-n;
        needL = n; needR = 0;
    else
        % forward stencil uses nodes: 0, +1, +2, ..., +n
        nodes = 0:1:n;
        needL = 0; needR = n;
    end

    % Build per-dimension pad widths [nl(1),nr(1),...,nl(D),nr(D)]
    pad = zeros(1, 2*D);
    % For periodic BC we can safely pad on both sides by n (simplifies wrap)
    if fd_mode == "wrap"
        pad(2*d-1) = n;  % left
        pad(2*d)   = n;  % right
    else
        pad(2*d-1) = needL;
        pad(2*d)   = needR;
    end

    % Prepare fdpad options
    opts = struct('mode', fd_mode);
    if fd_mode == "constant"
        opts.val = 0;
    end
    opts = namedargs2cell(opts);

    % Pad
    Xp = fdpad(X, pad, opts{:});

    % Coefficients for first derivative at x0 = 0 on the chosen nodes
    w = fdcoeffs(nodes);   % length n+1, sums to 0
    % Scale by 1/h for physical spacing
    w = w / h;

    % Build vectorized slices and accumulate
    sz = size(X);
    Y  = zeros(sz, 'like', X);

    % Base indices along d for the interior in the padded array
    base = (1:sz(d)) + pad(2*d-1);  % shift by left pad on axis d

    % Pre-build ":" subscripts
    S = repmat({':'}, 1, D);

    % Accumulate contributions for each stencil node
    for k = 1:numel(nodes)
        offset = nodes(k);   % negative for backward, positive for forward
        Sk = S;
        Sk{d} = base + offset;
        Y = Y + w(k) * Xp(Sk{:});
    end
end

function Y = ctrfd(X, n, d, h, bc)
% CTRFD Central finite-difference first derivative along dimension d.
%
%   Y = ctrfd(X, n, d, h, bc)
%     X    : input array (any dimensions)
%     n    : accuracy order (positive even integer, e.g., 2,4,...)
%     d    : dimension to differentiate along (1-based)
%     h    : grid spacing (>0)
%     bc   : boundary type: 'zero' | 'periodic' | 'reflect'
%
% Notes
% - Uses fdpad with per-side widths [.., n_left(d)=n, n_right(d)=n, ..].
% - Fully vectorized; only loops over stencil taps.

    arguments
        X {mustBeNumeric}
        n (1,1) {mustBeInteger, mustBePositive}
        d (1,1) {mustBeInteger, mustBePositive}
        h {mustBeNumeric, mustBeFinite, mustBePositive}
        bc {mustBeMember(bc, {'zero','periodic','reflect'})}
    end

    n = ceil(n/2);
    D = ndims(X);
    pad = zeros(1, 2*D);
    pad(2*d-1) = n;
    pad(2*d)   = n;

    switch bc
        case "zero"
            pad_opts = struct('mode',"constant", 'val', 0);
        case "periodic"
            pad_opts = struct('mode',"wrap");
        case "reflect"
            pad_opts = struct('mode',"edge");
    end
    pad_opts = namedargs2cell(pad_opts);

    Xp = fdpad(X, pad, pad_opts{:});

    nodes = -n:n;
    w = fdcoeffs(nodes) / h;

    sz = size(X);
    Y  = zeros(sz, 'like', X);

    base = (1:sz(d)) + pad(2*d-1);
    S = repmat({':'}, 1, D);

    for k = 1:numel(nodes)
        Sk = S;
        Sk{d} = base + nodes(k);
        Y = Y + w(k) * Xp(Sk{:});
    end
end

function Y = fwfd(X, n, d, h, bc)
% FWFD Forward finite-difference first derivative along dimension d.
%
%   Y = fwfd(X, n, d, h, bc) computes first derivative using forward
%   difference stencil (uses points at i, i+1, i+2, ..., i+n).
%
%     X    : input array (any dimensions)
%     n    : accuracy order (positive integer, e.g., 1,2,3,...)
%     d    : dimension to differentiate along (1-based)
%     h    : grid spacing (>0)
%     bc   : boundary type: 'zero' | 'periodic' | 'reflect'

    arguments
        X {mustBeNumeric}
        n (1,1) {mustBeInteger, mustBePositive}
        d (1,1) {mustBeInteger, mustBePositive}
        h {mustBeNumeric, mustBeFinite, mustBePositive}
        bc {mustBeMember(bc, {'zero','periodic','reflect'})}
    end

    D = ndims(X);
    pad = zeros(1, 2*D);
    pad(2*d-1) = 0;     % no left padding needed
    pad(2*d)   = n;     % right padding for forward stencil

    switch bc
        case "zero"
            pad_opts = struct('mode',"constant", 'val', 0);
        case "periodic"
            pad_opts = struct('mode',"wrap");
        case "reflect"
            pad_opts = struct('mode',"edge");
    end
    pad_opts = namedargs2cell(pad_opts);

    Xp = fdpad(X, pad, pad_opts{:});

    % Forward stencil nodes: 0, +1, +2, ..., +n
    nodes = 0:n;
    w = fdcoeffs(nodes) / h;

    sz = size(X);
    Y = zeros(sz, 'like', X);

    base = (1:sz(d)) + pad(2*d-1);
    S = repmat({':'}, 1, D);

    for k = 1:numel(nodes)
        Sk = S;
        Sk{d} = base + nodes(k);
        Y = Y + w(k) * Xp(Sk{:});
    end
end

function Y = bwfd(X, n, d, h, bc)
% BWFD Backward finite-difference first derivative along dimension d.
%
%   Y = bwfd(X, n, d, h, bc) computes first derivative using backward
%   difference stencil (uses points at i, i-1, i-2, ..., i-n).
%
%     X    : input array (any dimensions)
%     n    : accuracy order (positive integer, e.g., 1,2,3,...)
%     d    : dimension to differentiate along (1-based)
%     h    : grid spacing (>0)
%     bc   : boundary type: 'zero' | 'periodic' | 'reflect'

    arguments
        X {mustBeNumeric}
        n (1,1) {mustBeInteger, mustBePositive}
        d (1,1) {mustBeInteger, mustBePositive}
        h {mustBeNumeric, mustBeFinite, mustBePositive}
        bc {mustBeMember(bc, {'zero','periodic','reflect'})}
    end

    D = ndims(X);
    pad = zeros(1, 2*D);
    pad(2*d-1) = n;     % left padding for backward stencil
    pad(2*d)   = 0;     % no right padding needed

    switch bc
        case "zero"
            pad_opts = struct('mode',"constant", 'val', 0);
        case "periodic"
            pad_opts = struct('mode',"wrap");
        case "reflect"
            pad_opts = struct('mode',"edge");
    end
    pad_opts = namedargs2cell(pad_opts);

    Xp = fdpad(X, pad, pad_opts{:});

    % Backward stencil nodes: 0, -1, -2, ..., -n
    nodes = 0:-1:-n;
    w = fdcoeffs(nodes) / h;

    sz = size(X);
    Y = zeros(sz, 'like', X);

    base = (1:sz(d)) + pad(2*d-1);
    S = repmat({':'}, 1, D);

    for k = 1:numel(nodes)
        Sk = S;
        Sk{d} = base + nodes(k);
        Y = Y + w(k) * Xp(Sk{:});
    end
end

function Z = filtering(Z, grid, params)
% FILTERING Apply spectral filter to macro coefficients.
%
%   Z = filtering(Z, params) applies a spectral filter to the macro (Fourier)
%   coefficients to damp high-frequency modes and suppress oscillations.
%
%   Supported filter types:
%   - 'rot': Rational filter, sigma_k = 1/(1 + alpha * k^4)
%            where alpha = g(k_max) * dt^(1 - (k/k_max)^p)
%   - 'exp': Exponential filter, sigma_k = exp(-alpha*(k/k_c)^p)

    Nu = params.Nu;
    Nr = params.Nr;
    tp = params.filterType;
    s = params.filterStrength;
    p = params.filterOrder;
    c = params.filterCutoff;
    dt = grid.dt;

    kmax = ceil((Nu - 1) / 2);
    k = [0, repelem(1:kmax, 2)];
    k = k(1:Nu);
    if kmax == 0
        eta = 0;
    else
        eta = k / kmax;
    end

    switch lower(tp)
        case 'rot'
            k0 = 1; 
            alpha = s * (max(k0, kmax) - k0) .* dt.^(1-eta.^p);
            sigma = 1 ./ (1 + alpha .* k.^4);

        case 'exp'
            sigma = exp(-s * eta.^p);
            sigma(eta <= c) = 1;

        otherwise
            error('Unknown filter type: %s (supported: rot, exp)', filterType);
    end

    Z(:, 1:Nu) = Z(:, 1:Nu) .* reshape(sigma, [1, Nu]);
end

function Z = limiting(Z, basis, grid, params)
% LIMITING Apply slope/positivity limiter to control oscillations.
%
%   Z = limiting(Z, basis, grid, params) applies the Zhang-Shu limiter
%   to the URM decomposition Z = [U; R] to control numerical oscillations
%   and ensure physical constraints (e.g., positive density).

    Nu = params.Nu;
    Nr = params.Nr;
    tp = params.limiterType;
    tol = params.limiterEps;

    switch lower(tp)
        case 'cutoff'
            % Simple cutoff limiter
            if Nu > 0
                Z(:, 1) = max(0, Z(:, 1));
            end

        case 'ls-r'
            if Nu > 0
                Z = applyLinearScalingLimiter(Z, Nu, Nr, grid, params);
            end

        otherwise
            error('Unknown limiter type: %s (supported: cutoff, positive, ls-r, opt-r)', tp);
    end
end

function Z = applyZhangShuLimiter(Z, Nu, Nr, params)
% APPLYZHANGSHULIMITER Zhang-Shu positivity-preserving slope limiter.
%
%   Scales higher-order Fourier modes to ensure positive density while
%   preserving orthogonality of the decomposition.

    eps = params.limiterEps;
    Nx = size(Z, 1);

    if Nu <= 1
        return;  % No limiting needed for constant mode only
    end

    % Extract macro coefficients U
    U = Z(:, 1:Nu);

    % Compute cell average (zeroth moment): rho = sqrt(2π) * U(1)
    rho_avg = sqrt(2*pi) * U(:, 1);

    % Estimate maximum oscillation from higher modes
    oscillation = sum(abs(U(:, 2:Nu)), 2);

    % Minimum value of reconstructed distribution
    f_min = rho_avg - oscillation;

    % Find cells where limiting is needed
    needLimit = f_min < eps;

    if any(needLimit)
        % Compute limiting factor θ ∈ (0, 1]
        theta = ones(Nx, 1);
        idx = needLimit & (oscillation > 1e-14);
        theta(idx) = min(1, (rho_avg(idx) - eps) ./ oscillation(idx));
        theta = max(0, min(1, theta));

        % Apply limiter: scale higher modes by θ
        for k = 2:Nu
            U(:, k) = theta .* U(:, k);
        end

        Z(:, 1:Nu) = U;
    end

    % Final safety check
    minU0 = eps / sqrt(2*pi);
    Z(:, 1) = max(Z(:, 1), minU0);
end

function Z = applyLinearScalingLimiter(Z, Nu, Nr, grid, params)
% APPLYLINEARSCALINGLIMITER Linear scaling realizability limiter.
%
%   Implements the ls-r limiter from Laiu et al. (2019), equation (4.11).
%   For each cell, finds maximum α ∈ [0,1] such that:
%       Z_ls = [u; α*ũ]
%   satisfies the realizability conditions (C1)-(C6) from Theorem 4.1.
%
%   Reference: Laiu et al., "A Positive Asymptotic-Preserving Scheme for
%   Linear Kinetic Transport Equations", SIAM J. Sci. Comput. (2019).

    Nx = size(Z, 1);
    eps = params.limiterEps;
    dx = grid.dx;

    % Extract U = [u; ũ] for each cell
    U = Z(:, 1:Nu);

    % Compute coefficient vectors for realizability conditions
    [ax, ay, ax2, ay2, axy] = computeRealizabilityCoeffs(Nu, grid, params);

    % Apply limiter to each cell
    for i = 1:Nx
        u_i = U(i, :).';  % Column vector [Nu x 1]
        u0 = u_i(1);      % Zeroth moment (mean)

        % Skip if mean is already negative or near zero
        if u0 < eps
            u_i(2:end) = 0;
            U(i, :) = u_i.';
            continue;
        end

        u_tilde = u_i(2:end);  % Micro coefficients

        % Check if limiting is needed
        if checkRealizabilityConditions(u_i, ax, ay, ax2, ay2, axy, dx, eps)
            continue;  % Already realizable
        end

        % Binary search for maximum α ∈ [0,1]
        alpha_min = 0;
        alpha_max = 1;
        alpha = 1;

        % Binary search with tolerance
        max_iter = 50;
        tol_alpha = 1e-12;

        for iter = 1:max_iter
            alpha = (alpha_min + alpha_max) / 2;

            % Test scaled coefficients
            u_test = [u0; alpha * u_tilde];

            if checkRealizabilityConditions(u_test, ax, ay, ax2, ay2, axy, dx, eps)
                % Realizable, try larger alpha
                alpha_min = alpha;
            else
                % Not realizable, try smaller alpha
                alpha_max = alpha;
            end

            if alpha_max - alpha_min < tol_alpha
                break;
            end
        end

        % Apply final scaling
        alpha = alpha_min;
        u_i(2:end) = alpha * u_tilde;
        U(i, :) = u_i.';
    end

    Z(:, 1:Nu) = U;
end

function Z = applyOptimizationLimiter(Z, Nu, Nr, grid, params)
% APPLYOPTIMIZATIONLIMITER Optimization-based realizability limiter.
%
%   Implements the opt-r limiter from Laiu et al. (2019), equation (4.12).
%   For each cell, solves the constrained optimization problem:
%       minimize   (1/2)||ṽ - ũ||₂²
%       subject to [u; ṽ] satisfies (C1)-(C6)
%
%   Reference: Laiu et al., "A Positive Asymptotic-Preserving Scheme for
%   Linear Kinetic Transport Equations", SIAM J. Sci. Comput. (2019).

    Nx = size(Z, 1);
    eps = params.limiterEps;
    dx = grid.dx;

    % Extract U = [u; ũ] for each cell
    U = Z(:, 1:Nu);

    % Compute coefficient vectors for realizability conditions
    [ax, ay, ax2, ay2, axy] = computeRealizabilityCoeffs(Nu, grid, params);

    % Apply limiter to each cell
    for i = 1:Nx
        u_i = U(i, :).';  % Column vector [Nu x 1]
        u0 = u_i(1);      % Zeroth moment (mean)

        % Skip if mean is already negative or near zero
        if u0 < eps
            u_i(2:end) = 0;
            U(i, :) = u_i.';
            continue;
        end

        u_tilde = u_i(2:end);  % Micro coefficients [Ne x 1]
        Ne = length(u_tilde);

        % Check if limiting is needed
        if checkRealizabilityConditions(u_i, ax, ay, ax2, ay2, axy, dx, eps)
            continue;  % Already realizable
        end

        % Set up quadratic programming problem
        % minimize (1/2) v' * H * v - f' * v
        % subject to A * v <= b

        H = eye(Ne);
        f = u_tilde;

        % Build inequality constraints from realizability conditions
        [A_ineq, b_ineq] = buildRealizabilityConstraints(u0, ax, ay, ax2, ay2, axy, dx, eps);

        % Solve QP using quadprog
        options = optimoptions('quadprog', 'Display', 'off');

        try
            v_opt = quadprog(H, -f, A_ineq, b_ineq, [], [], [], [], [], options);

            if ~isempty(v_opt)
                u_i(2:end) = v_opt;
            else
                % Fallback to linear scaling if QP fails
                warning('QP solver failed at cell %d, using ls-r fallback', i);
                u_test = applyLinearScalingToCell(u_i, ax, ay, ax2, ay2, axy, dx, eps);
                u_i = u_test;
            end
        catch
            % Fallback to linear scaling if QP fails
            u_test = applyLinearScalingToCell(u_i, ax, ay, ax2, ay2, axy, dx, eps);
            u_i = u_test;
        end

        U(i, :) = u_i.';
    end

    Z(:, 1:Nu) = U;
end

function [ax, ay, ax2, ay2, axy] = computeRealizabilityCoeffs(Nu, grid, params)
% COMPUTEREALIZABILITYCOEFFS Compute coefficient vectors for realizability.
%
%   Computes the vectors ax, ay, ax2, ay2, axy used in realizability
%   conditions (C1)-(C6) from Theorem 4.1 of Laiu et al. (2019).
%
%   These are inner products of basis functions with angular moments:
%   ax = [⟨m²Ωₓ⟩, ⟨m̃mΩₓ⟩ᵀ]ᵀ
%   ay = [⟨m²Ωᵧ⟩, ⟨m̃mΩᵧ⟩ᵀ]ᵀ
%   ax2 = [1/3, ⟨m̃mΩ²ₓ⟩ᵀ]ᵀ
%   ay2 = [1/3, ⟨m̃mΩ²ᵧ⟩ᵀ]ᵀ
%   axy = [0, ⟨m̃mΩₓΩᵧ⟩ᵀ]ᵀ

    Nq = params.Nv_u;
    angle = grid.angle_u;

    % Velocity components
    vx = cos(angle);
    vy = sin(angle);

    % Quadrature weights
    w = grid.w_u;

    % Fourier basis functions
    E = computeFourierBasis(angle, params);  % [Nq x Nu]

    % Compute coefficient vectors
    ax = zeros(Nu, 1);
    ay = zeros(Nu, 1);
    ax2 = zeros(Nu, 1);
    ay2 = zeros(Nu, 1);
    axy = zeros(Nu, 1);

    for k = 1:Nu
        m_k = E(:, k);  % k-th basis function

        % ax(k) = ⟨m²Ωₓ⟩ = ∫ m_k² * vx dΩ
        ax(k) = sum(w .* (m_k.^2) .* vx);

        % ay(k) = ⟨m²Ωᵧ⟩ = ∫ m_k² * vy dΩ
        ay(k) = sum(w .* (m_k.^2) .* vy);

        % ax2(k) = ⟨m²Ω²ₓ⟩ = ∫ m_k² * vx² dΩ
        ax2(k) = sum(w .* (m_k.^2) .* (vx.^2));

        % ay2(k) = ⟨m²Ω²ᵧ⟩ = ∫ m_k² * vy² dΩ
        ay2(k) = sum(w .* (m_k.^2) .* (vy.^2));

        % axy(k) = ⟨m²ΩₓΩᵧ⟩ = ∫ m_k² * vx*vy dΩ
        axy(k) = sum(w .* (m_k.^2) .* vx .* vy);
    end
end

function ok = checkRealizabilityConditions(u, ax, ay, ax2, ay2, axy, dx, eps)
% CHECKREALIZABILITYCONDITIONS Check if conditions (C1)-(C6) are satisfied.
%
%   Returns true if the coefficient vector u = [u0; ũ] satisfies all six
%   realizability conditions from Theorem 4.1 of Laiu et al. (2019).

    u0 = u(1);

    % Evaluate linear combinations
    ax_u = ax.' * u;
    ay_u = ay.' * u;
    ax2_u = ax2.' * u;
    ay2_u = ay2.' * u;
    axy_u = axy.' * u;

    % (C1) u0 ± ax'*u >= 0
    if u0 + ax_u < -eps || u0 - ax_u < -eps
        ok = false;
        return;
    end

    % (C2) u0 ± ay'*u >= 0
    if u0 + ay_u < -eps || u0 - ay_u < -eps
        ok = false;
        return;
    end

    % (C3) u0 >= ax2'*u >= 0
    if ax2_u < -eps || u0 - ax2_u < -eps
        ok = false;
        return;
    end

    % (C4) u0 >= ay2'*u >= 0
    if ay2_u < -eps || u0 - ay2_u < -eps
        ok = false;
        return;
    end

    % (C5) u0 ± 2*axy'*u >= 0
    if u0 + 2*axy_u < -eps || u0 - 2*axy_u < -eps
        ok = false;
        return;
    end

    % (C6) (ax2'/Δx² ± 2axy'/(ΔxΔy) + ay2'/Δy²)*u >= 0
    dx1 = dx(1);
    dx2 = dx(2);

    coeff1 = ax2_u / (dx1^2) + 2*axy_u / (dx1*dx2) + ay2_u / (dx2^2);
    coeff2 = ax2_u / (dx1^2) - 2*axy_u / (dx1*dx2) + ay2_u / (dx2^2);

    if coeff1 < -eps || coeff2 < -eps
        ok = false;
        return;
    end

    ok = true;
end

function [A, b] = buildRealizabilityConstraints(u0, ax, ay, ax2, ay2, axy, dx, eps)
% BUILDREALIZABILITYCONSTRAINTS Build inequality constraints for QP.
%
%   Constructs linear inequality constraints A*v <= b encoding the
%   realizability conditions (C1)-(C6) for the micro coefficients v.

    Ne = length(ax) - 1;  % Number of micro coefficients

    % Extract micro parts of coefficient vectors
    ax_tilde = ax(2:end);
    ay_tilde = ay(2:end);
    ax2_tilde = ax2(2:end);
    ay2_tilde = ay2(2:end);
    axy_tilde = axy(2:end);

    % Scalar parts
    ax0 = ax(1);
    ay0 = ay(1);
    ax20 = ax2(1);
    ay20 = ay2(1);
    axy0 = axy(1);

    % Initialize constraint matrices
    A = [];
    b = [];

    % (C1) u0 ± (ax0 + ax_tilde'*v) >= 0
    % => ±ax_tilde'*v <= u0 ∓ ax0 - eps
    A = [A; ax_tilde.'];
    b = [b; u0 - ax0 - eps];
    A = [A; -ax_tilde.'];
    b = [b; u0 + ax0 - eps];

    % (C2) u0 ± (ay0 + ay_tilde'*v) >= 0
    % => ±ay_tilde'*v <= u0 ∓ ay0 - eps
    A = [A; ay_tilde.'];
    b = [b; u0 - ay0 - eps];
    A = [A; -ay_tilde.'];
    b = [b; u0 + ay0 - eps];

    % (C3) u0 >= ax20 + ax2_tilde'*v >= 0
    % => ax2_tilde'*v <= u0 - ax20 - eps
    % => -ax2_tilde'*v <= -ax20 - eps
    A = [A; ax2_tilde.'];
    b = [b; u0 - ax20 - eps];
    A = [A; -ax2_tilde.'];
    b = [b; -ax20 - eps];

    % (C4) u0 >= ay20 + ay2_tilde'*v >= 0
    % => ay2_tilde'*v <= u0 - ay20 - eps
    % => -ay2_tilde'*v <= -ay20 - eps
    A = [A; ay2_tilde.'];
    b = [b; u0 - ay20 - eps];
    A = [A; -ay2_tilde.'];
    b = [b; -ay20 - eps];

    % (C5) u0 ± 2*(axy0 + axy_tilde'*v) >= 0
    % => ±2*axy_tilde'*v <= u0 ∓ 2*axy0 - eps
    A = [A; 2*axy_tilde.'];
    b = [b; u0 - 2*axy0 - eps];
    A = [A; -2*axy_tilde.'];
    b = [b; u0 + 2*axy0 - eps];

    % (C6) (ax2/Δx² ± 2axy/(ΔxΔy) + ay2/Δy²) >= 0
    dx1 = dx(1);
    dx2 = dx(2);

    % Positive combination
    coeff = ax2_tilde / (dx1^2) + 2*axy_tilde / (dx1*dx2) + ay2_tilde / (dx2^2);
    const = ax20 / (dx1^2) + 2*axy0 / (dx1*dx2) + ay20 / (dx2^2);
    A = [A; -coeff.'];
    b = [b; -const - eps];

    % Negative combination
    coeff = ax2_tilde / (dx1^2) - 2*axy_tilde / (dx1*dx2) + ay2_tilde / (dx2^2);
    const = ax20 / (dx1^2) - 2*axy0 / (dx1*dx2) + ay20 / (dx2^2);
    A = [A; -coeff.'];
    b = [b; -const - eps];
end


function u_out = applyLinearScalingToCell(u, ax, ay, ax2, ay2, axy, dx, eps)
% APPLYLINEARSCALINGTOCELL Apply linear scaling limiter to a single cell.
%
%   Helper function used as fallback when optimization fails.

    u0 = u(1);
    u_tilde = u(2:end);

    % Binary search for maximum α ∈ [0,1]
    alpha_min = 0;
    alpha_max = 1;

    max_iter = 50;
    tol_alpha = 1e-12;

    for iter = 1:max_iter
        alpha = (alpha_min + alpha_max) / 2;

        % Test scaled coefficients
        u_test = [u0; alpha * u_tilde];

        if checkRealizabilityConditions(u_test, ax, ay, ax2, ay2, axy, dx, eps)
            alpha_min = alpha;
        else
            alpha_max = alpha;
        end

        if alpha_max - alpha_min < tol_alpha
            break;
        end
    end

    % Apply final scaling
    u_out = [u0; alpha_min * u_tilde];
end


function D = computeArtificialDissipation(U, dx, params)
% COMPUTEARTIFICIALDISSIPATION Compute artificial dissipation term.
%
%   D = computeArtificialDissipation(U, dx, params) computes the
%   artificial dissipation term from Laiu et al. (2019), equation (3.10):
%   -C_AD * sum_d (dx(d)^3 * delta_d^4) U for arbitrary dimensions.
%
%   For the case with sigma=0 and epsilon=1, C_AD = 1/(2-theta).

    arguments
        U double
        dx double {mustBeVector}
        params struct
    end

    theta = params.thetaMinmod;
    C_AD = 1 / (2 - theta);

    Nd = length(dx);
    D = zeros(size(U));

    for d = 1:Nd
        delta4_d = compute4thOrderDissipation(U, dx(d), theta, d);
        D = D - C_AD * dx(d)^3 * delta4_d;
    end
end


function delta4 = compute4thOrderDissipation(w, h, theta, dim)
% COMPUTE4THORDERDISSIPATION Compute 4th-order dissipation operator.
%
%   delta4 = compute4thOrderDissipation(w, h, theta, dim) computes the
%   4th-order dissipation operator delta^4 from equations (3.12)-(3.14)
%   in Laiu et al. (2019) for arbitrary dimensions.

    arguments
        w double
        h double {mustBePositive}
        theta double {mustBePositive}
        dim double {mustBePositive, mustBeInteger}
    end

    sz = size(w);
    N = sz(dim);

    w_plus = zeros(sz);
    w_minus = zeros(sz);

    idx_prev = repmat({':'}, 1, length(sz));
    idx_curr = repmat({':'}, 1, length(sz));
    idx_next = repmat({':'}, 1, length(sz));

    idx_prev{dim} = 1:N-2;
    idx_curr{dim} = 2:N-1;
    idx_next{dim} = 3:N;

    w_prev = w(idx_prev{:});
    w_curr = w(idx_curr{:});
    w_next = w(idx_next{:});

    slope_left = (w_curr - w_prev) / h;
    slope_right = (w_next - w_curr) / h;

    slope_limited = minmod(slope_right, slope_left, theta);

    w_plus(idx_curr{:}) = w_curr - 0.5 * h * slope_limited;
    w_minus(idx_curr{:}) = w_curr + 0.5 * h * slope_limited;

    idx_first = repmat({':'}, 1, length(sz));
    idx_first{dim} = 1;
    w_plus(idx_first{:}) = w(idx_first{:});
    w_minus(idx_first{:}) = w(idx_first{:});

    idx_last = repmat({':'}, 1, length(sz));
    idx_last{dim} = N;
    w_plus(idx_last{:}) = w(idx_last{:});
    w_minus(idx_last{:}) = w(idx_last{:});

    delta4 = zeros(sz);

    jump_right = w_plus(idx_next{:}) - w_minus(idx_curr{:});
    jump_left = w_plus(idx_curr{:}) - w_minus(idx_prev{:});

    delta4(idx_curr{:}) = (jump_right - jump_left) / h^4;
end


function s = minmod(a, b, theta)
% MINMODLIMITER Minmod slope limiter.
%
%   s = minmod(a, b, theta) computes the minmod-limited slope
%   from equation (3.14) in Laiu et al. (2019).

    arguments
        a double
        b double
        theta double {mustBePositive} = 1.5
    end

    s = zeros(size(a));

    c = (a + b) / 2;

    pos_mask = (a > 0) & (b > 0);
    neg_mask = (a < 0) & (b < 0);

    s(pos_mask) = min(theta * a(pos_mask),...
                      min(theta * b(pos_mask), c(pos_mask)));

    s(neg_mask) = max(theta * a(neg_mask),...
                      max(theta * b(neg_mask), c(neg_mask)));
end


function I = trapzNd(F, dx)
% TRAPZND Multidimensional trapezoidal integration.
%
%   I = trapzNd(F, dx) computes the integral of @a F over a rectangular
%   domain using the trapezoidal rule. @a dx is a vector containing grid
%   spacings for each dimension.

    Nd = length(dx);
    I = F;
    for d = 1:Nd
        I = trapz(I, d) * dx(d);
    end
end

function F_out = rbfInterpolate(R, angle_in, angle_out, params)
% RBFINTERPOLATE Interpolate residual from input grid to output grid using RBF.
%
%   F_out = rbfInterpolate(R, angle_in, angle_out, params) interpolates
%   residual values @a R from input angles @a angle_in to output angles
%   @a angle_out using radial basis function interpolation on S^1.
%
%   Input:
%   - R: [Nx x Nr] residual values at input angles
%   - angle_in: [Nr x 1] input angles
%   - angle_out: [Nv x 1] output angles
%   - params: parameter structure with rbfType and rbfEpsilon
%
%   Output:
%   - F_out: [Nx x Nv] interpolated values at output angles

    Nx = size(R, 1);
    Nr = length(angle_in);
    Nv = length(angle_out);

    % Build RBF kernel matrix A: A(i,j) = phi(||angle_in(i) - angle_in(j)||)
    A = zeros(Nr, Nr);
    for i = 1:Nr
        for j = 1:Nr
            dtheta = abs(angle_in(i) - angle_in(j));
            dtheta = min(dtheta, 2*pi - dtheta);
            A(i, j) = evaluateRBF(dtheta, params);
        end
    end

    % Build interpolation matrix B: B(i,j) = phi(||angle_out(i) - angle_in(j)||)
    B = zeros(Nv, Nr);
    for i = 1:Nv
        for j = 1:Nr
            dtheta = abs(angle_out(i) - angle_in(j));
            dtheta = min(dtheta, 2*pi - dtheta);
            B(i, j) = evaluateRBF(dtheta, params);
        end
    end

    % Check conditioning
    condA = cond(A);
    if condA > 1e12
        warning('RBF interpolation matrix is ill-conditioned (cond=%.2e). Consider changing rbfType or rbfEpsilon.', condA);
    end

    % Solve RBF interpolation system: A * W = R' => W = A \ R'
    % Then interpolate: F_out = B * W = B * (A \ R')
    % For multiple spatial points, solve: A * W = R' where R is [Nx x Nr]
    % Result: W is [Nr x Nx], F_out = (B * W)' = W' * B' is [Nx x Nv]

    % Add small regularization for stability
    reg = 1e-10 * trace(A) / Nr;
    A = A + reg * eye(Nr);

    % Solve: A * W = R' (W is [Nr x Nx])
    W = A \ R';

    % Interpolate: F_out = W' * B' = (B * W)'
    F_out = (B * W)';
end

function phi = evaluateRBF(r, params)
% EVALUATERBF Evaluate radial basis function at distance r.
%
%   phi = evaluateRBF(r, params) computes RBF value at distance @a r
%   using the RBF type and shape parameter specified in @a params.
%
%   Supported RBF types:
%   - 'wendland_c2': Wendland C^2 compactly supported RBF
%   - 'gaussian': Gaussian RBF exp(-epsilon^2 * r^2)
%   - 'multiquadric': Multiquadric RBF sqrt(1 + (epsilon*r)^2)

    epsilon = params.rbfEpsilon;

    switch lower(params.rbfType)
        case 'wendland_c2'
            % Wendland C^2: (1 - εr)₊⁴ (4εr + 1)
            s = epsilon * r;
            if s >= 1
                phi = 0;
            else
                phi = (1 - s)^4 * (4*s + 1);
            end

        case 'gaussian'
            % Gaussian RBF
            phi = exp(-(epsilon * r)^2);

        case 'multiquadric'
            % Multiquadric RBF
            phi = sqrt(1 + (epsilon * r)^2);

        otherwise
            error('Unknown RBF type: %s', params.rbfType);
    end
end
