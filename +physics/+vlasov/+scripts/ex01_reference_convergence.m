%==========================================================================
% FileName: ex01_reference_convergence.m
% Description:
%   Convergence against a single fine reference solution for the ex01ldh
%   Landau damping test.
%
%   1. Compute (or load) a reference solution with nx = refNx and
%      refModes Hermite modes.
%   2. Compute coarser solutions and measure
%          || f_(h,N) - f_ref ||_{L2}  (all Hermite modes, coarse
%                                       zero-padded),
%      split into modes k < kSplit and k >= kSplit, plus the errors on
%      Q (electric field) and Psi (potential).
%
%   Two studies are run for each tau0:
%     'space' : nModes fixed (even), h refined         -> spatial order
%     'joint' : nModes = nx (same label as the MATLAB runs) -> joint
%
%   Mesh convention (as everywhere in this module):
%       UniformGrid(xBBox, nx/2), i.e. nCells = nx/2.
%   With the default refNx = 512, refModes = 512 the reference has 256
%   cells and 512 modes, i.e. the "512 x 512" level of ex01_convergence.
%
%   Comparison across meshes: the coarse solution is evaluated on the
%   quadrature points of the reference mesh by splitting each coarse
%   reference cell [-1/2, 1/2] into r = nCellsRef/nCells sub-intervals
%   (nested uniform meshes, any power-of-two ratio). The physical
%   coordinates of both point sets are checked to coincide.
%
%   Reference quality: the reference itself carries a spatial error.
%   With refCheck = true it is estimated as
%       |f_ref - f_(refNx/2, refModes)| / (2^p - 1),  p = 2 (P1),
%   and rows whose error is within 10x this estimate are flagged.
%   Rule of thumb: keep the finest tested mesh at least 4x coarser than
%   the reference.
%==========================================================================

clc;
clear;
close all;

%% USER PARAMETERS
tau0List   = [10, 1e3, 1e5];           % collision times
refNx      = 512;                      % reference nx label (cells = nx/2)
refModes   = 512;                      % reference number of Hermite modes
refCheck   = true;                     % estimate the reference error
expectedP  = 2;                        % expected spatial order (P1)

study(1).name   = 'space';
study(1).nx     = [16, 32, 64, 128];
study(1).nModes = [16, 16, 16, 16];    % fixed, even

study(2).name   = 'joint';
study(2).nx     = [16, 32, 64, 128];
study(2).nModes = [16, 32, 64, 128];   % same labels as ex01_convergence

kSplit     = 16;                       % split index for reporting
tFinal     = 0.5;                        % final time
dt         = 1e-3;                     % time step (same for all runs)
delta      = 0.01;                     % perturbation amplitude
density    = 4;                        % quadrature density for errors
reuseRef   = true;                     % load reference from disk if present

outDir = fullfile(fileparts(mfilename('fullpath')), 'ex01ref');
if ~exist(outDir, 'dir'), mkdir(outDir); end

for s = 1:numel(study)
    assert(all(mod(refNx ./ study(s).nx, 1) == 0) && ...
        all(log2(refNx ./ study(s).nx) == round(log2(refNx ./ study(s).nx))), ...
        'refNx / nx must be a power of two (nested meshes).');
    assert(all(study(s).nModes <= refModes), ...
        'Tested nModes must not exceed refModes.');
end

%% COMMON CONFIGURATION (mirrors ex01_convergence.m, id = 'ex01ldh')
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config0 = physics.scheme.ConfigParser.parseFile(common);

config0.id = 'ex01ldh';
config0.nDims = 1;
config0.L = 4*pi;
config0.xBBox = [0, config0.L];
config0.kappa = 2*pi/config0.L*pi/3;
config0.T0 = 1;
config0.delta = delta;
config0.epsilon = 1;
config0.tFinal = tFinal;
config0.vBBox = [-3.2, 3.2];
config0.ckptTimeStamps = [];
config0.bc = [];
config0.exact = [];
config0.verbose = 0;
config0.cfl = [];
config0.dt = dt;

config0.phiInf = phiInf(config0);
config0.EInf = EInf(config0);
config0.cInf = cInf(config0);
config0.rhoInf = rhoInf(config0);
config0.sqrtRhoInf = sqrtRhoInf(config0);
config0.dSqrtRhoInf = dSqrtRhoInf(config0);
config0.ic = DInit(config0);

fprintf('[M] xBasisOrder = %d (%s), T = %g, dt = %g\n', ...
    config0.xBasisOrder, config0.xBasisType, tFinal, dt);
fprintf('[M] Reference: nx = %d (%d cells), %d modes\n', ...
    refNx, refNx/2, refModes);

%% MAIN LOOP
for tau0 = tau0List
    cfgT = config0;
    cfgT.tau0 = tau0;
    cfgT.tau = cfgT.tau0 * cfgT.epsilon;
    fprintf('\n================ tau0 = %g ================\n', tau0);

    % ---- Reference -------------------------------------------------------
    refFile = fullfile(outDir, sprintf( ...
        'ref_tau0%.0e_nx%d_modes%d_P%d_T%g_dt%.0e.mat', ...
        tau0, refNx, refModes, config0.xBasisOrder-1, tFinal, dt));
    if reuseRef && isfile(refFile)
        fprintf('[M] Loading reference from %s\n', refFile);
        ref = physics.vlasov.ApHermiteState.load(refFile);
    else
        ref = runOne(cfgT, refNx, refModes);
        ref.save(refFile);
    end

    refErr = NaN;
    if refCheck
        refHalf = runOne(cfgT, refNx/2, refModes);
        e = crossMeshModeErrors(ref, refHalf, 'D', density);
        refErr = norm(e) / (2^expectedP - 1);
        fprintf('[M] Estimated reference error ~ %.3e\n', refErr);
    end

    % ---- Studies ---------------------------------------------------------
    for s = 1:numel(study)
        nxs = study(s).nx;
        nms = study(s).nModes;
        n = numel(nxs);
        eAll = zeros(n, 1); eLow = zeros(n, 1); eHigh = zeros(n, 1);
        eQ = zeros(n, 1); ePsi = zeros(n, 1);
        for i = 1:n
            st = runOne(cfgT, nxs(i), nms(i));
            e = crossMeshModeErrors(ref, st, 'D', density);
            eAll(i) = norm(e);
            eLow(i) = norm(e(1:min(kSplit, end)));
            eHigh(i) = norm(e(kSplit+1:end));
            eQ(i) = norm(crossMeshModeErrors(ref, st, 'Q', density));
            ePsi(i) = norm(crossMeshModeErrors(ref, st, 'Psi', density));
        end

        nCells = nxs(:)/2;
        nearRef = eAll < 10*refErr;
        T = table(nxs(:), nCells, nms(:), eAll, orders(eAll), ...
            eLow, orders(eLow), eHigh, orders(eHigh), ...
            eQ, orders(eQ), ePsi, orders(ePsi), nearRef, ...
            'VariableNames', {'nx', 'nCells', 'nModes', 'D_all', 'O_all', ...
            sprintf('D_k_lt_%d', kSplit), 'O_low', ...
            sprintf('D_k_ge_%d', kSplit), 'O_high', ...
            'Q', 'O_Q', 'Psi', 'O_Psi', 'NearRefError'});
        fprintf('\n-- study ''%s'' (tau0 = %g) --\n', study(s).name, tau0);
        disp(T);
        if refCheck && any(nearRef)
            fprintf(['    Rows with NearRefError = true are within 10x the ', ...
                'estimated reference error;\n    their order is not ', ...
                'reliable. Use a finer reference or drop those rows.\n']);
        end
        writetable(T, fullfile(outDir, sprintf('%s_tau0%.0e_ref%dx%d.csv', ...
            study(s).name, tau0, refNx, refModes)));
    end
end

fprintf('\n[M] Results written to %s\n', outDir);

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function state = runOne(cfg, nx, nModes)
    cfg.nx = nx;
    cfg.nh = nModes;   % hermite(1, nh) builds nh modes
    fprintf('[M] Run nx = %4d (%4d cells), nModes = %4d ... ', ...
        nx, nx/2, nModes);
    tic;
    [scheme, state] = buildProblem(cfg);
    state = scheme.run(state);
    assert(abs(scheme.TDisc.Timeline.Now - cfg.tFinal) <= ...
        100*eps(max(1, cfg.tFinal)), 'Run did not reach the final time.');
    assert(size(state.Dofs.D, 2) == nModes, ...
        'Expected %d Hermite modes, got %d.', nModes, size(state.Dofs.D, 2));
    fprintf('%.1f s\n', toc);
end

function e = crossMeshModeErrors(ref, st, fieldName, density)
    % Per-component L2 error between a reference state and a state on a
    % nested coarser (or equal) uniform mesh. The coarse solution is
    % evaluated at the reference quadrature points; components missing in
    % either state are zero-padded.
    fineSpace = ref.XDisc;
    coarseSpace = st.XDisc;
    nf = fineSpace.NMeshElements;
    nc = coarseSpace.NMeshElements;
    r = nf / nc;
    assert(r >= 1 && r == round(r), 'Meshes are not nested.');

    metric = physics.analysis.FiniteElementAbsoluteMetric( ...
        fineSpace, reductions={'L2'});
    metric.setPoints(density);
    xRef = metric.Nodes;                       % 1 x nq on [-1/2, 1/2]

    % Sub-cell j of a coarse cell: [-1/2 + j/r, -1/2 + (j+1)/r]
    xSub = cell(1, r);
    for j = 0:r-1
        xSub{j+1} = -1/2 + (j + xRef + 1/2) / r;
    end
    xSub = [xSub{:}];

    % Check that the two point sets coincide physically and in order.
    Xf = fineSpace.Mesh.collocate(xRef);
    Xc = coarseSpace.Mesh.collocate(xSub);
    assert(numel(Xf) == numel(Xc) && ...
        max(abs(Xf(:) - Xc(:))) <= 1e-10 * max(1, max(abs(Xf(:)))), ...
        'Quadrature points of the nested meshes do not match.');

    U = fineSpace.eval(xRef, ref.Dofs.(fieldName));
    V = coarseSpace.eval(xSub, st.Dofs.(fieldName));
    n = max(size(U, 2), size(V, 2));
    U(:, end+1:n) = 0;
    V(:, end+1:n) = 0;

    e = metric.eval(U, V);
    e = e(1, :);
end

function o = orders(e)
    o = nan(size(e));
    for i = 2:numel(e)
        if e(i-1) > 0 && e(i) > 0
            o(i) = log2(e(i-1)/e(i));
        end
    end
end

function [scheme, state] = buildProblem(config)
    scheme = physics.vlasov.ApDghsScheme(config=config);

    vElement = approx.element.L2OrthotopeElement.hermite( ...
        1, config.nh, T=config.T0);
    vDisc = approx.space.SpectralSpace(vElement);

    switch lower(char(scheme.getConfig('xBasisType')))
        case 'modal'
            xElement = approx.element.BH1OrthotopeElement.modal( ...
                config.nDims, scheme.getConfig('xBasisOrder'), ...
                pattern=scheme.getConfig('xBasisPattern'));
        case 'nodal'
            xElement = approx.element.BH1OrthotopeElement.nodal( ...
                config.nDims, scheme.getConfig('xBasisOrder'));
        otherwise
            error('Unsupported xBasisType: %s', ...
                char(scheme.getConfig('xBasisType')));
    end

    xMesh = approx.mesh.UniformGrid(config.xBBox, config.nx/2);
    xDisc = approx.space.FiniteElementSpace(xMesh, xElement);
    state = physics.vlasov.ApHermiteState(xDisc, vDisc);
end

%% ---- Problem data (identical to ex01_time_convergence.m) --------------

function h = DInit(config)
    h = @(x, varargin) DInitImpl(config, x, varargin{:});
end

function D = DInitImpl(config, x, varargin)
    if isempty(varargin)
        nh = config.nh;
    else
        nh = varargin{1};
    end
    rhoInf = rhoInfImpl(config, x);
    D0 = (rhoInf + config.delta*cos(config.kappa*x)) ./ sqrt(rhoInf);
    D = [D0, zeros(1, size(x, 2)*(nh-1))];
end

function h = phiInf(config)
    h = @(x) phiInfImpl(config, x);
end

function phi = phiInfImpl(config, x)
    switch config.id
        case 'ex01ldh'
            phi = zeros(1, size(x, 2));
        case 'ex01ldi'
            phi = 0.2*sin(config.kappa*x);
        otherwise
            error('Unknown test id: %s', config.id);
    end
end

function h = rhoInf(config)
    h = @(x) rhoInfImpl(config, x);
end

function rho = rhoInfImpl(config, x)
    rho = config.cInf * exp(-phiInfImpl(config, x)/config.T0);
end

function h = sqrtRhoInf(config)
    h = @(x) sqrt(rhoInfImpl(config, x));
end

function h = EInf(config)
    syms xs;
    ESym = diff(-phiInfImpl(config, xs), xs);
    if ~has(ESym, xs)
        c = double(ESym);
        h = @(x) c*ones(1, size(x, 2));
    else
        h = matlabFunction(ESym, 'Vars', xs);
    end
end

function h = dSqrtRhoInf(config)
    syms xs;
    sqrtRhoSym = sqrt(config.cInf) * ...
        exp(-phiInfImpl(config, xs)/(2*config.T0));
    dSym = diff(sqrtRhoSym, xs);
    if ~has(dSym, xs)
        c = double(dSym);
        h = @(x) c*ones(1, size(x, 2));
    else
        h = matlabFunction(dSym, 'Vars', xs);
    end
end

function c = cInf(config)
    c = 2*config.L / integral( ...
        @(x) exp(-phiInfImpl(config, x)/config.T0), -config.L, config.L);
end
