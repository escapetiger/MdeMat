%==========================================================================
% FileName: ex01_order_diagnostics.m
% Description:
%   Diagnostics for the spatial order of the AP-DG-Hermite scheme on the
%   ex01ldh Landau damping test (same setting as ex01_convergence.m).
%
%   TEST 1 (parity)  : fixed number of Hermite modes, refine h only.
%                      Compare an even and an odd number of modes.
%                      Expected (P1, alternating fluxes): even -> 2,
%                      odd -> 1.
%   TEST 2 (joint)   : refine h and double the number of modes (as in
%                      ApHermiteState.refine). Report the full phase-space
%                      error and its split into modes k < kSplit and
%                      k >= kSplit.
%   TEST 3 (dt)      : on the finest levels of TEST 2, halve dt and
%                      measure how much the solution moves, compared with
%                      the spatial error. Also redo TEST 2 with dt ~ h.
%   TEST 4 (tail)    : L2 norm of each Hermite mode ||D_k|| at the final
%                      time for several meshes (filamentation tail).
%
%   All errors use consecutive levels (fine vs. coarse, factor 2 in h)
%   with the FiniteElementRichardsonMetric; missing coarse modes are
%   zero-padded, exactly as in Analyzer.addRelative.
%
%   Number of modes convention: hermite(1, nModes) builds nModes modes
%   D_0..D_{nModes-1}, i.e. N_H = nModes - 1.
%==========================================================================

clc;
clear;
close all;

%% USER PARAMETERS
runTest = struct('parity', true, 'joint', true, 'dt', true, 'tail', true);

nxList      = [16, 32, 64, 128, 256];   % spatial levels (factor 2)
parityModes = [16, 17];                 % fixed mode counts for TEST 1
nModes0     = 16;                       % coarsest mode count for TEST 2
kSplit      = 16;                       % split index for reporting
tau0        = 1e5;                      % collision time
delta       = 0.05;                     % perturbation amplitude
tFinal      = 1;                        % final time
dt0         = 1e-3;                     % base time step (as in config.txt)
nxDt        = nxList(end-1:end);        % levels used in TEST 3
nxTail      = [32, 64, 128];            % meshes used in TEST 4
nModesTail  = 64;                       % mode count used in TEST 4
density     = 4;                        % quadrature density for errors

outDir = fullfile(fileparts(mfilename('fullpath')), 'ex01diag');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% COMMON CONFIGURATION (mirrors ex01_convergence.m, id = 'ex01ldh')
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config0 = physics.scheme.ConfigParser.parseFile(common);

config0.id = 'ex01ldh';
config0.nDims = 1;
config0.L = 4*pi;
config0.xBBox = [0, config0.L];
config0.kappa = 2*pi/config0.L;
config0.T0 = 1;
config0.delta = delta;
config0.epsilon = 1;
config0.tau0 = tau0;
config0.tau = config0.tau0 * config0.epsilon;
config0.tFinal = tFinal;
config0.vBBox = [-3.2, 3.2];
config0.ckptTimeStamps = [];
config0.bc = [];
config0.exact = [];
config0.verbose = 0;
config0.cfl = [];
config0.dt = dt0;

config0.phiInf = phiInf(config0);
config0.EInf = EInf(config0);
config0.cInf = cInf(config0);
config0.rhoInf = rhoInf(config0);
config0.sqrtRhoInf = sqrtRhoInf(config0);
config0.dSqrtRhoInf = dSqrtRhoInf(config0);
config0.ic = DInit(config0);

fprintf('[M] xBasisOrder = %d (%s), tSplit = %s, dt0 = %g, tau0 = %g\n', ...
    config0.xBasisOrder, config0.xBasisType, config0.tSplitName, dt0, tau0);

% Cache of final states, keyed by (nx, nModes, dt), so that runs shared
% between tests are only computed once.
cache = containers.Map('KeyType', 'char', 'ValueType', 'any');

%% TEST 1: PARITY OF THE NUMBER OF MODES (h-refinement only)
if runTest.parity
    fprintf('\n===== TEST 1: fixed number of modes, refine h only =====\n');
    for nModes = parityModes
        nxs = nxList;
        nms = repmat(nModes, size(nxs));
        dts = repmat(dt0, size(nxs));
        [states, cache] = runLevels(config0, nxs, nms, dts, cache);
        T = levelTable(states, nxs, nms, dts, kSplit, density);
        fprintf('\n-- nModes = %d (N_H = %d, %s) --\n', nModes, nModes-1, ...
            ternary(mod(nModes, 2) == 0, 'even #modes', 'odd #modes'));
        disp(T);
        writetable(T, fullfile(outDir, sprintf('test1_parity_nModes%d.csv', nModes)));
    end
end

%% TEST 2: JOINT REFINEMENT (h and number of modes)
if runTest.joint
    fprintf('\n===== TEST 2: joint refinement (h, nModes doubled) =====\n');
    nxs = nxList;
    nms = nModes0 * 2.^(0:numel(nxs)-1);
    dts = repmat(dt0, size(nxs));
    [states, cache] = runLevels(config0, nxs, nms, dts, cache);
    T = levelTable(states, nxs, nms, dts, kSplit, density);
    disp(T);
    writetable(T, fullfile(outDir, 'test2_joint_dtFixed.csv'));
end

%% TEST 3: TIME-STEP SENSITIVITY ON THE FINEST LEVELS
if runTest.dt
    fprintf('\n===== TEST 3a: dt sensitivity on finest joint levels =====\n');
    rows = {};
    for nx = nxDt
        nModes = nModes0 * nx / nxList(1);
        dtsLoc = dt0 ./ [1, 2, 4];
        [st, cache] = runLevels(config0, repmat(nx, 1, 3), ...
            repmat(nModes, 1, 3), dtsLoc, cache);
        % Same mesh, same modes: plain per-mode L2 differences.
        e12 = sameSpaceModeErrors(st{1}, st{2}, 'D', density);
        e24 = sameSpaceModeErrors(st{2}, st{3}, 'D', density);
        rows(end+1, :) = {nx, nModes, dtsLoc(1), norm(e12), norm(e24), ...
            observedOrder(norm(e12), norm(e24)), ...
            norm(e12(1:min(kSplit, end))), norm(e12(kSplit+1:end))}; %#ok<AGROW>
    end
    T3 = cell2table(rows, 'VariableNames', {'nx', 'nModes', 'dt', ...
        'diff_dt_dt2', 'diff_dt2_dt4', 'TimeOrder', ...
        'diff_dt_dt2_kLow', 'diff_dt_dt2_kHigh'});
    disp(T3);
    fprintf(['    Compare diff_dt_dt2 with the spatial error of TEST 2 at the ', ...
        'same level.\n    If comparable, dt0 pollutes the spatial order.\n']);
    writetable(T3, fullfile(outDir, 'test3a_dt_sensitivity.csv'));

    fprintf('\n===== TEST 3b: joint refinement with dt ~ h =====\n');
    nxs = nxList;
    nms = nModes0 * 2.^(0:numel(nxs)-1);
    dts = dt0 * nxList(1) ./ nxs;   % dt halves with h
    [states, cache] = runLevels(config0, nxs, nms, dts, cache);
    T = levelTable(states, nxs, nms, dts, kSplit, density);
    disp(T);
    writetable(T, fullfile(outDir, 'test3b_joint_dtScaled.csv'));
end

%% TEST 4: HERMITE TAIL ||D_k|| AT FINAL TIME
if runTest.tail
    fprintf('\n===== TEST 4: Hermite spectrum of the numerical solution =====\n');
    nxs = nxTail;
    nms = repmat(nModesTail, size(nxs));
    dts = repmat(dt0, size(nxs));
    [states, cache] = runLevels(config0, nxs, nms, dts, cache);

    k = (0:nModesTail-1).';
    S = table(k);
    exact = delta * (config0.kappa*tFinal).^k ./ sqrt(factorial(k)); % free-streaming size
    S.freeStreamingEstimate = exact;
    figure('Name', 'Hermite tail'); hold on;
    for i = 1:numel(nxs)
        nk = modeNorms(states{i}, 'D', density);
        nk(1) = NaN;  % D_0 contains the equilibrium; not part of the tail
        S.(sprintf('nx%d', nxs(i))) = nk(:);
        semilogy(k, nk, 'o-', 'DisplayName', sprintf('nx = %d', nxs(i)));
    end
    semilogy(k(2:end), exact(2:end), 'k--', 'DisplayName', '\delta(\kappa t)^k/\surd k!');
    set(gca, 'YScale', 'log');
    xlabel('k'); ylabel('||D_k||_{L^2}'); legend('Location', 'best'); grid on;
    title(sprintf('nModes = %d, t = %g', nModesTail, tFinal));
    hold off;
    disp(S(1:8:end, :));
    fprintf(['    If ||D_k|| is far above the free-streaming estimate and ', ...
        'scales like h (not h^2),\n    the tail is spatial error carried ', ...
        'into the high modes.\n']);
    writetable(S, fullfile(outDir, 'test4_tail.csv'));
    savefig(fullfile(outDir, 'test4_tail.fig'));
end

fprintf('\n[M] Results written to %s\n', outDir);

%% ========================================================================
%% LOCAL FUNCTIONS
%% ========================================================================

function [states, cache] = runLevels(config0, nxs, nms, dts, cache)
    % Run (or fetch from cache) the final state for each (nx, nModes, dt).
    states = cell(numel(nxs), 1);
    for i = 1:numel(nxs)
        key = sprintf('nx%d_nm%d_dt%.6e', nxs(i), nms(i), dts(i));
        if isKey(cache, key)
            states{i} = cache(key);
            continue;
        end
        cfg = config0;
        cfg.nx = nxs(i);
        cfg.nh = nms(i);       % hermite(1, nh) builds nh modes
        cfg.dt = dts(i);
        cfg.cfl = [];
        fprintf('[M] Run nx = %4d, nModes = %4d, dt = %.3e ... ', ...
            nxs(i), nms(i), dts(i));
        tic;
        [scheme, state] = buildProblem(cfg);
        state = scheme.run(state);
        tReached = scheme.TDisc.Timeline.Now;
        assert(abs(tReached - cfg.tFinal) <= 100*eps(max(1, cfg.tFinal)), ...
            'Run did not reach the final time.');
        assert(size(state.Dofs.D, 2) == nms(i), ...
            'Expected %d Hermite modes, got %d.', nms(i), size(state.Dofs.D, 2));
        fprintf('%.1f s\n', toc);
        states{i} = state;
        cache(key) = state;
    end
end

function T = levelTable(states, nxs, nms, dts, kSplit, density)
    % Errors between consecutive levels, split by Hermite mode range.
    n = numel(states) - 1;
    eAll = zeros(n, 1); eLow = zeros(n, 1); eHigh = zeros(n, 1);
    ePsi = zeros(n, 1); eQ = zeros(n, 1);
    for i = 1:n
        e = richardsonModeErrors(states{i+1}, states{i}, 'D', density);
        eAll(i) = norm(e);
        eLow(i) = norm(e(1:min(kSplit, end)));
        eHigh(i) = norm(e(kSplit+1:end));
        ePsi(i) = norm(richardsonModeErrors(states{i+1}, states{i}, 'Psi', density));
        eQ(i) = norm(richardsonModeErrors(states{i+1}, states{i}, 'Q', density));
    end
    level = arrayfun(@(i) sprintf('(%d,%d)->(%d,%d)', nxs(i), nms(i), ...
        nxs(i+1), nms(i+1)), (1:n).', 'UniformOutput', false);
    T = table(level, dts(1:n).', eAll, orders(eAll), eLow, orders(eLow), ...
        eHigh, orders(eHigh), ePsi, orders(ePsi), eQ, orders(eQ), ...
        'VariableNames', {'level', 'dt', 'D_all', 'O_all', ...
        sprintf('D_k_lt_%d', kSplit), 'O_low', ...
        sprintf('D_k_ge_%d', kSplit), 'O_high', 'Psi', 'O_Psi', 'Q', 'O_Q'});
end

function e = richardsonModeErrors(fine, coarse, fieldName, density)
    % Per-component L2 error between fine and coarse (nested, factor 2),
    % coarse components zero-padded. Same computation as
    % Analyzer.addRelative with the 'L2' reduction.
    C = fine.Dofs.(fieldName);
    C0 = coarse.Dofs.(fieldName);
    metric = physics.analysis.FiniteElementRichardsonMetric( ...
        fine.XDisc, coarse.XDisc, reductions={'L2'});
    metric.setPoints(density);
    U = fine.XDisc.eval(metric.Nodes{1}, C(metric.NestLI, :));
    U0 = coarse.XDisc.eval(metric.Nodes{2}, C0);
    n = size(U, 2);
    assert(size(U0, 2) <= n, 'Coarse has more components than fine.');
    U0(:, size(U0, 2)+1:n) = 0;
    e = metric.eval(U, U0);
    e = e(1, :);
end

function e = sameSpaceModeErrors(a, b, fieldName, density)
    % Per-component L2 difference of two states on the same spaces.
    space = a.XDisc;
    metric = physics.analysis.FiniteElementAbsoluteMetric( ...
        space, reductions={'L2'});
    metric.setPoints(density);
    xRef = metric.Nodes;
    e = metric.eval(space.eval(xRef, a.Dofs.(fieldName)), ...
        space.eval(xRef, b.Dofs.(fieldName)));
    e = e(1, :);
end

function nk = modeNorms(state, fieldName, density)
    % ||D_k||_{L2(T)} for every component k.
    space = state.XDisc;
    metric = physics.analysis.FiniteElementAbsoluteMetric( ...
        space, reductions={'L2'});
    metric.setPoints(density);
    U = space.eval(metric.Nodes, state.Dofs.(fieldName));
    nk = metric.eval(U, zeros(size(U)));
    nk = nk(1, :);
end

function o = orders(e)
    o = nan(size(e));
    for i = 2:numel(e)
        o(i) = observedOrder(e(i-1), e(i));
    end
end

function q = observedOrder(e1, e2)
    if e1 > 0 && e2 > 0 && isfinite(e1) && isfinite(e2)
        q = log(e1/e2)/log(2);
    else
        q = NaN;
    end
end

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
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
