%==========================================================================
% FileName: ex01_time_convergence.m
% Description:
%   Temporal convergence test for the ex01ldh AP-DG-Hermite Vlasov problem.
%
%   The spatial mesh and Hermite truncation are fixed. Only dt is refined.
%   Temporal errors are estimated from consecutive solutions with dt and
%   dt/2, using the same final time and the same discrete spatial space.
%==========================================================================

clc;
clear;
close all;

%% PROBLEM CONFIGURATION
common = fullfile(fileparts(mfilename('fullpath')), 'config.txt');
config0 = physics.scheme.ConfigParser.parseFile(common);

% Match ex01ldh in ex01_convergence.m.
config0.id = 'ex01ldh';
config0.nDims = 1;
config0.L = 4*pi;
config0.xBBox = [0, config0.L];
config0.kappa = 2*pi/config0.L;
config0.T0 = 1;
config0.delta = 0.5;
config0.nx = 16;
config0.nh = 16;
config0.epsilon = 1;
config0.tau0 = 1e5;
config0.tau = config0.tau0 * config0.epsilon;
config0.tFinal = 1;
config0.vBBox = [-3.2, 3.2];
config0.ckptTimeStamps = [];
config0.bc = [];
config0.exact = [];
config0.verbose = 0;

% Disable CFL explicitly: only the fixed time step below controls dt.
config0.cfl = [];

% Each dt divides T=1 exactly. The last entry is an additional finer run
% needed to estimate the order for the smallest dt in dtList.
dtList = [0.02, 0.01, 0.005, 0.0025, 0.00125];
dtAll = [dtList, dtList(end)/2];

% Build the same equilibrium and initial data as ex01_convergence.m.
config0.phiInf = phiInf(config0);
config0.EInf = EInf(config0);
config0.cInf = cInf(config0);
config0.rhoInf = rhoInf(config0);
config0.sqrtRhoInf = sqrtRhoInf(config0);
config0.dSqrtRhoInf = dSqrtRhoInf(config0);
config0.ic = DInit(config0);

%% RUN WITH SUCCESSIVELY REFINED TIME STEPS
states = cell(numel(dtAll), 1);

for j = 1:numel(dtAll)
    cfg = config0;
    cfg.dt = dtAll(j);
    cfg.cfl = [];

    fprintf('\n[M] Time-convergence run %d/%d: dt = %.8g\n', ...
        j, numel(dtAll), dtAll(j));

    [scheme, state] = buildProblem(cfg);
    states{j} = scheme.run(state);

    % Check that each run reached the same final time.
    tReached = scheme.TDisc.Timeline.Now;
    tolT = 100*eps(max(1, cfg.tFinal));
    assert(abs(tReached - cfg.tFinal) <= tolT, ...
        'Time integration did not reach the prescribed final time.');
end

%% ESTIMATE TEMPORAL ERRORS AND ORDERS
% E(dt) = || U_dt - U_(dt/2) ||.
% If the method is in its asymptotic temporal regime, E(dt) = O(dt^q),
% and E(dt)/E(dt/2) should approach 2^q.
nDt = numel(dtList);
errD = zeros(nDt, 1);
errPsi = zeros(nDt, 1);
errQ = zeros(nDt, 1);

for j = 1:nDt
    errD(j) = spatialL2SequenceError(states{j}, states{j+1}, 'D');
    errPsi(j) = spatialL2SequenceError(states{j}, states{j+1}, 'Psi');
    errQ(j) = spatialL2SequenceError(states{j}, states{j+1}, 'Q');
end

orderD = nan(nDt, 1);
orderPsi = nan(nDt, 1);
orderQ = nan(nDt, 1);

for j = 1:nDt-1
    orderD(j) = observedOrder(errD(j), errD(j+1));
    orderPsi(j) = observedOrder(errPsi(j), errPsi(j+1));
    orderQ(j) = observedOrder(errQ(j), errQ(j+1));
end

result = table(dtList(:), errD, orderD, errPsi, orderPsi, errQ, orderQ, ...
    'VariableNames', {'dt', 'D_L2_sequence_Error', 'D_Order', ...
    'Psi_L2_Error', 'Psi_Order', 'Q_L2_Error', 'Q_Order'});

disp(result);

outFile = fullfile(fileparts(mfilename('fullpath')), ...
    'ex01_time_convergence.csv');
writetable(result, outFile);
fprintf('\n[M] Results written to %s\n', outFile);

%% LOCAL FUNCTIONS

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

function e = spatialL2SequenceError(state, ref, fieldName)
    % Compute spatial L2 error for each column, then aggregate columns in l2.
    % State and reference must use the same spatial and Hermite spaces.
    space = state.XDisc;
    C = state.Dofs.(fieldName);
    CRef = ref.Dofs.(fieldName);

    assert(size(C, 1) == size(CRef, 1) && ...
        size(C, 2) == size(CRef, 2), ...
        'State and reference dimensions differ for field %s.', fieldName);

    metric = physics.analysis.FiniteElementAbsoluteMetric( ...
        space, reductions={'L2'});
    metric.setPoints(8);
    xRef = metric.Nodes;

    U = space.eval(xRef, C);
    URef = space.eval(xRef, CRef);
    R = metric.eval(U, URef);

    % R contains one spatial L2 error for each component/Hermite mode.
    e = norm(R(1, :), 2);
end

function q = observedOrder(e1, e2)
    if e1 > 0 && e2 > 0 && isfinite(e1) && isfinite(e2)
        q = log(e1/e2)/log(2);
    else
        q = NaN;
    end
end

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
    % Preserve the normalization convention in the original example.
    c = 2*config.L / integral( ...
        @(x) exp(-phiInfImpl(config, x)/config.T0), -config.L, config.L);
end
