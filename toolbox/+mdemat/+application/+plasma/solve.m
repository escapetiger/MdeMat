function [result, scheme] = solve(record, options)
%SOLVE Advance the original self-consistent DG/Hermite plasma equations.
%   Both APDGHS and DGHS use their original initialization, Poisson coupling,
%   linear/collisional and nonlinear splitting, filtering and diagnostics.
%   observe and callback are optional explicit compute-time capabilities.
arguments
    record (1, 1) struct
    options.observe = []
    options.callback = []
    options.maxSteps (1, 1) {mustBePositive, mustBeInteger} = 1000000
end
fields = {'schemaVersion'; 'model'; 'spatialSpace'; 'velocitySpace'; ...
    'initial'; 'finalTime'; 'timeStep'; 'method'; 'config'};
assert(isequal(sort(fieldnames(record)), sort(fields)) && ...
    isequal(record.model, 'vlasov-poisson-fokker-planck') && ...
    isequal(record.schemaVersion, 1) && isstruct(record.config) && isscalar(record.config), ...
    'MdeMat:Plasma:InvalidProblem', 'Expected a plasma problem record.');
record = mdemat.application.plasma.problem(record.spatialSpace, record.velocitySpace, ...
    record.initial, record.finalTime, record.timeStep, method=record.method, config=record.config);
config = record.config;
config.ic = record.initial;
config.nh = record.velocitySpace.NDofs;
config.xBasisType = record.spatialSpace.Element.Approximator.Type;
config.xBasisOrder = record.spatialSpace.NLocalDofs;
config.tFinal = record.finalTime;
config.dt = record.timeStep;
config.cfl = [];
if ~isfield(config, 'T0')
    config.T0 = record.velocitySpace.Element.Approximator.Basis.Rhs.Temperature;
end
switch record.method
    case 'ap'
        scheme = mdemat.application.plasma.ApDghsScheme(config=config);
        state = mdemat.application.plasma.ApHermiteState(record.spatialSpace, record.velocitySpace);
    case 'dghs'
        scheme = mdemat.application.plasma.DghsScheme(config=config);
        state = mdemat.application.plasma.HermiteState(record.spatialSpace, record.velocitySpace);
end
result = mdemat.internal.execution.compute(scheme, state, observe=options.observe, ...
    callback=options.callback, maxSteps=options.maxSteps);
end
