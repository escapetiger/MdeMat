function [result, scheme] = solve(record, options)
%SOLVE Advance the complete macro/micro DG radiation model without observers.
%   record is produced by radiation.problem. All original decomposition,
%   scattering, absorption, source, boundary, filter and positivity kernels
%   remain available through record.config. No reference solution is needed.
arguments
    record (1, 1) struct
    options.observe = []
    options.callback = []
    options.maxSteps (1, 1) {mustBePositive, mustBeInteger} = 1000000
end
fields = {'schemaVersion'; 'model'; 'spatialSpace'; 'velocitySpace'; ...
    'initial'; 'finalTime'; 'timeStep'; 'config'};
assert(isequal(sort(fieldnames(record)), sort(fields)) && ...
    isequal(record.model, 'linear-radiation') && isequal(record.schemaVersion, 1) && ...
    isstruct(record.config) && isscalar(record.config), ...
    'MdeMat:Radiation:InvalidProblem', 'Expected a radiation problem record.');
record = mdemat.application.radiation.problem(record.spatialSpace, record.velocitySpace, ...
    record.initial, record.finalTime, record.timeStep, config=record.config);
config = record.config;
config.ic = record.initial;
config.tFinal = record.finalTime;
config.dt = record.timeStep;
config.cfl = [];
config.nDims = record.spatialSpace.NDims;
config.xBasisType = record.spatialSpace.Element.Approximator.Type;
if ~isfield(config, 'xBasisOrder')
    config.xBasisOrder = record.spatialSpace.NLocalDofs;
end
space = record.velocitySpace;
if space.NLhsDofs > 0
    geometry = space.Lhs.Element.Geometry;
    reduction = space.Lhs.Element.DimReduction;
else
    geometry = space.Rhs.Element.Geometry;
    reduction = space.Rhs.Element.DimReduction;
end
if ~isfield(config, 'vDimReduction')
    config.vDimReduction = reduction;
end
if ~isfield(config, 'E')
    config.E = 1/geometry.magnitude();
end
if space.NLhsDofs == 0
    config.nu = 0;
elseif config.nDims == 1
    config.nu = space.NLhsDofs;
elseif config.nDims == 2 && strcmp(config.vDimReduction, 'topology')
    config.nu = (space.NLhsDofs+1)/2;
else
    config.nu = sqrt(space.NLhsDofs);
end
config.nv = space.NRhsDofs;
assert(config.nu == fix(config.nu), 'MdeMat:Radiation:InvalidMacroSpace', ...
    'The macro space must contain complete degree blocks for this geometry.');
scheme = mdemat.application.radiation.MmDgScheme(config=config);
ownedVelocitySpace = mdemat.common.fem.SumSpace(space.Lhs, space.Rhs, ...
    lhsWeight=space.LhsWeight, rhsWeight=space.RhsWeight);
state = mdemat.application.radiation.MacroMicroState(record.spatialSpace, ownedVelocitySpace);
result = mdemat.internal.execution.compute(scheme, state, observe=options.observe, ...
    callback=options.callback, maxSteps=options.maxSteps);
end
