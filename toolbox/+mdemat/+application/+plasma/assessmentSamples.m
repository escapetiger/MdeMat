function samples = assessmentSamples(scheme, state, coarse, field, time)
%ASSESSMENTSAMPLES Reconstruct physical F on positive joint (x,v) quadrature.
%   Only scalar F is supplied here. D remains an explicit coefficient diagnostic.
%   The measure is dx dv on the FINITE Config.vBBox interval, never the full
%   velocity line or Maxwellian-weighted coefficient norm. Velocity quadrature
%   is Gauss-Legendre with max(8,2*max(number of modes)) points. Its adequacy
%   is an assessment input to check by quadrature refinement, not certification.
%
%   Each State uses its own velocity basis/temperature. AP reconstruction also
%   multiplies its own sqrtRhoInf. Interlevel comparison uses common physical
%   points without modal padding, even when the velocity bases differ.
arguments
    scheme (1,1) mdemat.application.plasma.VlasovScheme
    state (1,1) mdemat.common.fem.KineticState
    coarse = []
    field {mustBeTextScalar} = 'F'
    time {mustBeNumeric,mustBeScalarOrEmpty,mustBeReal,mustBeFinite} = []
end
if ~strcmp(field,'F')
    error('mdemat:plasma:assessmentSamples:Field','Physical assessment supplies scalar F only.');
end
checkState(state);
options = scheme.Observation.Analysis;
hasReference = isfield(options.Exacts,'F') && ~isempty(options.Exacts.F);
samples = [];
if ~hasReference && isempty(coarse), return; end
bounds = scheme.Config.vBBox;
if ~isnumeric(bounds) || ~isreal(bounds) || numel(bounds) ~= 2 || ...
        any(~isfinite(bounds),'all') || bounds(1) >= bounds(2)
    error('mdemat:plasma:assessmentSamples:VelocityDomain', ...
        'Specify an increasing finite one-dimensional velocity interval.');
end
bounds = double(bounds(:).');
modeCount = state.NVDofs;
if ~hasReference
    checkState(coarse);
    modeCount = max(modeCount,coarse.NVDofs);
end
velocityCount = max(8,2*modeCount);
rule = mdemat.common.quadrature.GaussLegendreRule(1);
[velocity,velocityWeights] = rule.generate(mdemat.data.geometry.Orthotope(bounds),velocityCount);
if hasReference
    plan = mdemat.common.fem.samplingPlan(state.XDisc,options.PointsPerDimension);
    values = reconstruct(state,plan.Nodes,velocity,[]);
    [x,v] = ndgrid(plan.Coordinates(1,:),velocity);
    reference = mdemat.application.plasma.referenceDistribution( ...
        options.Exacts.F,x(:).',v(:).',time);
    comparison = 'exact';
    referenceName = 'supplied physical F';
else
    plan = mdemat.common.fem.nestedSamplingPlan( ...
        state.XDisc,coarse.XDisc,options.PointsPerDimension);
    values = reconstruct(state,plan.FineNodes,velocity,plan.FineDofIndices);
    reference = reconstruct(coarse,plan.CoarseNodes,velocity,[]);
    comparison = 'interlevel';
    referenceName = 'coarse physical F';
end
if any(~isfinite(values),'all') || any(~isfinite(reference),'all')
    error('mdemat:plasma:assessmentSamples:NonFinite','Physical assessment samples must be finite.');
end
samples = struct('Values',reshape(values,1,[]),'Reference',reshape(reference,1,[]), ...
    'Weights',kron(velocityWeights(:),plan.Weights(:)).','ComponentIds',1, ...
    'Comparison',comparison,'ReferenceId',sprintf('%s; truncated v=[%.17g,%.17g]', ...
    referenceName,bounds(1),bounds(2)), ...
    'VelocityBounds',bounds,'VelocityNodes',velocity,'SpatialCoordinates',plan.Coordinates, ...
    'Measure','dx dv on the explicitly truncated velocity interval');
end

function checkState(state)
if ~(isa(state,'mdemat.application.plasma.HermiteState') || ...
        isa(state,'mdemat.application.plasma.ApHermiteState')) || ~isscalar(state) || ...
        state.XDisc.NDims ~= 1 || state.VDisc.NDims ~= 1
    error('mdemat:plasma:assessmentSamples:State','Use a one-dimensional Hermite/AP State.');
end
end

function values = reconstruct(state,nodes,velocity,rows)
coefficients = state.Dofs.D;
if ~isempty(rows), coefficients = coefficients(rows,:); end
coefficients = state.XDisc.eval(nodes,coefficients);
if isa(state,'mdemat.application.plasma.ApHermiteState')
    scale = state.Coefs.sqrtRhoInf;
    if ~isempty(rows), scale = scale(rows,:); end
    scale = state.XDisc.eval(nodes,scale);
    if size(scale,2) ~= 1
        error('mdemat:plasma:assessmentSamples:Scale','sqrtRhoInf must be one scalar spatial field.');
    end
    coefficients = coefficients.*scale;
end
values = state.VDisc.eval(velocity,coefficients.');
end
