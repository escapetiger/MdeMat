function values = sampleField(space, coefficients, plan)
%SAMPLEFIELD Evaluate coefficients into components-by-points numeric samples.
%   VALUES = sampleField(SPACE,C,PLAN) uses a samplingPlan's reference nodes.
%   It does not reconstruct kinetic representations or infer component names.
arguments
    space (1,1) mdemat.common.fem.FiniteElementSpace
    coefficients {mustBeNumeric}
    plan (1,1) struct
end
if ~isfield(plan,'Nodes')
    error('mdemat:analysis:InvalidSampling', 'Sampling plan does not match the element count.');
end
if isfield(plan,'ElementCount') && plan.ElementCount ~= space.NMeshElements
    error('mdemat:analysis:InvalidSampling', 'Sampling plan does not match the element count.');
end
values = space.eval(plan.Nodes, coefficients).';
end
