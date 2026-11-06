function plan = samplingPlan(space, pointsPerDimension)
%SAMPLINGPLAN Quadrature coordinates and physical weights as plain values.
%   P = samplingPlan(SPACE,NP) accepts a positive integer per reference axis;
%   a scalar repeats on all axes. Ordering is quadrature-point then element,
%   exactly as SPACE.eval. Coordinates is dimension-by-points, Weights is a
%   row including abs(det J). DG traces are retained, never uniqued.
arguments
    space (1,1) mdemat.common.fem.FiniteElementSpace
    pointsPerDimension {mustBeNumeric, mustBeReal, mustBeVector, mustBeInteger, ...
        mustBePositive, mustBeFinite}
end
dimensions = space.NDims;
counts = double(pointsPerDimension(:).');
if isscalar(counts), counts = repmat(counts, 1, dimensions); end
if numel(counts) ~= dimensions || ...
        ~isa(space.Element.Geometry, 'mdemat.data.geometry.Orthotope')
    error('mdemat:analysis:InvalidSampling', 'Use one quadrature count per orthotope axis.');
end
integrator = space.Element.Volume.copy();
integrator.setPoints(space.Element.Geometry, counts);
determinants = abs(space.Mesh.computeElementJacobianDeterminants());
if isscalar(determinants), determinants = repmat(determinants, space.NMeshElements, 1); end
plan = struct('Nodes', integrator.Nodes, 'ReferenceNodes', integrator.Nodes, ...
    'Coordinates', space.Mesh.collocate(integrator.Nodes), ...
    'Weights', kron(determinants(:), integrator.Weights(:)).', ...
    'PointsPerDimension', counts, 'ElementCount', space.NMeshElements);
end
