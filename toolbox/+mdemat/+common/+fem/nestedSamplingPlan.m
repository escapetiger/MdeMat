function plan = nestedSamplingPlan(fine, coarse, pointsPerDimension)
%NESTEDSAMPLINGPLAN Common physical quadrature on nested tensor grids.
%   P = nestedSamplingPlan(FINE,COARSE,NP) supports an integer refinement ratio
%   per axis, including equal grids. NP is the quadrature count per fine-cell
%   axis. FineNodes/CoarseNodes sample the same physical points, ordered by
%   fine cell within coarse cell, then coarse cell. FineDofIndices reorders
%   the fine coefficient rows for this ordering. Weights includes coarse-cell
%   Jacobians and subcell scaling. This computes no Richardson correction.
arguments
    fine (1,1) mdemat.common.fem.FiniteElementSpace
    coarse (1,1) mdemat.common.fem.FiniteElementSpace
    pointsPerDimension {mustBeNumeric, mustBeReal, mustBeVector, mustBeInteger, ...
        mustBePositive, mustBeFinite}
end
if ~isa(fine.Mesh,'mdemat.data.mesh.Grid') || ~isa(coarse.Mesh,'mdemat.data.mesh.Grid') || ...
        fine.NDims ~= coarse.NDims || ...
        ~isa(fine.Element.Geometry,'mdemat.data.geometry.Orthotope') || ...
        ~isa(coarse.Element.Geometry,'mdemat.data.geometry.Orthotope')
    error('mdemat:analysis:NonNestedMesh', 'Nested sampling requires equal-dimensional tensor grids.');
end
ratio = fine.Mesh.Resolution./coarse.Mesh.Resolution;
if any(ratio < 1 | ratio ~= fix(ratio))
    error('mdemat:analysis:NonNestedMesh', 'Fine resolution must be an integer multiple per axis.');
end
for dimension = 1:fine.NDims
    axis = coarse.Mesh.Vertices{dimension};
    expected = zeros(1, (numel(axis)-1)*ratio(dimension)+1);
    for cellIndex = 1:numel(axis)-1
        expected((cellIndex-1)*ratio(dimension)+(1:ratio(dimension))) = ...
            axis(cellIndex)+(axis(cellIndex+1)-axis(cellIndex))*(0:ratio(dimension)-1)/ratio(dimension);
    end
    expected(end) = axis(end);
    actual = fine.Mesh.Vertices{dimension};
    if numel(actual) ~= numel(expected) || ...
            any(abs(actual-expected) > 64*eps(max(1,max(abs(expected)))))
        error('mdemat:analysis:NonNestedMesh', 'Fine vertices must subdivide the supplied coarse cells.');
    end
end
finePlan = mdemat.common.fem.samplingPlan(fine, pointsPerDimension);
nodes = finePlan.ReferenceNodes;
weights = fine.Element.Volume.copy();
weights.setPoints(fine.Element.Geometry, finePlan.PointsPerDimension);
referenceWeights = weights.Weights(:).';
childIndexer = mdemat.data.index.MultiIndexer(shape=ratio);
children = childIndexer.linearToMulti((1:prod(ratio)).');
coarseIndexer = mdemat.data.index.MultiIndexer(shape=coarse.Mesh.Resolution);
fineIndexer = mdemat.data.index.MultiIndexer(shape=fine.Mesh.Resolution);
coarseCells = coarseIndexer.linearToMulti((1:coarse.NMeshElements).');
coarseNodes = zeros(fine.NDims, size(nodes,2)*prod(ratio));
elementIndices = zeros(fine.NMeshElements,1);
% Reference orthotopes have axis intervals supplied by the geometry itself.
bounds = fine.Element.Geometry.BBox;
bounds = bounds(:);
lower = bounds(1:2:end);
width = bounds(2:2:end)-bounds(1:2:end);
for child = 1:size(children,1)
    columns = (child-1)*size(nodes,2)+(1:size(nodes,2));
    coarseNodes(:,columns) = lower + ((children(child,:).'-1).*width+nodes-lower)./ratio(:);
end
for element = 1:coarse.NMeshElements
    multi = (coarseCells(element,:)-1).*ratio+children;
    elementIndices((element-1)*prod(ratio)+(1:prod(ratio))) = fineIndexer.multiToLinear(multi);
end
local = fine.NLocalDofs;
dofIndices = reshape((elementIndices(:).'-1)*local+(1:local).',[],1);
determinants = abs(coarse.Mesh.computeElementJacobianDeterminants());
if isscalar(determinants), determinants = repmat(determinants,coarse.NMeshElements,1); end
physicalWeights = kron(determinants(:), repmat(referenceWeights(:),prod(ratio),1)/prod(ratio));
coordinates = coarse.Mesh.collocate(coarseNodes);
plan = struct('FineNodes',nodes,'CoarseNodes',coarseNodes, ...
    'FineDofIndices',dofIndices,'FineElementIndices',elementIndices, ...
    'Weights',physicalWeights.','Coordinates',coordinates,'RefinementRatio',ratio, ...
    'PointsPerDimension',finePlan.PointsPerDimension);
end
