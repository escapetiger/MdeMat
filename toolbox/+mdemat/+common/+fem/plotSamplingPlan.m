function plan = plotSamplingPlan(space, pointsPerDimension)
%PLOTSAMPLINGPLAN Describe interior display samples without graphics state.
%   PLAN = plotSamplingPlan(SPACE,N) uses N points in each reference direction;
%   a vector supplies one count per direction. Nodes are dimension-by-points.
%   Coordinates follows SPACE.eval's element ordering. For a verified Grid
%   tensor product, Axes and Permutation reorder those samples into ndgrid
%   column-major order. Repeated/trace coordinates remain an explicit point
%   cloud rather than being merged or averaged.
arguments
    space (1, 1) mdemat.common.fem.FiniteElementSpace
    pointsPerDimension {mustBeNumeric, mustBeReal, mustBeVector, ...
        mustBePositive, mustBeInteger, mustBeFinite}
end
geometry = space.Element.Geometry;
if ~isa(geometry, 'mdemat.data.geometry.Orthotope')
    error('mdemat:fem:plotSamplingPlan:Geometry', ...
        'Interior tensor display sampling requires an orthotope reference geometry.');
end
nd = space.NDims;
counts = double(pointsPerDimension(:).');
if isscalar(counts), counts = repmat(counts, 1, nd); end
if numel(counts) ~= nd
    error('mdemat:fem:plotSamplingPlan:Dimension', ...
        'Supply one point count or one per spatial direction.');
end
axes = cell(1, nd);
for direction = 1:nd
    fractions = (1:counts(direction))/(counts(direction)+1);
    axes{direction} = geometry.Lower(direction) + ...
        (geometry.Upper(direction)-geometry.Lower(direction))*fractions;
end
grids = cell(1, nd);
[grids{:}] = ndgrid(axes{:});
nodes = zeros(nd, prod(counts));
for direction = 1:nd, nodes(direction, :) = grids{direction}(:).'; end
coordinates = reshape(space.Mesh.collocate(nodes), nd, []);
plan = struct('Nodes', nodes, 'Coordinates', coordinates, ...
    'Permutation', 1:size(coordinates, 2));
if ~isa(space.Mesh, 'mdemat.data.mesh.Grid'), return; end

% Grid topology supplies the tensor claim; the point-to-grid bijection is
% checked independently before using a compact Axes representation.
physicalAxes = cell(1, nd);
indices = cell(1, nd);
shape = zeros(1, nd);
for direction = 1:nd
    physicalAxes{direction} = unique(coordinates(direction, :), 'sorted');
    shape(direction) = numel(physicalAxes{direction});
    [found, indices{direction}] = ismember(coordinates(direction, :), physicalAxes{direction});
    if ~all(found), return; end
end
if prod(shape) ~= size(coordinates, 2), return; end
linear = ones(1, size(coordinates, 2));
stride = 1;
for direction = 1:nd
    linear = linear + (indices{direction}-1)*stride;
    stride = stride*shape(direction);
end
if numel(unique(linear)) ~= numel(linear), return; end
[~, permutation] = sort(linear);
plan.Axes = physicalAxes;
plan.Permutation = permutation;
end
