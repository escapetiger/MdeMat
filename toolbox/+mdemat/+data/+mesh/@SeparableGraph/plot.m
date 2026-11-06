function h = plot(graph, options)
%PLOT Draw bounded tensor factors; preserve endpoints and induced edges.
%   Projection selects distinct original dimensions. Sampling never invents
%   an edge between retained vertices. Too small an endpoint budget is rejected.
arguments
    graph (1, 1) mdemat.data.mesh.SeparableGraph
    options.axes = []
    options.vertexColor = 'b'
    options.edgeColor = 'k'
    options.vertexSize (1, 1) double {mustBePositive, mustBeFinite} = 36
    options.edgeWidth (1, 1) double {mustBePositive, mustBeFinite} = 1
    options.showVertices (1, 1) logical = true
    options.showEdges (1, 1) logical = true
    options.maxVertices (1, 1) double {mustBePositive, mustBeInteger, mustBeFinite} = 10000
    options.maxEdges (1, 1) double {mustBeNonnegative} = 100000
    options.projectionDims double {mustBePositive, mustBeInteger, mustBeFinite} = []
    options.replace logical {mustBeScalarOrEmpty} = []
end
if ~isempty(options.projectionDims) && ~isvector(options.projectionDims)
    error('mdemat:data:mesh:SeparableGraph:InvalidProjection', 'Projection dimensions must be a vector.');
end
dims = options.projectionDims(:).';
if isempty(dims), dims = 1:min(graph.NDims, 3); end
if numel(dims) > 3 || numel(unique(dims)) ~= numel(dims) || any(dims > graph.NDims)
    error('mdemat:data:mesh:SeparableGraph:InvalidProjection', 'Choose one to three distinct existing dimensions.');
end
v = graph.Vertices; e = graph.Edges; v = v(dims); e = e(dims);
counts = cellfun(@numel, v); minimum = min(2, counts);
if prod(minimum) > options.maxVertices
    error('mdemat:data:mesh:SeparableGraph:VertexBudgetTooSmall', 'The budget cannot retain all factor endpoints.');
end
original = counts;
if prod(counts) > options.maxVertices
    ratio = (options.maxVertices / prod(counts)) ^ (1 / numel(counts));
    counts = max(minimum, floor(counts * ratio));
    while prod(counts) > options.maxVertices
        candidates = find(counts > minimum);
        [~, index] = max(counts(candidates)); d = candidates(index);
        counts(d) = counts(d) - 1;
    end
end
for d = 1:numel(v)
    indices = unique(round(linspace(1, original(d), counts(d))));
    v{d} = v{d}(indices);
    if isempty(e{d})
        e{d} = zeros(0, 2);
    else
        [~, a] = ismember(e{d}(:, 1), indices); [~, b] = ismember(e{d}(:, 2), indices);
        keep = a > 0 & b > 0;
        e{d} = [a(keep), b(keep)];
    end
end
counts = cellfun(@numel, v);
edgeCount = sum(cellfun(@(a) size(a, 1), e) .* (prod(counts) ./ counts));
if options.showEdges && edgeCount > options.maxEdges
    error('mdemat:data:mesh:SeparableGraph:EdgeBudgetExceeded', 'Retained induced edges exceed maxEdges.');
end
nodes = cell(size(v)); [nodes{:}] = ndgrid(v{:});
vertices = cell2mat(cellfun(@(x) x(:), nodes, 'UniformOutput', false));
edges = zeros(0, 2);
if options.showEdges
    edges = zeros(edgeCount, 2); offset = 0; strides = [1, cumprod(counts(1:end-1))];
    for d = 1:numel(v)
        other = setdiff(1:numel(v), d); indices = cell(1, numel(other));
        if isempty(other)
            base = 1;
        else
            args = arrayfun(@(n) 1:n, counts(other), 'UniformOutput', false);
            [indices{:}] = ndgrid(args{:}); base = ones(prod(counts(other)), 1);
            for j = 1:numel(other), base = base + (indices{j}(:) - 1) * strides(other(j)); end
        end
        a = base(:) + (e{d}(:, 1).' - 1) * strides(d);
        b = base(:) + (e{d}(:, 2).' - 1) * strides(d);
        batch = [a(:), b(:)];
        edges(offset + (1:size(batch, 1)), :) = batch; offset = offset + size(batch, 1);
    end
end
% Ordinary and factor graphs share one numeric drawing implementation.
h = mdemat.data.mesh.internal.plotGraph(vertices, edges, axes=options.axes, vertexColor=options.vertexColor, edgeColor=options.edgeColor, ...
    vertexSize=options.vertexSize, edgeWidth=options.edgeWidth, showVertices=options.showVertices, ...
    showEdges=options.showEdges, showReference=false, replace=options.replace, maxEdges=options.maxEdges);
if isempty(options.axes), ax = gca; else, ax = options.axes; end
title(ax, sprintf('%dD SeparableGraph: %d sampled vertices; projected dims %s', ...
    graph.NDims, size(vertices, 1), mat2str(dims)));
end
