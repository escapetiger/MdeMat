function h = plot(graph, varargin)
%PLOT Draw vertices and edges; return graphics without modifying the graph.
%   plot(graph, axes=ax, ...) is additive unless replace=true. With omitted
%   axes, MATLAB hold-off replacement is retained. See mesh.internal.plotGraph
%   for vertex/edge colors, sizes, visibility and the edge budget.
if ~isscalar(graph)
    error('mdemat:data:mesh:Graph:InvalidGraph', 'Supply a scalar Graph.');
end
h = mdemat.data.mesh.internal.plotGraph(graph.Vertices, graph.Edges, varargin{:});
end
