function h = drawEdges(ax, vertices, edges, color, width)
%DRAWEDGES Draw a graph edge batch once using NaN-separated segments.
h = gobjects(0, 1);
if isempty(edges), return; end
nd = size(vertices, 2);
coordinates = cell(1, max(2, nd));
for d = 1:nd
    coordinates{d} = [vertices(edges(:, 1), d).'; ...
        vertices(edges(:, 2), d).'; nan(1, size(edges, 1))];
    coordinates{d} = coordinates{d}(:);
end
if nd == 1
    coordinates{2} = repmat([0; 0; NaN], size(edges, 1), 1);
end
if nd < 3
    h = plot(ax, coordinates{1}, coordinates{2}, '-', 'Color', color, 'LineWidth', width);
else
    h = plot3(ax, coordinates{1}, coordinates{2}, coordinates{3}, '-', ...
        'Color', color, 'LineWidth', width);
end
end
