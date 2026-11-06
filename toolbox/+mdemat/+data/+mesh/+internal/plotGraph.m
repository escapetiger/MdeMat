function h = plotGraph(v, e, options)
%PLOT Draw a graph without attaching graphics state to its data.
%   H contains Vertices/Edges/Reference graphics. Explicit axes are additive
%   unless replace=true. With no axes, ordinary hold-off replacement applies.
arguments
    v {mustBeNumeric}
    e {mustBeNumeric}
    options.axes = []
    options.vertexColor = 'b'
    options.edgeColor = 'k'
    options.vertexSize (1, 1) double {mustBePositive, mustBeFinite} = 36
    options.edgeWidth (1, 1) double {mustBePositive, mustBeFinite} = 1
    options.showVertices (1, 1) logical = true
    options.showEdges (1, 1) logical = true
    options.showReference (1, 1) logical = true
    options.replace logical {mustBeScalarOrEmpty} = []
    options.maxEdges (1, 1) double {mustBeNonnegative} = Inf
end
h = struct('Vertices', gobjects(0, 1), 'Edges', gobjects(0, 1), 'Reference', gobjects(0, 1));
if isempty(v) && isempty(e), return; end
if ~isnumeric(v) || ~isreal(v) || ~ismatrix(v) || any(~isfinite(v(:))) || ...
        ~ismember(size(v, 2), 1:3)
    error('mdemat:data:mesh:Graph:UnsupportedDimension', 'Plot finite real vertex matrices in one to three dimensions.');
end
if ~isnumeric(e) || ~isreal(e) || ~ismatrix(e) || ...
        ~isempty(e) && size(e, 2) ~= 2 || ...
        any(~isfinite(e(:)) | e(:) ~= fix(e(:)) | e(:) < 1 | e(:) > size(v, 1))
    error('mdemat:data:mesh:Graph:InvalidEdges', 'Edges must contain valid vertex pairs.');
end
if options.showEdges && size(e, 1) > options.maxEdges
    error('mdemat:data:mesh:Graph:EdgeBudgetExceeded', 'The requested graph exceeds maxEdges.');
end
explicit = ~isempty(options.axes);
if explicit, ax = options.axes; else, ax = gca; end
if ~isscalar(ax) || ~isgraphics(ax, 'axes')
    error('mdemat:data:mesh:Graph:InvalidAxes', 'Supply one valid axes.');
end
next = ax.NextPlot; previous = get(groot, 'CurrentFigure');
guard = onCleanup(@() restoreAxes(ax, next, previous));
try
    replace = options.replace;
    if isempty(replace), replace = ~explicit && ~ishold(ax); end
    if replace, cla(ax); end
    ax.NextPlot = 'add';
    if options.showEdges && ~strcmp(options.edgeColor, 'none')
        h.Edges = mdemat.data.mesh.internal.drawEdges(ax, v, e, options.edgeColor, options.edgeWidth);
    end
    if options.showVertices && ~strcmp(options.vertexColor, 'none')
        switch size(v, 2)
            case 1, h.Vertices = scatter(ax, v(:, 1), zeros(size(v, 1), 1), options.vertexSize, options.vertexColor, 'filled');
            case 2, h.Vertices = scatter(ax, v(:, 1), v(:, 2), options.vertexSize, options.vertexColor, 'filled');
            case 3, h.Vertices = scatter3(ax, v(:, 1), v(:, 2), v(:, 3), options.vertexSize, options.vertexColor, 'filled');
        end
    end
    if size(v, 2) == 1
        ylim(ax, [-0.5, 0.5]);
        if options.showReference && (~options.showEdges || strcmp(options.edgeColor, 'none'))
            range = xlim(ax);
            h.Reference = plot(ax, range, [0, 0], '-', 'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.5);
        end
    end
    xlabel(ax, 'X');
    if size(v, 2) >= 2, ylabel(ax, 'Y'); end
    if size(v, 2) == 3, zlabel(ax, 'Z'); view(ax, 3); end
    axis(ax, 'equal'); grid(ax, 'on');
    title(ax, sprintf('%dD Graph: %d vertices, %d edges', size(v, 2), size(v, 1), size(e, 1)));
catch exception
    parts = struct2cell(h);
    for k = 1:numel(parts), delete(parts{k}(isgraphics(parts{k}))); end
    rethrow(exception);
end
end

function restoreAxes(ax, next, previous)
if isgraphics(ax), ax.NextPlot = next; end
if isempty(previous) || isgraphics(previous, 'figure'), set(groot, 'CurrentFigure', previous); end
end
