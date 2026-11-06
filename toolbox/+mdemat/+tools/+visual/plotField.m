function artists = plotField(axes, datasets, field, componentId, options)
%PLOTFIELD Draw an explicitly selected field without owning the supplied axes.
%   Values are components-by-points. Tensor Axes describe sampling nodes in
%   ndgrid order; nonuniform 2-D nodes use interpolated surfaces. Point views
%   are never silently interpolated. NaN denotes missing values; Inf is rejected.
%   The caller owns returned artists. Failure removes only newly made artists.
arguments
    axes {mustBeA(axes, 'matlab.graphics.axis.Axes')}
    datasets struct
    field (1, 1) string
    componentId (1, 1) {mustBeNumeric, mustBePositive, mustBeInteger, mustBeFinite}
    options.mode (1, 1) string {mustBeMember(options.mode, ["1d", "2d", "3d", "slice1d", "points"])} = "1d"
    options.style (1, 1) struct = mdemat.tools.visual.plotStyle()
    options.datasetIds string = string(fieldnames(datasets)).'
    options.label (1, 1) string = ""
end
ids = options.datasetIds(:).';
[group, ~] = mdemat.tools.visual.selectComponent(datasets, field, componentId, datasetIds=ids);
shared = any(options.mode == ["1d", "slice1d", "points"]);
expected = numel(ids);
if shared, expected = 1; end
if numel(axes) ~= expected || any(~isgraphics(axes, 'axes'))
    error('mdemat:tools:visual:InvalidAxes', 'Axes must match the fixed dataset bindings.');
end
style = options.style;
required = fieldnames(mdemat.tools.visual.plotStyle());
if ~all(isfield(style, required))
    error('mdemat:tools:visual:InvalidStyle', 'Supply a complete plotStyle value.');
end
artists = gobjects(0, 1);
previousFigure = get(groot, 'CurrentFigure');
nextPlot = arrayfun(@(ax) string(ax.NextPlot), axes);
guard = onCleanup(@() restoreAxes(axes, nextPlot, previousFigure));
try
    for ax = axes(:).', ax.NextPlot = 'add'; end
    limits = [Inf, -Inf];
    for index = 1:numel(ids)
        id = ids(index);
        ax = axes(min(index, numel(axes)));
        if shared, ax = axes(1); end
        if ~isfield(group, id)
            if ~shared, title(ax, id + ": no data", 'Interpreter', 'none'); end
            continue;
        end
        entry = group.(id);
        view = entry.Data.(field);
        values = displayValues(view.Values, style.ComplexDisplay);
        finite = values(isfinite(values));
        if ~isempty(finite)
            limits = [min(limits(1), min(finite)), max(limits(2), max(finite))];
        end
        kind = "line";
        if strcmp(entry.Type, 'numeric'), kind = "scatter"; end
        curve = mdemat.tools.visual.curveStyle(style, index, index, kind=kind);
        lineStyle = curve.LineStyle;
        color = curve.Color;
        native = namedargs2cell(curve);
        switch options.mode
            case '1d'
                if isfield(view, 'Axes')
                    if numel(view.Axes) ~= 1, incompatible(); end
                    x = view.Axes{1};
                else
                    if size(view.Coordinates, 1) ~= 1, incompatible(); end
                    x = view.Coordinates(1, :);
                    if ~strcmp(lineStyle, 'none') && ~isfield(view, 'Connectivity')
                        error('mdemat:tools:visual:PointOrderRequired', ...
                            'Point-view lines require explicit Connectivity.');
                    end
                end
                y = values;
                if isfield(view, 'Connectivity') && ~strcmp(lineStyle, 'none')
                    e = view.Connectivity;
                    if size(e, 2) ~= 2, incompatible(); end
                    x = [reshape(x(e.'), 2, []); nan(1, size(e, 1))];
                    y = [reshape(y(e.'), 2, []); nan(1, size(e, 1))];
                end
                h = plot(ax, x(:), y(:), native{:}, 'DisplayName', strrep(id, '_', '-'));
                artists = [artists; h(:)]; %#ok<AGROW>
            case '2d'
                if isfield(view, 'Axes')
                    if numel(view.Axes) ~= 2, incompatible(); end
                    [X, Y] = ndgrid(view.Axes{:});
                    u = reshape(values, size(X));
                    if all(size(X) > 1)
                        h = surface(ax, X, Y, zeros(size(X)), u, ...
                            'FaceColor', 'interp', 'EdgeColor', 'none');
                    else
                        h = scatter(ax, X(:), Y(:), 36, u(:), 'filled');
                    end
                else
                    x = view.Coordinates;
                    if size(x, 1) ~= 2 || ~isfield(view, 'Connectivity') || ...
                            size(view.Connectivity, 2) < 3, incompatible(); end
                    h = patch(ax, 'Faces', view.Connectivity, 'Vertices', x.', ...
                        'FaceVertexCData', values(:), 'FaceColor', 'interp', 'EdgeColor', 'none');
                end
                artists = [artists; h(:)]; %#ok<AGROW>
                viewAxis(ax, 2, style);
            case {'3d', 'slice1d'}
                if ~isfield(view, 'Axes') || ~ismember(numel(view.Axes), [2, 3]), incompatible(); end
                [coordinates, u] = tensor(view, values);
                if any(cellfun(@numel, coordinates) < 2)
                    error('mdemat:tools:visual:InsufficientGrid', 'Interpolation needs two nodes on every axis.');
                end
                if options.mode == "3d"
                    if numel(coordinates) ~= 3, incompatible(); end
                    mid = cellfun(@(a) (a(1) + a(end)) / 2, coordinates);
                    h = slice(ax, coordinates{1}, coordinates{2}, coordinates{3}, ...
                        permute(u, [2, 1, 3]), mid(1), mid(2), mid(3));
                    artists = [artists; h(:)]; %#ok<AGROW>
                    set(h, 'EdgeColor', 'none');
                    viewAxis(ax, 3, style);
                else
                    lower = cellfun(@(a) a(1), coordinates);
                    upper = cellfun(@(a) a(end), coordinates);
                    diagonalLength = norm(upper - lower);
                    t = linspace(0, 1, style.SlicePoints);
                    points = lower(:) + (upper(:) - lower(:)) * t;
                    queries = num2cell(points, 2);
                    y = interpn(coordinates{:}, u, queries{:}, 'linear');
                    r = (t - 0.5) * diagonalLength;
                    h = plot(ax, r, y, native{:}, 'DisplayName', strrep(id, '_', '-'));
                    artists = [artists; h(:)]; %#ok<AGROW>
                end
            case 'points'
                if isfield(view, 'Axes')
                    if numel(view.Axes) > 3, incompatible(); end
                    nodes = cell(size(view.Axes)); [nodes{:}] = ndgrid(view.Axes{:});
                    x = cell2mat(cellfun(@(a) a(:).', nodes(:), 'UniformOutput', false));
                else
                    x = view.Coordinates;
                end
                switch size(x, 1)
                    case 1
                        h = scatter(ax, x(1, :), values, 36, color, 'DisplayName', strrep(id, '_', '-'));
                    case 2
                        h = scatter(ax, x(1, :), x(2, :), 36, values, 'filled', 'DisplayName', strrep(id, '_', '-'));
                    case 3
                        h = scatter3(ax, x(1, :), x(2, :), x(3, :), 36, values, 'filled', 'DisplayName', strrep(id, '_', '-'));
                    otherwise
                        incompatible();
                end
                artists = [artists; h(:)]; %#ok<AGROW>
        end
        if ~shared
            title(ax, strrep(id, '_', '-') + ": " + field + string(componentId), ...
                'FontSize', style.TitleFontSize, 'Interpreter', 'none');
        end
    end
    for ax = axes(:).'
        ax.FontSize = style.FontSize;
        ax.XScale = 'linear'; ax.YScale = 'linear';
        if style.LogX, ax.XScale = 'log'; end
        if style.LogY, ax.YScale = 'log'; end
        xlabel(ax, 'x', 'FontSize', style.AxisLabelFontSize);
        label = options.label;
        if strlength(label) == 0, label = field + string(componentId); end
        if shared
            ylabel(ax, label, 'FontSize', style.AxisLabelFontSize, 'Interpreter', 'none');
        else
            ylabel(ax, 'y', 'FontSize', style.AxisLabelFontSize);
        end
        if options.mode == "slice1d", xlabel(ax, 'signed diagonal distance'); end
        if any(options.mode == ["2d", "3d", "points"])
            colormap(ax, style.ColorMap);
            previousBars = findall(ancestor(ax, 'figure'), 'Type', 'colorbar');
            bar = colorbar(ax);
            if ~any(bar == previousBars(:))
                setappdata(bar, 'MdeMatOwnerAxes', ax);
                artists = [artists; bar]; %#ok<AGROW>
            end
            if all(isfinite(limits))
                clim(ax, colorRange(limits, style));
            end
        end
    end
    entries = artists(~isgraphics(artists, 'colorbar'));
    if shared
        if isempty(entries), title(axes(1), 'no data'); else, title(axes(1), ''); end
    end
    if shared && ~isempty(entries) && isempty(axes(1).Legend)
        leg = legend(axes(1), entries, 'Location', style.LegendLocation, 'Interpreter', 'none');
        setappdata(leg, 'MdeMatOwnerAxes', axes(1));
        artists = [artists; leg];
        leg.FontSize = style.LegendFontSize;
    end
catch exception
    delete(artists(isgraphics(artists)));
    rethrow(exception);
end
end

function restoreAxes(axes, nextPlot, previousFigure)
for j = 1:numel(axes)
    if isgraphics(axes(j)), axes(j).NextPlot = char(nextPlot(j)); end
end
if isempty(previousFigure) || isgraphics(previousFigure, 'figure')
    set(groot, 'CurrentFigure', previousFigure);
end
end

function values = displayValues(values, mode)
if ~isreal(values)
    switch mode
        case 'real', values = real(values);
        case 'imag', values = imag(values);
        case 'abs', values = abs(values);
        case 'phase', values = angle(values);
        otherwise
            error('mdemat:tools:visual:ComplexDisplayRequired', 'Choose an explicit complex display.');
    end
end
if issparse(values), values = full(values); end
if any(isinf(values(:)))
    error('mdemat:tools:visual:InfiniteValues', 'Infinite values cannot define a physical display range.');
end
end

function [coordinates, u] = tensor(view, values)
coordinates = view.Axes;
n = cellfun(@numel, coordinates);
u = reshape(full(values), n);
for d = 1:numel(coordinates)
    if coordinates{d}(1) > coordinates{d}(end)
        coordinates{d} = flip(coordinates{d});
        u = flip(u, d);
    end
end
end

function viewAxis(ax, dimension, style)
axis(ax, 'equal', 'tight');
xlabel(ax, 'x', 'FontSize', style.AxisLabelFontSize);
ylabel(ax, 'y', 'FontSize', style.AxisLabelFontSize);
if dimension == 3, zlabel(ax, 'z'); view(ax, 3); else, view(ax, 2); end
end

function limits = colorRange(range, style)
if style.SymmetricColorLimits && range(1) < 0 && range(2) > 0
    limits = [-max(abs(range)), max(abs(range))];
    return;
end
scale = max(abs(range));
width = range(2) - range(1);
if ~isfinite(width), width = scale; end
padding = style.ColorPadding * width;
if padding == 0, padding = max(style.ColorPadding * scale, 0.1); end
limits = [max(-realmax, range(1) - padding), min(realmax, range(2) + padding)];
end

function incompatible()
error('mdemat:tools:visual:IncompatibleView', 'This display mode requires explicit compatible geometry.');
end
