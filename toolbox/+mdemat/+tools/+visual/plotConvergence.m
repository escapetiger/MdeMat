function [artists, fit] = plotConvergence(ax, scales, errors, options)
%PLOTCONVERGENCE Draw scalar error data and an explicitly normalized order line.
%   referenceOrder=q draws E(k)*(h/h(k))^q over the inclusive fitRange, anchored
%   at its first positive assessed error. Otherwise fitted=true draws the log
%   least-squares line. The fitted Order and LogAmplitude are returned in FIT.
%   Zero/NaN errors are retained by the caller but omitted from this log view;
%   negative/Inf errors are rejected. No State, table or model is interpreted.
%   Only new artists are returned. Temporary hold/current-figure state is kept.
arguments
    ax (1, 1) matlab.graphics.axis.Axes
    scales double {mustBeReal, mustBeVector, mustBeFinite, mustBePositive}
    errors double {mustBeReal, mustBeVector}
    options.referenceOrder double {mustBeReal, mustBeFinite, mustBeScalarOrEmpty} = []
    options.fitted (1, 1) logical = false
    options.fitRange double = []
    options.style (1, 1) struct = mdemat.tools.visual.plotStyle()
    options.colorIndex (1, 1) double {mustBeInteger, mustBePositive, mustBeFinite} = 1
    options.label (1, 1) string = "Error"
end
if ~isgraphics(ax, 'axes') || numel(scales) ~= numel(errors) || ...
        any(isinf(errors) | errors < 0, 'all')
    error('mdemat:tools:visual:InvalidConvergenceData', ...
        'Use matching positive finite scales and nonnegative finite or NaN errors.');
end
h = scales(:); e = errors(:); valid = isfinite(e) & e > 0;
range = options.fitRange;
if isempty(range), range = [1, numel(h)]; end
if ~isreal(range) || numel(range) ~= 2 || any(~isfinite(range) | range ~= fix(range)) || ...
        range(1) < 1 || range(2) > numel(h) || range(1) > range(2)
    error('mdemat:tools:visual:InvalidFitRange', 'fitRange is an inclusive [first last] index pair.');
end
indices = (range(1):range(2)).'; indices = indices(valid(indices));
fit = struct('Order', NaN, 'LogAmplitude', NaN, 'FitIndices', indices, ...
    'ReferenceOrder', NaN, 'Anchor', [NaN, NaN]);
artists = struct('Curve', gobjects(0, 1), 'Reference', gobjects(0, 1), 'Legend', gobjects(0, 1));
x = []; offset = NaN;
if numel(indices) >= 2
    x = logRatio(h(indices), h(indices(1)));
    y = logRatio(e(indices), e(indices(1)));
    centered = x - mean(x); denominator = sum(centered .^ 2);
    if denominator > 0
        fit.Order = sum(centered .* (y - mean(y))) / denominator;
        offset = mean(y) - fit.Order * mean(x);
        fit.LogAmplitude = log(e(indices(1))) + offset - fit.Order * log(h(indices(1)));
    end
end
reference = [];
if ~isempty(options.referenceOrder) || options.fitted
    if ~isfinite(fit.Order)
        error('mdemat:tools:visual:InsufficientFitPoints', ...
            'An order line requires positive assessed errors at two distinct scales.');
    end
    order = options.referenceOrder;
    if isempty(order)
        order = fit.Order;
    else
        offset = 0; % A prescribed order is normalized to an actual data point.
    end
    anchorLog = log(e(indices(1))) + offset;
    reference = exp(anchorLog + order * x);
    if any(~isfinite(reference) | reference <= 0)
        error('mdemat:tools:visual:UnrepresentableReference', ...
            'The requested order line is not representable as positive finite values.');
    end
    fit.ReferenceOrder = order;
    fit.Anchor = [h(indices(1)), exp(anchorLog)];
end
if ~any(valid), return; end
nextPlot = ax.NextPlot; previous = get(groot, 'CurrentFigure');
guard = onCleanup(@() restoreAxes(ax, nextPlot, previous));
try
    ax.NextPlot = 'add';
    style = options.style;
    attributes = mdemat.tools.visual.curveStyle(style, options.colorIndex, options.colorIndex);
    attributes.Marker = style.ScatterScheme{mod(options.colorIndex - 1, numel(style.ScatterScheme)) + 1};
    args = namedargs2cell(attributes);
    artists.Curve = plot(ax, h(valid), e(valid), args{:}, 'DisplayName', options.label);
    if ~isempty(reference)
        artists.Reference = plot(ax, h(indices), reference, '--', 'Color', attributes.Color, ...
            'LineWidth', style.LineWidth, 'DisplayName', "reference p=" + string(fit.ReferenceOrder));
    end
    ax.XScale = 'log'; ax.YScale = 'log'; ax.FontSize = style.FontSize;
    xlabel(ax, 'Resolution', 'FontSize', style.AxisLabelFontSize);
    ylabel(ax, 'Error', 'FontSize', style.AxisLabelFontSize);
    grid(ax, 'on');
    if isempty(ax.Legend)
        lines = [artists.Curve; artists.Reference];
        artists.Legend = legend(ax, lines, 'Location', style.LegendLocation, 'Interpreter', 'none');
        setappdata(artists.Legend, 'MdeMatOwnerAxes', ax);
        artists.Legend.FontSize = style.LegendFontSize;
    end
catch exception
    values = struct2cell(artists);
    for k = 1:numel(values), delete(values{k}(isgraphics(values{k}))); end
    rethrow(exception);
end
end

function result = logRatio(a, b)
% Avoid log(a)-log(b) cancellation near b and a/b overflow far from b.
result = log(a) - log(b);
near = abs(a - b) <= b / 2;
result(near) = log1p((a(near) - b) / b);
end

function restoreAxes(ax, nextPlot, previous)
if isgraphics(ax, 'axes'), ax.NextPlot = nextPlot; end
if isempty(previous) || isgraphics(previous, 'figure'), set(groot, 'CurrentFigure', previous); end
end
