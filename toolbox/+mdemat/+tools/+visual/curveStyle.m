function attributes = curveStyle(style, colorIndex, patternIndex, options)
%CURVESTYLE Select native graphics properties from a plotting style value.
%   namedargs2cell(curveStyle(...)) can be passed to plot. Dataset positions
%   should be fixed across frames so missing curves do not change their color.
arguments
    style (1, 1) struct
    colorIndex (1, 1) double {mustBePositive, mustBeInteger, mustBeFinite} = 1
    patternIndex (1, 1) double {mustBePositive, mustBeInteger, mustBeFinite} = 1
    options.kind (1, 1) string {mustBeMember(options.kind, ["line", "scatter"])} = "line"
end
if ~all(isfield(style, {'Palette', 'LineScheme', 'ScatterScheme', 'LineWidth'})) || ...
        ~iscell(style.Palette) || isempty(style.Palette) || ...
        ~iscell(style.LineScheme) || isempty(style.LineScheme) || ...
        ~iscell(style.ScatterScheme) || isempty(style.ScatterScheme)
    error('mdemat:tools:visual:InvalidStyle', 'Supply a nonempty validated plotStyle value.');
end
color = style.Palette{mod(colorIndex - 1, numel(style.Palette)) + 1};
line = 'none'; marker = 'none';
if options.kind == "line"
    line = style.LineScheme{mod(patternIndex - 1, numel(style.LineScheme)) + 1};
else
    marker = style.ScatterScheme{mod(patternIndex - 1, numel(style.ScatterScheme)) + 1};
end
attributes = struct('Color', color, 'LineStyle', line, 'Marker', marker, 'LineWidth', style.LineWidth);
end
