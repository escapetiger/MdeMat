function state = createLayout(figure, spec)
%CREATELAYOUT Allocate only a new owned container; roll back partial creation.
container = uipanel('Parent', figure, 'Units', 'normalized', ...
    'Position', [0, 0, 1, 1], 'BorderType', 'none');
try
    layout = tiledlayout(container, spec.Size(1), spec.Size(2), ...
        'TileSpacing', spec.TileSpacing, 'Padding', spec.Padding);
    axes = cell(1, numel(spec.Panels));
    for p = 1:numel(spec.Panels)
        panel = spec.Panels(p);
        axes{p} = gobjects(1, numel(panel.Tiles));
        for a = 1:numel(panel.Tiles)
            axes{p}(a) = nexttile(layout, panel.Tiles(a), panel.Spans(a, :));
        end
    end
    state = struct('Container', container, 'Layout', layout, 'Axes', {axes});
catch exception
    if isgraphics(container), delete(container); end
    rethrow(exception);
end
end
