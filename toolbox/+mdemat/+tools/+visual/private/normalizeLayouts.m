function layouts = normalizeLayouts(input)
%NORMALIZELAYOUTS Validate all bindings and tile occupancy before allocating.
layouts = struct('Name', {}, 'Size', {}, 'Panels', {}, 'Figure', {}, ...
    'TileSpacing', {}, 'Padding', {});
if ~isstruct(input) || ~(isempty(input) || isvector(input))
    error('mdemat:tools:visual:InvalidLayout', 'Layouts must be a struct vector.');
end
names = strings(1, 0);
for k = 1:numel(input)
    item = input(k);
    if ~all(isfield(item, {'Name', 'Size', 'Panels'}))
        invalid('Each layout needs Name, Size and Panels.');
    end
    name = string(item.Name);
    gridSize = item.Size;
    if ~isscalar(name) || ismissing(name) || ~isvarname(char(name)) || any(names == name)
        invalid('Layout names must be unique MATLAB identifiers.');
    end
    if ~isnumeric(gridSize) || ~isequal(gridSize(:).', gridSize) || numel(gridSize) ~= 2 || ...
            any(~isfinite(gridSize) | gridSize < 1 | gridSize ~= fix(gridSize))
        invalid('Size must be [positiveRows positiveColumns].');
    end
    if ~isstruct(item.Panels) || isempty(item.Panels) || ~isvector(item.Panels)
        invalid('Panels must be a nonempty struct vector.');
    end
    panels = struct('Mode', {}, 'Field', {}, 'ComponentId', {}, 'DatasetIds', {}, ...
        'Tiles', {}, 'Spans', {}, 'Label', {}, 'RenderFn', {});
    occupied = false(gridSize(1), gridSize(2));
    for j = 1:numel(item.Panels)
        p = item.Panels(j);
        mode = string(value(p, 'Mode', '1d'));
        field = string(value(p, 'Field', ''));
        component = value(p, 'ComponentId', 1);
        ids = string(value(p, 'DatasetIds', strings(1, 0)));
        ids = ids(:).';
        fn = value(p, 'RenderFn', []);
        label = string(value(p, 'Label', ''));
        if ~isscalar(label) || ismissing(label), invalid('Label must be one text scalar.'); end
        if ~isscalar(mode) || ~any(mode == ["1d", "2d", "3d", "slice1d", "points", "custom"]) || ...
                ~isscalar(field) || ismissing(field) || ~isnumeric(component) || ~isreal(component) || ~isscalar(component) || ...
                ~isfinite(component) || component < 1 || component ~= fix(component) || ...
                numel(unique(ids)) ~= numel(ids) || any(ismissing(ids)) || ...
                any(~arrayfun(@(s) isvarname(char(s)), ids)) || ...
                ~(isempty(fn) || isa(fn, 'function_handle') && isscalar(fn))
            invalid('Panel mode, field, component, dataset IDs or callback is invalid.');
        end
        if isempty(fn) && (strlength(field) == 0 || isempty(ids) || mode == "custom")
            invalid('Built-in panels need a field and fixed DatasetIds.');
        end
        count = 1;
        if any(mode == ["2d", "3d"]), count = numel(ids); end
        tiles = value(p, 'Tiles', []);
        if ~isempty(fn)
            count = max(1, count);
            if ~isempty(tiles), count = numel(tiles); end
        end
        spans = value(p, 'Spans', [1, 1]);
        if isequal(size(spans), [1, 2]), spans = repmat(spans, count, 1); end
        if ~isnumeric(spans) || ~isequal(size(spans), [count, 2]) || ...
                any(~isfinite(spans(:)) | spans(:) < 1 | spans(:) ~= fix(spans(:)))
            invalid('Spans must be one row or one [rows columns] row per axis.');
        end
        if ~isempty(tiles) && (~isnumeric(tiles) || numel(tiles) ~= count || ...
                any(~isfinite(tiles(:)) | tiles(:) < 1 | tiles(:) ~= fix(tiles(:))))
            invalid('Tiles must contain one tile start per axis.');
        end
        explicit = ~isempty(tiles);
        if ~explicit, tiles = zeros(1, count); else, tiles = tiles(:).'; end
        for n = 1:count
            candidates = 1:prod(gridSize);
            if explicit, candidates = tiles(n); end
            found = false;
            for t = candidates
                row = floor((t - 1) / gridSize(2)) + 1;
                column = mod(t - 1, gridSize(2)) + 1;
                rr = row + (0:spans(n, 1)-1);
                cc = column + (0:spans(n, 2)-1);
                if rr(end) <= gridSize(1) && cc(end) <= gridSize(2) && ~any(occupied(rr, cc), 'all')
                    occupied(rr, cc) = true; tiles(n) = t; found = true; break;
                end
            end
            if ~found, invalid('Panel tiles overlap or exceed the complete layout.'); end
        end
        panels(end + 1) = struct('Mode', mode, 'Field', field, 'ComponentId', component, ...
            'DatasetIds', ids, 'Tiles', tiles, 'Spans', spans, ...
            'Label', label, 'RenderFn', fn); %#ok<AGROW>
    end
    figure = value(item, 'Figure', []);
    if ~isempty(figure) && (~isscalar(figure) || ~isgraphics(figure, 'figure'))
        invalid('Figure must be an existing scalar figure or empty.');
    end
    spacing = string(value(item, 'TileSpacing', 'loose'));
    padding = string(value(item, 'Padding', 'loose'));
    if ~isscalar(spacing) || ismissing(spacing) || ~isscalar(padding) || ismissing(padding)
        invalid('TileSpacing and Padding must be text scalars.');
    end
    if spacing == "normal", spacing = "loose"; end
    if padding == "normal", padding = "loose"; end
    if padding == "none", padding = "tight"; end
    if ~isscalar(spacing) || ~any(spacing == ["loose", "compact", "tight", "none"]) || ...
            ~isscalar(padding) || ~any(padding == ["loose", "compact", "tight"])
        invalid('TileSpacing or Padding is invalid.');
    end
    layouts(end + 1) = struct('Name', name, 'Size', gridSize, 'Panels', panels, 'Figure', figure, ...
        'TileSpacing', spacing, 'Padding', padding); %#ok<AGROW>
    names(end + 1) = name; %#ok<AGROW>
end
end

function result = value(s, name, fallback)
result = fallback;
if isfield(s, name) && ~isempty(s.(name)), result = s.(name); end
end


function invalid(message)
error('mdemat:tools:visual:InvalidLayout', '%s', message);
end
