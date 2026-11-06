function validateDatasets(datasets)
%VALIDATEDATASETS Validate value samples without changing their order or types.
if ~isstruct(datasets) || ~isscalar(datasets)
    error('mdemat:tools:visual:InvalidDatasets', 'Datasets must be a scalar struct.');
end
names = fieldnames(datasets);
for k = 1:numel(names)
    entry = datasets.(names{k});
    if ~isvarname(names{k}) || ~isstruct(entry) || ~isscalar(entry) || ...
            ~isequal(sort(fieldnames(entry)), {'Data'; 'Type'}) || ...
            ~(ischar(entry.Type) && isrow(entry.Type) || ...
              isstring(entry.Type) && isscalar(entry.Type)) || ...
            ~any(strcmp(entry.Type, {'fixed', 'exact', 'numeric'})) || ...
            ~isstruct(entry.Data) || ~isscalar(entry.Data)
        error('mdemat:tools:visual:InvalidDataset', ...
            'Dataset %s must contain Type and a scalar Data struct.', names{k});
    end
    fields = fieldnames(entry.Data);
    for j = 1:numel(fields)
        view = entry.Data.(fields{j});
        if ~isstruct(view) || ~isscalar(view) || ...
                ~all(isfield(view, {'Values', 'ComponentIds'})) || ...
                ~isnumeric(view.Values) || ~ismatrix(view.Values)
            error('mdemat:tools:visual:InvalidView', ...
                'Each field must contain numeric C-by-N Values and ComponentIds.');
        end
        [nc, np] = size(view.Values);
        ids = view.ComponentIds;
        if ~isnumeric(ids) || ~isreal(ids) || size(ids, 1) ~= 1 || ...
                numel(ids) ~= nc || any(~isfinite(ids) | ids < 1 | ids ~= fix(ids)) || ...
                numel(unique(ids)) ~= nc
            error('mdemat:tools:visual:InvalidComponentIds', ...
                'ComponentIds must identify the C rows uniquely.');
        end
        tensor = isfield(view, 'Axes');
        points = isfield(view, 'Coordinates');
        if tensor == points
            error('mdemat:tools:visual:InvalidCoordinates', ...
                'Supply exactly one of Axes or Coordinates.');
        end
        if tensor
            axes = view.Axes;
            if ~iscell(axes) || ~isrow(axes) || isempty(axes)
                error('mdemat:tools:visual:InvalidAxes', 'Axes must be a nonempty row cell.');
            end
            count = 1;
            for d = 1:numel(axes)
                a = axes{d};
                if ~isnumeric(a) || ~isreal(a) || ~isvector(a) || isempty(a) || ...
                        any(~isfinite(a(:))) || ...
                        ~(all(diff(a(:)) > 0) || all(diff(a(:)) < 0))
                    error('mdemat:tools:visual:InvalidAxes', ...
                        'Tensor axes must be finite, real and strictly monotonic.');
                end
                count = count * numel(a);
            end
            if count ~= np
                error('mdemat:tools:visual:CoordinateSizeMismatch', ...
                    'Values must follow the full ndgrid product of Axes.');
            end
        else
            x = view.Coordinates;
            if ~isnumeric(x) || ~isreal(x) || ~ismatrix(x) || ...
                    size(x, 1) < 1 || size(x, 2) ~= np || any(~isfinite(x(:)))
                error('mdemat:tools:visual:InvalidCoordinates', ...
                    'Coordinates must be finite real dimension-by-N points.');
            end
        end
        if isfield(view, 'Connectivity')
            e = view.Connectivity;
            if ~isnumeric(e) || ~isreal(e) || ~ismatrix(e) || ...
                    ~ismember(size(e, 2), [2, 3, 4]) || ...
                    any(~isfinite(e(:)) | e(:) < 1 | e(:) > np | e(:) ~= fix(e(:)))
                error('mdemat:tools:visual:InvalidConnectivity', ...
                    'Connectivity must contain valid point indices.');
            end
        end
        for metadata = {'ElementIds', 'TraceIds'}
            if isfield(view, metadata{1})
                ids = view.(metadata{1});
                if ~isnumeric(ids) || ~isreal(ids) || numel(ids) ~= np || ...
                        any(~isfinite(ids(:)) | ids(:) < 1 | ids(:) ~= fix(ids(:)))
                    error('mdemat:tools:visual:InvalidPointIds', ...
                        'Element and trace IDs must identify every point.');
                end
            end
        end
    end
end
end
