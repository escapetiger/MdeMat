function [group, missing] = selectComponent(datasets, field, componentId, options)
%SELECTCOMPONENT Select original component IDs from value samples.
%   GROUP preserves dataset keys, point order and the per-field view layout.
%   Missing fields/components are omitted, or rejected for required datasets.
arguments
    datasets struct
    field (1, 1) string
    componentId (1, 1) {mustBeNumeric, mustBePositive, mustBeInteger, mustBeFinite}
    options.datasetIds string = string(fieldnames(datasets)).'
    options.required string = strings(1, 0)
end
validateDatasets(datasets);
ids = options.datasetIds(:).';
if numel(unique(ids)) ~= numel(ids) || any(ismissing(ids)) || ...
        any(~arrayfun(@(id) isvarname(char(id)), ids))
    error('mdemat:tools:visual:InvalidDatasetIds', 'Dataset IDs must be unique.');
end
if any(~ismember(options.required, ids))
    error('mdemat:tools:visual:MissingSamples', 'A required dataset is not in the requested bindings.');
end
group = struct();
missing = struct('DatasetId', {}, 'Reason', {});
for id = ids
    name = char(id);
    reason = "";
    if ~isfield(datasets, name)
        reason = "dataset";
    elseif ~isfield(datasets.(name).Data, field)
        reason = "field";
    else
        entry = datasets.(name);
        view = entry.Data.(field);
        row = find(view.ComponentIds == componentId, 1);
        if isempty(row)
            reason = "component";
        else
            view.Values = view.Values(row, :);
            view.ComponentIds = componentId;
            entry.Data = struct();
            entry.Data.(field) = view;
            group.(name) = entry;
        end
    end
    if strlength(reason) > 0
        missing(end + 1) = struct('DatasetId', id, 'Reason', reason); %#ok<AGROW>
    end
end
if ~isempty(missing) && any(ismember([missing.DatasetId], options.required))
    error('mdemat:tools:visual:MissingSamples', 'A required dataset, field or component is missing.');
end
end
