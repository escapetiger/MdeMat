function datasets = sampleObservation(scheme, state, time, includeNumeric)
%SAMPLEOBSERVATION Sample an explicitly requested plotting frame in execution.
%   Geometry/FEM/physical reference interpretation belongs here and in the
%   model's samplePlotFields/samplePlotReferences methods, never in tools.
arguments
    scheme mdemat.internal.execution.Scheme
    state mdemat.common.fem.SpatialState
    time = []
    includeNumeric (1, 1) logical = true
end
options = scheme.Observation.Plot;
datasets = struct();
if isempty(options.PointsPerDimension) || any(options.PointsPerDimension <= 0)
    return;
end
plan = mdemat.common.fem.plotSamplingPlan(state.XDisc, options.PointsPerDimension);
types = options.DatasetTypes;
if isempty(fieldnames(types)), types.NUMERIC = 'numeric'; end
names = fieldnames(types);
numeric = []; exact = [];
for index = 1:numel(names)
    name = names{index};
    type = char(types.(name));
    switch type
        case 'numeric'
            if ~includeNumeric, continue; end
            if isempty(numeric), numeric = scheme.samplePlotFields(state, plan); end
            fields = numeric;
        case 'exact'
            if isempty(exact), exact = scheme.samplePlotReferences(state, plan, time); end
            fields = exact;
        case 'fixed'
            if isfield(options.FixedData, name)
                datasets.(name) = options.FixedData.(name);
            end
            continue;
        otherwise
            error('mdemat:execution:sampleObservation:DatasetType', ...
                'Dataset %s has an unsupported source type.', name);
    end
    datasets.(name) = struct('Type', type, 'Data', fields);
end
end
