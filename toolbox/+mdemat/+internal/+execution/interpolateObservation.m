function datasets = interpolateObservation(previous, current, fraction)
%INTERPOLATEOBSERVATION Interpolate paired numerical samples, preserving IDs.
%   Exact and fixed entries are supplied at the requested time by the owner.
%   Only numerical values are interpolated. Coordinates/selection/topology
%   must match; a mesh or component change cannot reuse the previous frame.
arguments
    previous (1, 1) struct
    current (1, 1) struct
    fraction (1, 1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
end
if fraction > 1
    error('mdemat:execution:interpolateObservation:Time', 'Interpolation is bounded by paired times.');
end
datasets = current;
names = fieldnames(current);
for index = 1:numel(names)
    name = names{index};
    if ~strcmp(current.(name).Type, 'numeric'), continue; end
    if ~isfield(previous, name)
        error('mdemat:execution:interpolateObservation:Identity', 'A numerical dataset is missing.');
    end
    fields = fieldnames(current.(name).Data);
    for field = 1:numel(fields)
        key = fields{field};
        a = previous.(name).Data.(key); b = current.(name).Data.(key);
        if ~isequal(size(a.Values), size(b.Values)) || ...
                ~isequaln(rmfield(a, 'Values'), rmfield(b, 'Values'))
            error('mdemat:execution:interpolateObservation:Identity', ...
                'Paired field geometry and component IDs must agree.');
        end
        datasets.(name).Data.(key).Values = (1-fraction)*a.Values + fraction*b.Values;
    end
end
end
