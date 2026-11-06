function fields = plotReferenceFields(space, exacts, components, plan, time)
%PLOTREFERENCEFIELDS Adapt spatial references to the same explicit value views.
arguments
    space (1, 1) mdemat.common.fem.FiniteElementSpace
    exacts (1, 1) struct
    components (1, 1) struct
    plan (1, 1) struct
    time = []
end
fields = struct();
names = fieldnames(components);
for index = 1:numel(names)
    name = names{index}; ids = components.(name);
    if isempty(ids) || ~isfield(exacts, name) || isempty(exacts.(name)), continue; end
    values = mdemat.common.fem.referenceValues(space, exacts.(name), plan.Nodes, time);
    ids = double(ids(:).');
    if any(ids > size(values, 1))
        error('mdemat:execution:plotReferenceFields:Components', ...
            'Reference field %s does not supply the selected components.', name);
    end
    view = struct('Values', values(ids, :), 'ComponentIds', ids);
    if isfield(plan, 'Axes')
        view.Axes = plan.Axes;
        view.Values = view.Values(:, plan.Permutation);
    else
        view.Coordinates = plan.Coordinates;
    end
    fields.(name) = view;
end
end
