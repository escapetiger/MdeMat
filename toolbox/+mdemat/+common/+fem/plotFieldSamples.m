function fields = plotFieldSamples(space, dofs, components, plan)
%PLOTFIELDSAMPLES Convert explicitly selected FE coefficients to value views.
%   Components contains original column IDs. Each returned field has Values
%   [components,points] and matching ComponentIds; tensor point permutations
%   are applied to values and geometry together.
arguments
    space (1, 1) mdemat.common.fem.FiniteElementSpace
    dofs (1, 1) struct
    components (1, 1) struct
    plan (1, 1) struct
end
fields = struct();
names = fieldnames(components);
for index = 1:numel(names)
    name = names{index};
    ids = components.(name);
    if isempty(ids) || ~isfield(dofs, name) || isempty(dofs.(name)), continue; end
    mustBeNumeric(ids); mustBeVector(ids); mustBePositive(ids);
    mustBeInteger(ids); mustBeFinite(ids);
    ids = double(ids(:).');
    if numel(unique(ids)) ~= numel(ids) || any(ids > size(dofs.(name), 2))
        error('mdemat:fem:plotFieldSamples:Components', ...
            'Selected IDs must be unique columns present in field %s.', name);
    end
    values = space.eval(plan.Nodes, dofs.(name)(:, ids)).';
    view = struct('Values', values, 'ComponentIds', ids);
    if isfield(plan, 'Axes')
        view.Axes = plan.Axes;
        view.Values = values(:, plan.Permutation);
    else
        view.Coordinates = plan.Coordinates;
    end
    fields.(name) = view;
end
end
