function tables = convergenceTables(summary)
%CONVERGENCETABLES Format series without losing component/reference identity.
%   UserData retains the full identity and resolved IDs at every level. Table
%   keys are readable labels, made unique and bounded by namelengthmax; they
%   are never decoded to recover an identity. NComponents tracks mode cutoff.
tables = struct();
for index = 1:numel(summary)
    series = summary(index);
    identity = rmfield(series.Identity, 'SpatialNorm');
    name = findIdentity(tables, identity);
    if isempty(name)
        name = uniqueName(proposedName(identity), fieldnames(tables));
        resolved = series.ResolvedComponentIds;
        tables.(name) = table(series.Scales(:), cellfun(@numel, resolved(:)), ...
            'VariableNames', {'Scale', 'NComponents'});
        tables.(name).Properties.UserData = struct('Identity', identity, ...
            'ResolvedComponentIds', {resolved(:)});
    elseif ~isequal(tables.(name).Scale, series.Scales(:)) || ...
            ~isequal(tables.(name).Properties.UserData.ResolvedComponentIds, ...
            series.ResolvedComponentIds(:))
        error('mdemat:execution:convergenceTables:SeriesMismatch', ...
            'Norm columns of one identity must use matching scales and component selections.');
    end
    normName = char(series.Identity.SpatialNorm);
    tables.(name).([normName, '_Error']) = series.Errors(:);
    rates = NaN(size(series.Errors(:)));
    rates(2:end) = series.LocalRates(:);
    tables.(name).([normName, '_Order']) = rates;
end
end

function name = findIdentity(tables, identity)
name = '';
keys = fieldnames(tables);
for index = 1:numel(keys)
    if isequaln(tables.(keys{index}).Properties.UserData.Identity, identity)
        name = keys{index};
        return;
    end
end
end

function name = proposedName(identity)
field = char(identity.Field);
ids = identity.ComponentIds;
if strcmp(identity.ComponentNorm, 'sequence')
    if ischar(ids)
        suffix = '_all';
    else
        suffix = sprintf('_%d', ids);
    end
    name = [field, suffix, '_sequence'];
else
    name = sprintf('%s%d', field, ids);
end
if numel(name) > namelengthmax
    name = [field, '_sequence'];
end
name = matlab.lang.makeValidName(name);
end

function name = uniqueName(base, existing)
maximum = namelengthmax;
name = base(1:min(numel(base), maximum));
index = 1;
while ismember(name, existing)
    index = index + 1;
    suffix = sprintf('_%d', index);
    name = [base(1:min(numel(base), maximum-numel(suffix))), suffix];
end
end
