function summaries = summarizeConvergence(observations)
%SUMMARIZECONVERGENCE Group explicit scalar observations without running a study.
%   Each observation has Scale, Field, ComponentIds, SpatialNorm, ComponentNorm,
%   Comparison, Reference and Error. Identity fields must agree within a series;
%   Scale is strictly decreasing in observation order. No missing level is filled.
%   ComponentIds='all' denotes the sequence over every component at each level;
%   it requires ComponentNorm='sequence' and ResolvedComponentIds=1:N for that
%   observation. Numeric selections retain their ordered IDs as their identity.
%   Summaries retain ResolvedComponentIds as a cell per level, alongside
%   Identity, Scales, Errors, LocalRates, Defined and Finite.
%   A single observation is retained with empty rates. No State, sampler, exact
%   callback, mesh, or graphics resource is accepted or inferred here.
%   Unassessed observations use Error=NaN and retain undefined rates.
arguments
    observations struct
end
identityFields = {'Field','ComponentIds','SpatialNorm','ComponentNorm','Comparison','Reference'};
required = [{'Scale','Error'}, identityFields];
summaries = struct('Identity', {}, 'Scales', {}, 'Errors', {}, ...
    'ResolvedComponentIds', {}, 'LocalRates', {}, 'Defined', {}, 'Finite', {});
for index = 1:numel(observations)
    observation = observations(index);
    if ~all(isfield(observation, required)) || ...
            ~isnumeric(observation.Scale) || ~isscalar(observation.Scale) || ...
            ~isreal(observation.Scale) || ~isfinite(observation.Scale) || observation.Scale <= 0 || ...
            ~isnumeric(observation.Error) || ~isscalar(observation.Error) || ...
            ~isreal(observation.Error) || isinf(observation.Error) || observation.Error < 0
        error('mdemat:analysis:InvalidObservation', 'Use finite scalar scales/errors and explicit identity.');
    end
    identity = rmfield(observation, setdiff(fieldnames(observation), identityFields));
    for field = ["Field", "SpatialNorm", "ComponentNorm", "Comparison", "Reference"]
        value = identity.(field);
        if ~((ischar(value) && (isrow(value) || isempty(value))) || ...
                (isstring(value) && isscalar(value) && ~ismissing(value)))
            error('mdemat:analysis:InvalidObservation', 'Identity text must be scalar.');
        end
        identity.(field) = char(value);
    end
    if ~ismember(identity.SpatialNorm, {'L1','L2','LInf','max'}) || ...
            ~ismember(identity.ComponentNorm, {'none','L1','L2','LInf','sequence'}) || ...
            ~ismember(identity.Comparison, {'exact','interlevel','unassessed'})
        error('mdemat:analysis:InvalidObservation', 'Unknown norm or comparison identity.');
    end
    ids = identity.ComponentIds;
    allComponents = ((ischar(ids) && isrow(ids)) || ...
        (isstring(ids) && isscalar(ids) && ~ismissing(ids))) && strcmpi(ids, 'all');
    if allComponents
        if ~strcmp(identity.ComponentNorm, 'sequence') || ...
                ~isfield(observation, 'ResolvedComponentIds')
            error('mdemat:analysis:InvalidObservation', ...
                '''all'' requires a sequence and explicit resolved component IDs.');
        end
        identity.ComponentIds = 'all';
        resolved = componentIds(observation.ResolvedComponentIds);
        if ~isequal(resolved, 1:numel(resolved))
            error('mdemat:analysis:InvalidObservation', ...
                '''all'' must resolve to every component in order, 1:N.');
        end
    else
        identity.ComponentIds = componentIds(ids);
        resolved = identity.ComponentIds;
        if isfield(observation, 'ResolvedComponentIds') && ...
                ~isequal(componentIds(observation.ResolvedComponentIds), resolved)
            error('mdemat:analysis:InvalidObservation', ...
                'Resolved component IDs must match a numeric selection in order.');
        end
    end
    if isempty(identity.Field) || ...
            (strcmp(identity.ComponentNorm,'none') && numel(identity.ComponentIds) ~= 1) || ...
            (strcmp(identity.Comparison,'unassessed') ~= isnan(observation.Error))
        error('mdemat:analysis:InvalidObservation','Observation identity and assessed status disagree.');
    end
    identity = orderfields(identity);
    match = find(arrayfun(@(item) isequal(item.Identity, identity), summaries), 1);
    if isempty(match)
        match = numel(summaries)+1;
        summaries(match) = struct('Identity', identity, 'Scales', [], 'Errors', [], ...
            'ResolvedComponentIds', {cell(0,1)}, 'LocalRates', [], ...
            'Defined', false(0,1), 'Finite', false(0,1));
    end
    summaries(match).Scales(end+1,1) = double(observation.Scale);
    summaries(match).Errors(end+1,1) = double(observation.Error);
    summaries(match).ResolvedComponentIds{end+1,1} = resolved;
end
for index = 1:numel(summaries)
    if numel(summaries(index).Scales) < 2, continue; end
    if strcmp(summaries(index).Identity.Comparison,'unassessed')
        if any(diff(summaries(index).Scales) >= 0)
            error('mdemat:analysis:InvalidRefinement','Observation scales must strictly decrease.');
        end
        count = numel(summaries(index).Scales)-1;
        summaries(index).LocalRates = NaN(count,1);
        summaries(index).Defined = false(count,1);
        summaries(index).Finite = false(count,1);
        continue;
    end
    rates = mdemat.tools.analyze.convergenceRates( ...
        summaries(index).Scales, summaries(index).Errors);
    summaries(index).LocalRates = rates.LocalRates;
    summaries(index).Defined = rates.Defined;
    summaries(index).Finite = rates.Finite;
end
end

function ids = componentIds(ids)
if ~isnumeric(ids) || ~isreal(ids) || ~isvector(ids) || isempty(ids) || ...
        any(~isfinite(ids) | ids < 1 | ids ~= fix(ids)) || numel(unique(ids)) ~= numel(ids)
    error('mdemat:analysis:InvalidObservation', 'Component IDs must be unique positive integers.');
end
ids = double(ids(:).');
end
