function mustBePenaltyType(x)
% MUSTBEPENALTYTYPE Validate penalty type input.

if ~isequal(size(x), [2, 2])
    error('Penalty types are incomplete.');
end

y = {{'boundary', 'boundaryw', 'all', 'allw'}, {'central', 'right', 'left'}; ...
    {'central', 'right', 'left'}, {}};

for i = 1:numel(x)
    if ~(ismember(x{i}, y{i}) || isempty(x{i}))
        error('Penalty type %s is not supported.', x{i});
    end
end
end