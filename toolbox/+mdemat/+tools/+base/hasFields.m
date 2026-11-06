function [tf, missingField] = hasFields(value, names)
%HASFIELDS Check every requested field and return the first missing name.
%   Names are a text scalar or an array of text scalars. The result is always
%   logical scalar, independent of output count. An empty list is true;
%   empty names and nonscalar text entries are rejected. Missing names are
%   visited in column-major order. Struct arrays have a common field set.
arguments
    value {mustBeA(value, 'struct')}
    names
end
if isempty(names) && ~ischar(names)
    tf = true;
    missingField = '';
    return;
end
if ischar(names)
    names = {names};
elseif isstring(names)
    if any(ismissing(names), 'all')
        error('mdemat:tools:base:hasFields:InvalidInput', 'Field names cannot be missing.');
    end
    names = cellstr(names(:));
end
if ~iscell(names) || ~all(cellfun(@validName, names), 'all')
    error('mdemat:tools:base:hasFields:InvalidInput', ...
        'Field names must be nonempty text scalars.');
end
names = cellfun(@char, names(:), 'UniformOutput', false);
missing = ~isfield(value, names);
tf = ~any(missing);
index = find(missing, 1);
missingField = '';
if ~isempty(index), missingField = names{index}; end
end

function tf = validName(name)
tf = (ischar(name) && isrow(name) && ~isempty(name)) || ...
    (isstring(name) && isscalar(name) && ~ismissing(name) && strlength(name) > 0);
end
