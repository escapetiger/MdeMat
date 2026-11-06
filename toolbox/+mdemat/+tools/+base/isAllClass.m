function tf = isAllClass(values, classes)
%ISALLCLASS Test every array element against an allowed class list using isa.
%   Cell contents are checked rather than their containers. Allowed classes
%   are text scalars; subclasses match. Empty arrays satisfy the predicate.
%   Heterogeneous object arrays retain per-element class semantics.
if ischar(classes)
    classes = {classes};
elseif isstring(classes)
    if any(ismissing(classes), 'all')
        error('mdemat:tools:base:isAllClass:InvalidInput', 'Class names cannot be missing.');
    end
    classes = cellstr(classes(:));
end
if ~iscell(classes) || ~all(cellfun(@(name) ...
        (ischar(name) && isrow(name) && ~isempty(name)) || ...
        (isstring(name) && isscalar(name) && ~ismissing(name) && strlength(name) > 0), ...
        classes), 'all')
    error('mdemat:tools:base:isAllClass:InvalidInput', ...
        'Allowed classes must be text scalars or a list of text scalars.');
end
classes = cellfun(@char, classes(:), 'UniformOutput', false);
tf = true;
if isempty(values), return; end
if ~iscell(values) && ~isobject(values)
    tf = any(cellfun(@(name) isa(values, name), classes));
    return;
end
for index = 1:numel(values)
    if iscell(values), value = values{index}; else, value = values(index); end
    if ~any(cellfun(@(name) isa(value, name), classes))
        tf = false;
        return;
    end
end
end
