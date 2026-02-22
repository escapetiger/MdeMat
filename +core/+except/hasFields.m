function [tf, missingField] = hasFields(s, f)
% HASFIELDS Check if a struct @a s has all specified fields @a f and return
% first missing one.
%
%   tf = hasFields(s, f) returns only the logical result without missing
%   field.
%
%   [tf, missingField] = hasFields(s, f) returns true if struct @a s
%   contains all fields specified in @a f, and also returns the first
%   missing field name.

core.except.assert(isstruct(s), 'InvalidInput', ...
    'Input s must be a structure');

isValidF = ischar(f) || isstring(f) || iscell(f) || isempty(f);
core.except.assert(isValidF, 'InvalidInput', ...
    'Input f must be char, string, or cell array');

if isempty(f)
    tf = true;
    if nargout > 1
        missingField = '';
    end
    return;
end

if ischar(f) || isstring(f)
    f = {char(f)};
end

if iscell(f)
    validCellContents = all(cellfun(@(x) ischar(x) || isstring(x), f));
    core.except.assert(validCellContents, 'InvalidInput', ...
        'All elements in cell array f must be char or string');
end

if nargout > 1
    missingField = '';
    for iField = 1:length(f)
        if ~isfield(s, f{iField})
            missingField = f{iField};
            break;
        end
    end
    tf = isempty(missingField);
else
    tf = all(isfield(s, f));
end

end