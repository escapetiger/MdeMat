function tf = hasFields(s, f)
% HASFIELDS Check if a struct @a s has all specified fields @a f.
%
%   tf = hasFields(s, f) returns true if struct @a s contains all fields
%   specified in @a f.
%
% Inputs:
%   s - Structure to check
%   f - Field name (string/char) or cell array of field names
%
% Outputs:
%   tf - True if all fields exist, false otherwise
%
% Examples:
%   % Single field check
%   s = struct('name', 'John', 'age', 30);
%   tf = hasFields(s, 'name'); % Returns true
%
%   % Multiple fields check
%   tf = hasFields(s, {'name', 'age'}); % Returns true
%
%   % Empty field list
%   tf = hasFields(s, {}); % Returns true

if ~isstruct(s)
    tf = false;
    return;
end

if isempty(f)
    tf = true;
    return;
end

if ischar(f) || isstring(f)
    f = {f};
end

if ~iscell(f)
    tf = false;
    return;
end

tf = all(isfield(s, f));
end