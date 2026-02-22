function tf = isAllClass(x, cls)
% ISALLCLASS Check if all elements of @a x belong to @a cls.
%
%   tf = isAllClass(x, cls) returns true if all elements in @a x belong to
%   one of the classes specified in @a cls.
%
% Inputs:
%   x - Element(s) (object array or cell array)
%   cls - Class name(s) (string or cell array of strings)
%
% Outputs:
%   tf - True only if all elements belong to the specified class(es)
%
% Examples:
%   % Check if all elements are numeric
%   tf = isAllClass([1, 2, 3], 'double');
%
%   % Check multiple allowed classes
%   mixedData = {1, 2.5, int32(3)};
%   tf = isAllClass(mixedData, {'double', 'int32'});
%
%   % Empty array returns true
%   tf = isAllClass([], 'any_class');
%
% See Also:
%   core.validate.isAllSameClass

if isempty(x)
    tf = true;
    return;
end

if ~iscell(cls), cls = {cls}; end

tf = true;
for j = 1:length(x)
    tf0 = false;
    for i = 1:length(cls)
        if iscell(x)
            y = x{j};
        else
            y = x(j);
        end
        tf0 = tf0 | isa(y, cls{i});
        if tf0, break; end
    end
    tf = tf & tf0;
    if ~tf, return; end
end

end