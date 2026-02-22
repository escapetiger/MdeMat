function s = toString(x)
% TOSTRING Convert a variable to a formatted string.
%
%   s = toString(x) returns a formatted string representing runtime
%   variable @a x.
%
% Inputs:
%   x - Variable
%
% Outputs:
%   s - String representation
%
% Examples:
%   % Basic usage
%   s = core.data.toString([]); % Returns '[]'
%   s = core.data.toString(1); % Returns '1'
%   s = core.data.toString(true); % Returns 'true'
%   s = core.data.toString('F'); % Returns '''F'''
%   s = core.data.toString({1, 2, 3}); % Returns 'cell'
%   s = core.data.toString(struct('a', 1)); % Returns 'struct'
%   s = core.data.toString(@(z) sin(z)); % Returns '@(z)sin(z)'

if isempty(x)
    s = '[]';
elseif isnumeric(x) && isscalar(x)
    s = num2str(x);
elseif islogical(x) && isscalar(x)
    if x
        s = 'true';
    else
        s = 'false';
    end
elseif ischar(x)
    s = ['''', x, ''''];
elseif iscell(x) && length(x) <= 3
    s = 'cell';
elseif isstruct(x)
    s = 'struct';
elseif isa(x, 'function_handle')
    s = func2str(x);
else
    s = class(x);
end

end