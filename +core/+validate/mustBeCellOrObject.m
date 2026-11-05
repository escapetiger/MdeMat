function mustBeCellOrObject(x)
% MUSTBECELLOROBJ Validate that @a x is a cell or an object.
%
%   mustBeCellOrObject(x) throws an error if @a x is neither a cell nor
%   an object.
%
% Inputs:
%   x - Value to validate
%
% Outputs:
%   NULL
%
% Examples:
%   % Valid object input
%   obj = struct('field', 1);
%   mustBeCellOrObject(obj); % Valid
%   
%   % Valid cell input
%   mustBeCellOrObject({1, 2, 3}); % Valid
%   
%   % Invalid input throws error
%   % mustBeCellOrObject(123); % Error
%
% See Also:
%   core.validate.mustBeCellOrNumeric

core.except.assert(isempty(x) || (iscell(x) && isobject(x)), ...
    'InxidInput', 'Value must be either a cell or object type.')

end