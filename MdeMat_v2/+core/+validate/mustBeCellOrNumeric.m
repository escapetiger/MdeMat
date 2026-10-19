function mustBeCellOrNumeric(x)
% MUSTBECELLORNUMERIC Validate that @a x is a cell or a number.
%
%   mustBeCellOrNumeric(x) throws an error if @a x is neither a cell nor a
%   number.
%
% Inputs:
%   x - Value to validate
%
% Outputs:
%   NULL
%
% Examples:
%   % Valid numeric input
%   mustBeCellOrNumeric([1, 2, 3]); % Valid
%   
%   % Valid cell input
%   mustBeCellOrNumeric({1, 'test', []}); % Valid
%   
%   % Invalid input throws error
%   % mustBeCellOrNumeric('abc'); % Error
%
% See Also:
%   core.validate.mustBeCellOrObject

core.except.assert(isempty(x) || (iscell(x) && isnumeric(x)), ...
    'InxidInput', 'Value must be either a cell or numeric type.')

end