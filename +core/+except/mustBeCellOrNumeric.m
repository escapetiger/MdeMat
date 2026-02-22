function mustBeCellOrNumeric(x)
% MUSTBECELLORNUMERIC Validate that @a x is a cell or a number.
%
%   mustBeCellOrNumeric(x) throws an error if @a x is neither a cell nor a
%   number.
%
% See also:
%   core.except.mustBeCellOrObject

core.except.assert(iscell(x) || isnumeric(x), ...
    'InvalidInput', 'Value must be either a cell or numeric type.')

end