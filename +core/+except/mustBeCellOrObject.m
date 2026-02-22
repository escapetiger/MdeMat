function mustBeCellOrObject(x)
% MUSTBECELLOROBJ Validate that @a x is a cell or an object.
%
%   mustBeCellOrObject(x) throws an error if @a x is neither a cell nor
%   an object.
%
% See also:
%   core.except.mustBeCellOrNumeric

core.except.assert(iscell(x) || isobject(x), ...
    'InvalidInput', 'Value must be either a cell or object type.')

end