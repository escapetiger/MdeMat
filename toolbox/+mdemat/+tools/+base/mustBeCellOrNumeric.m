function mustBeCellOrNumeric(x)
%MUSTBECELLORNUMERIC Validate the outer container as cell or numeric.
%
%   mustBeCellOrNumeric(x) throws an error if x is neither a cell nor a
%   numeric array. Cell contents, shape and finiteness are not checked.
%
% See also:
%   mdemat.tools.base.mustBeCellOrObject

if ~iscell(x) && ~isnumeric(x)
    error('mdemat:tools:base:mustBeCellOrNumeric:InvalidInput', ...
        'Value must be a cell or numeric array.');
end

end
