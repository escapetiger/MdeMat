function mustBeCellOrObject(x)
%MUSTBECELLOROBJECT Validate the outer container as cell or object.
%
%   mustBeCellOrObject(x) throws an error if x is neither a cell nor
%   an object array. Cell contents and element classes are not checked.
%
% See also:
%   mdemat.tools.base.mustBeCellOrNumeric

if ~iscell(x) && ~isobject(x)
    error('mdemat:tools:base:mustBeCellOrObject:InvalidInput', ...
        'Value must be a cell or object array.');
end

end
