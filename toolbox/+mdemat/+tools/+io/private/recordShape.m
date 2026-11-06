function [shape, count] = recordShape(shape, limits, sparseArray)
if nargin < 3, sparseArray = false; end
if ~isnumeric(shape) || ~isreal(shape) || ~isvector(shape) || ...
        numel(shape) < 2 || ...
        any(~isfinite(shape) | shape < 0 | shape ~= fix(shape) | shape > flintmax)
    error('mdemat:io:InvalidRecord','Invalid record dimensions.');
end
shape = double(shape(:).');
if numel(shape) > limits.MaxDimensions
    error('mdemat:io:ResourceLimit','Record dimensions exceed the dimensionality budget.');
end
count = prod(shape);
if ~sparseArray && (any(shape > limits.MaxElements) || count > limits.MaxElements)
    error('mdemat:io:ResourceLimit','Record dimensions exceed the element budget.');
end
end
