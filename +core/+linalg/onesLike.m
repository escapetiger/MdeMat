function Y = onesLike(X, d)
% ONESLIKE Create array of ones with size similar to input.
%
%   Y = onesLike(X) returns an array of the same size as @a X filled with
%   ones.
%
%   Y = onesLike(X, d) returns an array matching size along the dimension
%   @a d, which is filled with ones.
%
%   This function creates an array filled with ones that matches the
%   dimensions of @a X or a specific dimension @a d. When @a d exceeds
%   the number of dimensions in @a X, returns a 1x1 array of ones.
%
% See also:
%   core.linalg.zerosLike

arguments
    X
    d {mustBeInteger, mustBePositive} = []
end

if isempty(d)
    n = size(X);
else
    if d > ndims(X)
        n = [1, 1];
    else
        n = [1, size(X, d)];
    end
end

Y = ones(n);
end