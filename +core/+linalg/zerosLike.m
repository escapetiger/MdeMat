function Y = zerosLike(X, d)
% ZEROSLIKE Create array of zeros with size similar to input.
%
%   Y = zerosLike(X) returns an array of the same size as @a X filled with
%   zeros.
%
%   Y = zerosLike(X, d) returns an array matching size along the dimension
%   @a d, which is filled with zeros.
%
%   This function creates an array filled with zeros that matches the
%   dimensions of @a X or a specific dimension @a d. When @a d exceeds
%   the number of dimensions in @a X, returns a 1x1 array of zeros.
%
% See also:
%   core.linalg.onesLike

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

Y = zeros(n);
end