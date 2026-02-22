function Y = onesLike(X, d)
% ONESLIKE Create array of ones with size similar to input..
%
%   Y = onesLike(X) returns an array of the same size as @a X filled with
%   ones.
%
%   Y = onesLike(X, d) returns an array matching size along the dimension
%   @a d, which is filled with ones.
%
% Inputs:
%   X - Input array
%   d - Dimension index (optional)
%
% Outputs:
%   Y - Output array filled with ones
%
% Examples:
%   % Create ones array matching 3D array
%   A = rand(3, 4, 5);
%   Y1 = onesLike(A); % Returns 3x4x5 array of ones
%   
%   % Create ones array matching one dimension
%   Y2 = onesLike(A, 2); % Returns 1x4 array of ones
%
% See Also:
%   core.linalg.zerosLike

if nargin < 2 || isempty(d)
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