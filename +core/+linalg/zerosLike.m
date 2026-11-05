function Y = zerosLike(X, d)
% ZEROSLIKE Create array of zeros with size similar to input..
%
%   Y = zerosLike(X) returns an array of the same size as @a X filled with
%   zeros.
%
%   Y = zerosLike(X, d) returns an array matching size along the dimension
%   @a d, which is filled with zeros.
%
% Inputs:
%   X - Input array
%   d - Dimension index (optional)
%
% Outputs:
%   Y - Output array filled with zeros
%
% Examples:
%   % Create zeros array matching 3D array
%   A = rand(3, 4, 5);
%   Y1 = zerosLike(A); % Returns 3x4x5 array of zeros
%   
%   % Create zeros array matching one dimension
%   Y2 = zerosLike(A, 2); % Returns 1x4 array of zeros
%
% See Also:
%   core.linalg.onesLike

if nargin < 2 || isempty(d)
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