function X = zeroing(X, tol)
% ZEROING Set array entries below threshold to zero.
%
%   X = zeroing(X) replaces values in array @a X smaller than the default
%   tolerance with zeros.
%
%   X = zeroing(X, tol) replaces values in array @a X smaller than the
%   specified tolerance @a tol with zeros.
%
% Inputs:
%   X   - Input array
%   tol - Tolerance threshold (optional, default: eps*1e3)
%
% Outputs:
%   X - Output array with small values replaced by zeros
%
% Examples:
%   % Basic usage with default tolerance
%   A = [1, 1e-16, -1e-15, 2];
%   cleanA = zeroing(A); % Returns [1, 0, 0, 2]
%   
%   % With custom tolerance
%   B = [10, 0.01, 0.001, 0.0001, 5];
%   cleanB = zeroing(B, 0.001); % Returns [10, 0.01, 0.001, 0, 5]

if nargin < 2 || isempty(tol), tol = eps * 1e3; end
X(abs(X) < tol) = 0;
end