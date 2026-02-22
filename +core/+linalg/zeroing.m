function X = zeroing(X, tol)
% ZEROING Set array entries below threshold to zero.
%
%   X = zeroing(X) replaces values in array @a X smaller than the default
%   tolerance with zeros.
%
%   X = zeroing(X, tol) replaces values in array @a X smaller than the
%   specified tolerance @a tol with zeros.
%
%   This function removes numerical noise by setting values with absolute
%   magnitude below tolerance @a tol to zero. Uses default tolerance of
%   eps*1e3 when @a tol is not specified. Operation is applied
%   element-wise to the entire array @a X.

if nargin < 2 || isempty(tol), tol = eps * 1e3; end
X(abs(X) < tol) = 0;
end