function A = unfold(T, n, d)
% UNFOLD Tensor matricization.
%
%   A = unfold(T, n, d) returns a matrix @a A obtained by unfolding a
%   tensor @a T of shape @a n along the dimension @a d.
%
% Inputs:
%   T - Tensor of shape @a n
%   n - Tensor shape
%   d - Mode along which to unfold
%
% Outputs:
%   A - Unfolded matrix of shape [n(d), prod(n)/n(d)]
%
% Examples:
%   % Unfold 3D tensor along dimension 2
%   T = rand(3, 4, 5);
%   A = unfold(T, size(T), 2); % A will be a 4x15 matrix
%
% See Also:
%   core.linalg.fold

p = [d, setdiff(1:numel(n), d)];
T = permute(T, p);
A = reshape(T, n(d), []);
end