function A = unfold(T, n, d)
% UNFOLD Tensor matricization.
%
%   A = unfold(T, n, d) returns a matrix @a A obtained by unfolding a
%   tensor @a T of shape @a n along the dimension @a d.
%
%   This function performs tensor matricization by reshaping tensor @a T
%   into matrix form. Dimension @a d becomes the rows, while all other
%   dimensions are flattened into columns. The resulting matrix has
%   dimensions [n(d), prod(n)/n(d)].
%
% See also:
%   core.linalg.fold

p = [d, setdiff(1:numel(n), d)];
T = permute(T, p);
A = reshape(T, n(d), []);
end