function T = fold(A, n, d)
% FOLD Inverse of unfold.
%
%   T = fold(A, n, d) returns a tensor @a T of shape @a n obtained by
%   folding a matrix @a A along the dimension @a d.
%
%   This function reverses the unfold operation by reshaping matrix @a A
%   back into tensor form with shape @a n. The folding occurs along
%   dimension @a d, which determines the reshaping pattern. Matrix @a A
%   should have dimensions [n(d), prod(n)/n(d)].
%
% See also:
%   core.linalg.unfold

p = [d, setdiff(1:numel(n), d)];
T = reshape(A, n(p));
q = zeros(1, numel(n));
q(p) = 1:numel(n);
T = permute(T, q);
end