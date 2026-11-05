function T = fold(A, n, d)
% FOLD Inverse of unfold.
%
%   T = fold(A, n, d) returns a tensor @a T of shape @a n obtained by
%   folding a matrix @a A along the dimension @a d.
%
% Inputs:
%   A - Matrix of shape [n(d), prod(n)/n(d)]
%   n - Tensor shape
%   d - Mode along which to fold
%
% Outputs:
%   T - Folded tensor of shape @a n
%
% Examples:
%   % Unfold T along dimension 2, then fold it back
%   T = rand(3, 4, 5);
%   n = size(T);
%   A = core.linalg.unfold(T, n, 2);
%   TT = core.linalg.fold(A, n, 2); % TT equals to T
%
% See Also:
%   core.linalg.unfold

p = [d, setdiff(1:numel(n), d)];
T = reshape(A, n(p));
q = zeros(1, numel(n));
q(p) = 1:numel(n);
T = permute(T, q);
end