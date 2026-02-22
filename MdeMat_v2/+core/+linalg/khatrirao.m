function X = khatrirao(U, varargin)
% KHATRIRAO Compute Khatri-Rao product of input matrices.
%
%   X = khatrirao(A, B) computes column-wise Kronecker product of matrices
%   @a A and @a B. For matrices A (I-by-K) and B (J-by-K), result is
%   (I*J)-by-K matrix where each column k is Kronecker product of k-th
%   columns of A and B.
%
% Inputs:
%   U - Cell array of matrices or head input matrix
%   varargin - Tailed input matrices
%
% Outputs:
%   X - Khatri-Rao product matrix
%
% Examples:
%   % Two-matrix Khatri-Rao product
%   A = [1 2; 3 4];    % 2-by-2 matrix
%   B = [5 6; 7 8];    % 2-by-2 matrix
%   X = khatrirao(A, B);
%
%   % Sparse matrix Khatri-Rao product
%   As = sparse([1 0; 3 4]);
%   Bs = sparse([5 0; 0 8]);
%   X = khatrirao(As, Bs);
%
% See Also:
%   core.linalg.kronecker

if ~iscell(U), U = [{U}, varargin]; end
X = U{end};
J = size(X, 1);
for n = length(U) - 1:-1:1
    A = U{n};
    I = size(A, 1);
    X = arrayfun(@(k) kron(A(:,k), X(:,k)), 1:size(U{end}, 2), 'Un', 0);
    X = horzcat(X{:});
    J = I * J;
end
end