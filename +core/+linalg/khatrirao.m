function X = khatrirao(U, varargin)
% KHATRIRAO Compute Khatri-Rao product of input matrices.
%
%   X = khatrirao(A, B) computes column-wise Kronecker product of matrices
%   @a A and @a B. For matrices A (I-by-K) and B (J-by-K), result is
%   (I*J)-by-K matrix where each column k is Kronecker product of k-th
%   columns of A and B.
%
%   X = khatrirao(U) computes Khatri-Rao product when @a U is a cell
%   array of matrices.
%
%   The Khatri-Rao product computes column-wise Kronecker products,
%   preserving the number of columns while expanding rows. All input
%   matrices must have the same number of columns. Operation proceeds
%   from right to left through the input matrices.
%
% See also:
%   core.linalg.kronecker

arguments
    U
end
arguments (Repeating)
    varargin
end

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