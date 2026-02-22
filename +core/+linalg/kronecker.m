function K = kronecker(varargin)
% KRONECKER Kronecker product of multiple matrices.
%
%   K = kronecker(A1, A2, ...) computes Kronecker product of all input
%   matrices, with operation carried out from right to left. Extends
%   MATLAB's built-in kron to handle multiple matrices.
%
%   K = kronecker(As) computes Kronecker product when input is a cell
%   array @a As of matrices.
%
%   This function extends MATLAB's built-in kron function to handle
%   multiple matrices efficiently. The operation is computed from right
%   to left, so kronecker(A, B, C) equals kron(A, kron(B, C)).
%
% See also:
%   core.linalg.khatrirao

if iscell(varargin{1})
    As = varargin{1};
else
    As = varargin;
end
K = As{end};
for i = length(As) - 1:-1:1
    K = kron(As{i}, K);
end
end