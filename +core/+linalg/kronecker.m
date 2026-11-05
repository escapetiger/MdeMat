function K = kronecker(varargin)
% KRONECKER Kronecker product of multiple matrices.
%
%   K = kronecker(A1, A2, ...) computes Kronecker product of all input
%   matrices, with operation carried out from right to left. Extends
%   MATLAB's built-in kron to handle multiple As.
%
% Inputs:
%   varargin - Input matrices
%
% Outputs:
%   K - Kronecker product
%
% Examples:
%   % Basic usage
%   A = [1 2; 3 4];
%   B = [5 6; 7 8];
%   K = kronecker(A, B); % K will be a 4x4 matrix
%
%   % Using cell array
%   K = kronecker({A, B, eye(2)}); % K will be an 8x8 matrix
%
% See Also:
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