function A = sparseFromTriplets(T, m, n)
% SPARSEFROMTRIPLETS Create sparse matrix from triplet format.
%
%   A = sparseFromTriplets(T, m, n) converts a triplet list
%   (row, column, value) into a sparse matrix, automatically
%   filtering out invalid indices (zeros or negative values).
%
%   This function creates a sparse matrix from triplet format where @a T
%   contains [row, column, value] triplets. Invalid indices (zero or
%   negative row/column indices) are automatically filtered out before
%   creating the sparse matrix of size @a m by @a n.

k = find((T(:, 1) > 0) & (T(:, 2) > 0) & (abs(T(:, 3)) > 1e-14));
A = sparse(T(k, 1), T(k, 2), T(k, 3), m, n);
end