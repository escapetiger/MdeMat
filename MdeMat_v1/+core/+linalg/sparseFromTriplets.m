function A = sparseFromTriplets(T, m, n)
% SPARSEFROMTRIPLETS Create sparse matrix from triplet format.
%
%   A = sparseFromTriplets(T, m, n) converts a triplet list
%   (row, column, value) into a sparse matrix, automatically
%   filtering out invalid indices (zeros or negative values).
%
% Inputs:
%   T - Triplet matrix with columns [row, col, val]
%   m - Number of rows (positive integer)
%   n - Number of columns (positive integer)
%
% Outputs:
%   A - Sparse matrix of size m×n

k = find((T(:, 1) > 0) & (T(:, 2) > 0));
A = sparse(T(k, 1), T(k, 2), T(k, 3), m, n);
end