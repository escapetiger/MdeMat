function [A, b] = create_amg_test_matrix(n)
% CREATE_AMG_TEST_MATRIX Creates a 2D Poisson matrix ideal for AMG testing.
%
%   [A, b] = create_amg_test_matrix(n) creates an (n-2)^2 x (n-2)^2 sparse
%   matrix A and right-hand side b from discretizing the 2D Poisson equation:
%   
%   -∇²u = f  on Ω = (0,1) x (0,1)
%    u = 0    on ∂Ω
%
%   This produces a symmetric positive definite matrix with excellent 
%   AMG convergence properties.

if nargin < 1
    n = 50;  % Default grid size
end

% Interior grid points
ni = n - 2;
h = 1 / (n - 1);

% Create 2D Laplacian using finite differences
% -u_{i-1,j} - u_{i+1,j} - u_{i,j-1} - u_{i,j+1} + 4*u_{i,j} = h²*f_{i,j}

% Number interior points
N = ni * ni;

% Build matrix using 5-point stencil
row = zeros(5*N, 1);
col = zeros(5*N, 1);
val = zeros(5*N, 1);
idx = 0;

for j = 1:ni
    for i = 1:ni
        % Global index for point (i,j)
        k = (j-1)*ni + i;
        
        % Center point: 4/h²
        idx = idx + 1;
        row(idx) = k;
        col(idx) = k;
        val(idx) = 4 / h^2;
        
        % Left neighbor: -1/h²
        if i > 1
            idx = idx + 1;
            row(idx) = k;
            col(idx) = k - 1;
            val(idx) = -1 / h^2;
        end
        
        % Right neighbor: -1/h²
        if i < ni
            idx = idx + 1;
            row(idx) = k;
            col(idx) = k + 1;
            val(idx) = -1 / h^2;
        end
        
        % Bottom neighbor: -1/h²
        if j > 1
            idx = idx + 1;
            row(idx) = k;
            col(idx) = k - ni;
            val(idx) = -1 / h^2;
        end
        
        % Top neighbor: -1/h²
        if j < ni
            idx = idx + 1;
            row(idx) = k;
            col(idx) = k + ni;
            val(idx) = -1 / h^2;
        end
    end
end

% Trim arrays
row = row(1:idx);
col = col(1:idx);
val = val(1:idx);

% Create sparse matrix
A = sparse(row, col, val, N, N);

% Create right-hand side: f(x,y) = 2π²sin(πx)sin(πy)
% This gives exact solution u(x,y) = sin(πx)sin(πy)
b = zeros(N, 1);
for j = 1:ni
    for i = 1:ni
        x = i * h;
        y = j * h;
        k = (j-1)*ni + i;
        b(k) = 2 * pi^2 * sin(pi*x) * sin(pi*y) * h^2;
    end
end

fprintf('Created %dx%d Poisson matrix:\n', N, N);
fprintf('  Grid: %dx%d interior points\n', ni, ni);
fprintf('  h = %.6f\n', h);
fprintf('  nnz = %d\n', nnz(A));
fprintf('  Condition number ≈ O(h^-2) = O(%.0f)\n', 1/h^2);
fprintf('  Matrix properties: SPD, M-matrix, 5-point stencil\n');
fprintf('  Expected AMG convergence: Excellent\n');

end