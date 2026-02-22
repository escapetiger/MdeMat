function [A, b] = create_nonsymmetric_matrix(n, Pe)
% CREATE_NONSYMMETRIC_MATRIX Creates a 2D advection-diffusion matrix.
%
%   [A, b] = create_nonsymmetric_matrix(n, Pe) creates an (n-2)^2 x (n-2)^2
%   sparse matrix A and right-hand side b from discretizing:
%   
%   -ε∇²u + β·∇u = f  on Ω = (0,1) x (0,1)
%    u = 0            on ∂Ω
%
%   where β = [1, 1] (advection field) and Pe is the Péclet number.
%   High Pe → advection-dominated (non-symmetric, harder for AMG)
%   Low Pe → diffusion-dominated (nearly symmetric, good for AMG)

if nargin < 1
    n = 50;  % Default grid size
end
if nargin < 2
    Pe = 100;  % Default Péclet number (high = advection-dominated)
end

% Interior grid points
ni = n - 2;
h = 1 / (n - 1);

% Péclet number controls advection vs diffusion balance
epsilon = 1 / Pe;  % Diffusion coefficient
beta = [1, 1];     % Advection velocity

% Number interior points
N = ni * ni;

% Build matrix using upwind finite differences for advection
row = zeros(5*N, 1);
col = zeros(5*N, 1);
val = zeros(5*N, 1);
idx = 0;

fprintf('Creating %dx%d advection-diffusion matrix:\n', N, N);
fprintf('  Grid: %dx%d, h=%.6f\n', ni, ni, h);
fprintf('  Péclet number: %.1f\n', Pe);
fprintf('  Regime: %s\n', Pe > 10);

for j = 1:ni
    for i = 1:ni
        % Global index for point (i,j)
        k = (j-1)*ni + i;
        
        % Diffusion part: -ε∇²u (same as Laplacian)
        % Center point: 4ε/h²
        idx = idx + 1;
        row(idx) = k;
        col(idx) = k;
        val(idx) = 4 * epsilon / h^2;
        
        % Diffusion: Left neighbor
        if i > 1
            idx = idx + 1;
            row(idx) = k;
            col(idx) = k - 1;
            val(idx) = -epsilon / h^2;
        end
        
        % Diffusion: Right neighbor
        if i < ni
            idx = idx + 1;
            row(idx) = k;
            col(idx) = k + 1;
            val(idx) = -epsilon / h^2;
        end
        
        % Diffusion: Bottom neighbor
        if j > 1
            idx = idx + 1;
            row(idx) = k;
            col(idx) = k - ni;
            val(idx) = -epsilon / h^2;
        end
        
        % Diffusion: Top neighbor
        if j < ni
            idx = idx + 1;
            row(idx) = k;
            col(idx) = k + ni;
            val(idx) = -epsilon / h^2;
        end
        
        % Advection part: β·∇u using upwind differencing
        % β₁ ∂u/∂x (x-direction advection)
        if beta(1) > 0
            % Upwind: use left difference
            if i > 1
                % Current point coefficient
                val(find(row(1:idx) == k & col(1:idx) == k, 1)) = ...
                    val(find(row(1:idx) == k & col(1:idx) == k, 1)) + beta(1)/h;
                % Left neighbor coefficient  
                left_idx = find(row(1:idx) == k & col(1:idx) == k-1, 1);
                if ~isempty(left_idx)
                    val(left_idx) = val(left_idx) - beta(1)/h;
                else
                    idx = idx + 1;
                    row(idx) = k;
                    col(idx) = k - 1;
                    val(idx) = -beta(1)/h;
                end
            end
        else
            % Downwind: use right difference
            if i < ni
                % Current point coefficient
                val(find(row(1:idx) == k & col(1:idx) == k, 1)) = ...
                    val(find(row(1:idx) == k & col(1:idx) == k, 1)) - beta(1)/h;
                % Right neighbor coefficient
                right_idx = find(row(1:idx) == k & col(1:idx) == k+1, 1);
                if ~isempty(right_idx)
                    val(right_idx) = val(right_idx) + beta(1)/h;
                else
                    idx = idx + 1;
                    row(idx) = k;
                    col(idx) = k + 1;
                    val(idx) = beta(1)/h;
                end
            end
        end
        
        % β₂ ∂u/∂y (y-direction advection)
        if beta(2) > 0
            % Upwind: use bottom difference
            if j > 1
                % Current point coefficient
                val(find(row(1:idx) == k & col(1:idx) == k, 1)) = ...
                    val(find(row(1:idx) == k & col(1:idx) == k, 1)) + beta(2)/h;
                % Bottom neighbor coefficient
                bottom_idx = find(row(1:idx) == k & col(1:idx) == k-ni, 1);
                if ~isempty(bottom_idx)
                    val(bottom_idx) = val(bottom_idx) - beta(2)/h;
                else
                    idx = idx + 1;
                    row(idx) = k;
                    col(idx) = k - ni;
                    val(idx) = -beta(2)/h;
                end
            end
        else
            % Downwind: use top difference
            if j < ni
                % Current point coefficient
                val(find(row(1:idx) == k & col(1:idx) == k, 1)) = ...
                    val(find(row(1:idx) == k & col(1:idx) == k, 1)) - beta(2)/h;
                % Top neighbor coefficient
                top_idx = find(row(1:idx) == k & col(1:idx) == k+ni, 1);
                if ~isempty(top_idx)
                    val(top_idx) = val(top_idx) + beta(2)/h;
                else
                    idx = idx + 1;
                    row(idx) = k;
                    col(idx) = k + ni;
                    val(idx) = beta(2)/h;
                end
            end
        end
    end
end

% Trim arrays
row = row(1:idx);
col = col(1:idx);
val = val(1:idx);

% Create sparse matrix
A = sparse(row, col, val, N, N);

% Create right-hand side: f(x,y) = 1 (constant source)
b = ones(N, 1) * h^2;

% Check non-symmetry
symmetry_measure = norm(A - A', 'fro') / norm(A, 'fro');
fprintf('  Matrix properties:\n');
fprintf('    nnz = %d\n', nnz(A));
fprintf('    Non-symmetry measure: %.2e\n', symmetry_measure);
fprintf('    %s\n', symmetry_measure < 1e-12);
fprintf('  Expected behavior:\n');
if Pe > 50
    fprintf('    AMG: Poor convergence (advection-dominated)\n');
    fprintf('    ILU: Should work well\n');
else
    fprintf('    AMG: Moderate convergence\n');
    fprintf('    ILU: Excellent performance\n');
end

end