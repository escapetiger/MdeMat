
% Sample vertices on unit square
V = [0,0; 1,0; 0,1; 1,1; 0.5,0.5];

% Periodic shifts for 4 edges (left-right, right-left, bottom-top, top-bottom)
periodicShifts = [
    1, 0;   % shift right -> left
   -1, 0;   % shift left -> right
    0, 1;   % shift top -> bottom
    0, -1;  % shift bottom -> top
];

tol = 1e-10;

vertexToVertexTable = buildVertexToVertexTable(V, periodicShifts, tol);

% Display pairs
[i,j] = find(vertexToVertexTable);
disp([i,j])


function vertexToVertexTable = buildVertexToVertexTable(V, periodicShifts, tol)
% Build a sparse binary vertex-to-vertex table indicating periodic pairs.
%
% Inputs:
%   V - nVertices x dim matrix of vertex coordinates
%   periodicShifts - mShifts x dim matrix of periodic shift vectors
%   tol - tolerance for matching vertices (e.g., 1e-10)
%
% Output:
%   vertexToVertexTable - nVertices x nVertices sparse binary matrix
%                         vertexToVertexTable(i,j) = 1 means vertex i
%                         corresponds periodically to vertex j

nVertices = size(V,1);
vertexToVertexTable = sparse(nVertices, nVertices);

for k = 1:size(periodicShifts,1)
    shift = periodicShifts(k, :);
    
    % Shift vertices by negative shift (bring vertices from periodic neighbor domain back)
    Vshifted = V - shift;
    
    % Compute pairwise distances between original and shifted vertices
    D = pdist2(V, Vshifted);
    
    % Find matching pairs (distance less than tolerance)
    [I, J] = find(D < tol);
    
    % Update sparse adjacency matrix for periodic vertex pairs
    % Set vertexToVertexTable(i,j) = 1 if vertex i matches shifted vertex j
    % Also mark symmetric to keep undirected mapping
    vertexToVertexTable = vertexToVertexTable + sparse(I, J, 1, nVertices, nVertices);
    vertexToVertexTable = vertexToVertexTable + sparse(J, I, 1, nVertices, nVertices);
end

% Ensure binary matrix (in case of multiple shifts overlap)
vertexToVertexTable = spones(vertexToVertexTable);

end
