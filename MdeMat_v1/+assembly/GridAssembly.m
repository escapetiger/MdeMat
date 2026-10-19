classdef GridAssembly < fem.assembly.Assembly
    % GRIDASSEMBLY Grid-based assembly for structured meshes.
    %
    %   GridAssembly provides assembly methods for finite element
    %   operators on structured grid meshes. This class extends the base
    %   Assembly to handle grid-specific operations and boundary condition
    %   treatments.
    %
    %   The assembly supports different boundary condition types and
    %   provides specialized methods for handling Dirichlet boundary
    %   conditions on structured grids.
    %
    % Examples:
    %   % Create grid assembly with periodic boundaries
    %   assembly = GridAssembly(fe, gridMesh, op, 0);
    %   
    %   % Create grid assembly with Dirichlet boundaries
    %   assembly = GridAssembly(fe, gridMesh, op, 1);
    %
    % See also:
    %   fem.assembly.Assembly, fem.assembly.EulerianGridAssembly

    properties
        bcType % Boundary condition type (0=periodic, 1=Dirichlet)
    end

    methods
        function obj = GridAssembly(fe, mesh, op, bcType)
            % GRIDASSEMBLY Constructor for GridAssembly.
            %
            %   obj = GridAssembly(fe, mesh, op, bcType) creates an
            %   assembly for structured grid meshes with specified boundary
            %   condition treatment.
            %
            % Inputs:
            %   fe - Finite element object
            %   mesh - Grid mesh object (must be approx.mesh.Grid)
            %   op - Finite element operator object
            %   bcType - Boundary condition type (0=periodic, 1=Dirichlet)
            %
            % Outputs:
            %   obj - Constructed GridAssembly object

            core.except.assert(isa(mesh, 'approx.mesh.Grid'), ...
                'InvalidInput', ...
                'Finite element space must be based on grid.')

            obj@fem.assembly.Assembly(fe, mesh, op);
            obj.bcType = bcType;
        end
    end
    
    methods (Access = protected)
        function T = assembleLeftDirichletBoundaryConditionTerms(obj, d, c, f, varargin)
            % ASSEMBLELEFTDIRICHLETBOUNDARYCONDITIONTERMS Assemble left
            % boundary Dirichlet condition terms.
            %
            %   T = assembleLeftDirichletBoundaryConditionTerms(obj, d, c,
            %   f) assembles the contribution from Dirichlet boundary
            %   conditions on the left boundary of the specified dimension.
            %
            % Inputs:
            %   obj - The GridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector
            %   f - Boundary condition function handle
            %   varargin - Additional arguments for boundary function
            %
            % Outputs:
            %   T - Triplet matrix for boundary condition terms

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 1);
            X = obj.mesh.collocate(obj.fe.fluxData(2*d-1).nodes, m);

            T = zeros(np^2*ne, 3);

            F = G.fluxData(1, 2*d-1).vector;
            U = f(X, varargin{:});
            U = reshape(U, size(F, 2), []);
            Q = P.project(-F*U);
            nq = numel(Q);
            k = 1:nq;
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, (1:np).'), [], 1);
            T(k, 2) = ones(size(T(k, 1)));
            T(k, 3) = kron(c(l), ones(np, 1)) .* Q(:);
        end

        function T = assembleRightDirichletBoundaryConditionTerms(obj, d, c, f, varargin)
            % ASSEMBLYIGHTDIRICHLETBOUNDARYCONDITIONTERMS Assemble right
            % boundary Dirichlet condition terms.
            %
            %   T = assembleRightDirichletBoundaryConditionTerms(obj, d, c,
            %   f) assembles the contribution from Dirichlet boundary
            %   conditions on the right boundary of the specified
            %   dimension.
            %
            % Inputs:
            %   obj - The GridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector
            %   f - Boundary condition function handle
            %   varargin - Additional arguments for boundary function
            %
            % Outputs:
            %   T - Triplet matrix for boundary condition terms

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 1);
            X = obj.mesh.collocate(obj.fe.fluxData(2*d).nodes, m);

            T = zeros(np^2*ne, 3);

            F = G.fluxData(1, 2*d).vector;
            U = f(X, varargin{:});
            U = reshape(U, size(F, 2), []);
            Q = P.project(F*U);
            nq = numel(Q);
            k = 1:nq;
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, (1:np).'), [], 1);
            T(k, 2) = ones(size(T(k, 1)));
            T(k, 3) = kron(c(l), ones(np, 1)) .* Q(:);
        end

        function T = assembleLeftDirichletBoundaryJumpTerms(obj, d, c, f, g, varargin)
            % ASSEMBLELEFTDIRICHLETBOUNDARYJUMPTERMS Assemble left boundary
            % jump terms for Dirichlet conditions.
            %
            %   T = assembleLeftDirichletBoundaryJumpTerms(obj, d, c, f, g)
            %   assembles the contribution from jump terms in Dirichlet 
            %   boundary conditions on the left boundary.
            %
            % Inputs:
            %   obj - The GridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector
            %   f - Boundary condition function handle
            %   g - Jump function handle
            %   varargin - Additional arguments
            %
            % Outputs:
            %   T - Triplet matrix for boundary jump terms

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 1);
            X = obj.mesh.collocate(obj.fe.fluxData(2*d-1).nodes, m);

            T = zeros(np^2*ne, 3);

            F = G.fluxData(1, 2*d-1).vector;
            U = f(X, varargin{:});
            U = reshape(U, size(F, 2), []);
            V = g(X).';
            V = reshape(V, size(F, 2), []);
            V = V(:, l);
            Q = P.project(-F*(V - U));
            nq = numel(Q);
            k = 1:nq;
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, (1:np).'), [], 1);
            T(k, 2) = ones(size(T(k, 1)));
            T(k, 3) = kron(c(l), ones(np, 1)) .* Q(:);
        end

        function T = assembleRightDirichletBoundaryJumpTerms(obj, d, c, f, g, varargin)
            % ASSEMBLYIGHTDIRICHLETBOUNDARYJUMPTERMS Assemble right
            % boundary jump terms for Dirichlet conditions.
            %
            %   T = assembleRightDirichletBoundaryJumpTerms(obj, d, c, f,
            %   g) assembles the contribution from jump terms in Dirichlet
            %   boundary conditions on the right boundary.
            %
            % Inputs:
            %   obj - The GridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector
            %   f - Boundary condition function handle
            %   g - Jump function handle
            %   varargin - Additional arguments
            %
            % Outputs:
            %   T - Triplet matrix for boundary jump terms

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 1);
            X = obj.mesh.collocate(obj.fe.fluxData(2*d).nodes, m);

            T = zeros(np^2*ne, 3);

            F = G.fluxData(1, 2*d).vector;
            U = f(X, varargin{:});
            U = reshape(U, size(F, 2), []);
            V = g(X).';
            V = reshape(V, size(F, 2), []);
            V = V(:, l);
            Q = P.project(F*(U - V));
            nq = numel(Q);
            k = 1:nq;
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, (1:np).'), [], 1);
            T(k, 2) = ones(size(T(k, 1)));
            T(k, 3) = kron(c(l), ones(np, 1)) .* Q(:);
        end
    end
end