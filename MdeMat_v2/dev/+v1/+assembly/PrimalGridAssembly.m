classdef PrimalGridAssembly < fem.assembly.MixedGridAssembly
    % PRIMALGRIDASSEMBLY Grid-based assembly for primal variables.
    %
    %   PrimalGridAssembly implements Dirichlet boundary treatment for
    %   primal variables in mixed finite element systems. This class
    %   provides specific boundary flux implementations optimized for
    %   the primal variable formulation in mixed methods.
    %
    %   The assembly handles the primary variable in mixed formulations,
    %   typically representing the main unknown quantity such as pressure
    %   in flow problems or electric potential in electromagnetic problems.
    %
    % Examples:
    %   % Create primal assembly with right flux and Dirichlet BC
    %   assembly = PrimalGridAssembly(fe, mesh, op, 1, 2);
    %   
    %   % Assemble divergence operator
    %   divMatrix = assembly.divergence([1.0, 0.5]);
    %
    %   % Assemble partial derivative in specific dimension
    %   partialMatrix = assembly.partial(1, 1.5);
    %
    % See also:
    %   fem.assembly.MixedGridAssembly, fem.assembly.AuxiliaryGridAssembly
    
    methods (Access = protected)
        function T = assembleInnerLeftFluxDirichletBoundaryTerms(obj, d, c)
            % ASSEMBLEINNERLEFTFLUXDIRICHLETBOUNDARYTERMS Assemble inner
            % left flux Dirichlet boundary terms.
            %
            %   T = assembleInnerLeftFluxDirichletBoundaryTerms(obj, d, c)
            %   computes the inner contribution to a left flux along the
            %   specified dimension on boundaries with Dirichlet boundary
            %   conditions for primal variables.
            %
            % Inputs:
            %   obj - The PrimalGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix [row, col, value] for boundary contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 1);

            T = zeros(np^2*ne, 3);
            
            F = G.fluxData(1, 2*d).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l), v(:));
        end

        function T = assembleOuterLeftFluxDirichletBoundaryTerms(obj, d, c)
            % ASSEMBLEOUTERLEFTFLUXDIRICHLETBOUNDARYTERMS Assemble outer
            % left flux Dirichlet boundary terms.
            %
            %   T = assembleOuterLeftFluxDirichletBoundaryTerms(obj, d, c)
            %   computes the outer contribution to a left flux along the
            %   specified dimension on boundaries with Dirichlet boundary
            %   conditions for primal variables.
            %
            % Inputs:
            %   obj - The PrimalGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for boundary contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.mesh.indexer.multiToLinear(m, 1);
            m(:, d) = m(:, d) - 1;
            l2 = obj.mesh.indexer.multiToLinear(m, 1);

            T = zeros(np^2*ne, 3);
            
            F = G.fluxData(2, 2*d-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l1), v(:));
        end
    
        function T = assembleInnerRightFluxDirichletBoundaryTerms(obj, d, c)
            % ASSEMBLEINNERRIGHTFLUXDIRICHLETBOUNDARYTERMS Assemble inner
            % right flux Dirichlet boundary terms.
            %
            %   T = assembleInnerRightFluxDirichletBoundaryTerms(obj, d, c)
            %   computes the inner contribution to a right flux along the
            %   specified dimension on boundaries with Dirichlet boundary
            %   conditions for primal variables.
            %
            % Inputs:
            %   obj - The PrimalGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for boundary contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 1);

            T = zeros(np^2*ne, 3);
            
            F = G.fluxData(1, 2*d-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l), v(:));
        end

        function T = assembleOuterRightFluxDirichletBoundaryTerms(obj, d, c)
            % ASSEMBLEOUTERRIGHTFLUXDIRICHLETBOUNDARYTERMS Assemble outer
            % right flux Dirichlet boundary terms.
            %
            %   T = assembleOuterRightFluxDirichletBoundaryTerms(obj, d, c)
            %   computes the outer contribution to a right flux along the
            %   specified dimension on boundaries with Dirichlet boundary
            %   conditions for primal variables.
            %
            % Inputs:
            %   obj - The PrimalGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for boundary contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.mesh.indexer.multiToLinear(m, 1);
            m(:, d) = m(:, d) + 1;
            l2 = obj.mesh.indexer.multiToLinear(m, 1);

            T = zeros(np^2*ne, 3);
            
            F = G.fluxData(2, 2*d).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l1), v(:));
        end
    end
end