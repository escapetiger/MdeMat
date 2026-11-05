classdef AuxiliaryGridAssembly < fem.assembly.MixedGridAssembly
    % AUXILIARYGRIDASSEMBLY Grid-based assembly for auxiliary variables.
    %
    %   AuxiliaryGridAssembly implements Dirichlet boundary treatment for
    %   auxiliary variables in mixed finite element systems. This class
    %   handles coupling terms between primal and auxiliary variables with
    %   both interior and extraction contributions.
    %
    %   The assembly specializes in handling jump terms and boundary
    %   conditions specific to auxiliary variable formulations in mixed
    %   methods and discontinuous Galerkin schemes.
    %
    % Examples:
    %   % Create auxiliary assembly with left flux and Dirichlet BC
    %   assembly = AuxiliaryGridAssembly(fe, mesh, op, 1, 1);
    %   
    %   % Assemble outer jump operator
    %   jumpMatrix = assembly.outerJump(1, [1.0, 0.5]);
    %
    %   % Use in divergence computation
    %   divMatrix = assembly.divergence([1.0, 0.5]);
    %
    % See also:
    %   fem.assembly.MixedGridAssembly, fem.assembly.PrimalGridAssembly
    
    methods
        function A = outerJump(obj, d, c)
            % OUTERJUMP Assemble outer jump operator for auxiliary
            % variables.
            %
            %   A = outerJump(obj, d, c) creates the outer jump operator
            %   matrix for auxiliary variables along the specified
            %   dimension. This handles jump discontinuities at element
            %   boundaries for mixed formulations.
            %
            % Inputs:
            %   obj - The AuxiliaryGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector or scalar
            %
            % Outputs:
            %   A - Assembled outer jump sparse matrix

            ne = obj.nTotalElements;

            %< ------------------------------------------------------------
            %< COEFFICIENTS
            %< ------------------------------------------------------------
            if isa(obj.mesh, 'approx.mesh.UniformGrid')
                h = obj.mesh.spacings{d};
                c = ones(1, ne) .* c ./ h;
            end

            if isa(obj.mesh, 'approx.mesh.NonuniformGrid')
                h = obj.mesh.spacings;
                [h{:}] = ndgrid(h{:});
                c = c ./ h{d}(:);
            end

            c = c(:);

            %< ------------------------------------------------------------
            %< FLUX CONTRIBUTION
            %< ------------------------------------------------------------
            if isa(obj.fe, 'fem.element.L2FiniteElement')
                switch obj.fluxType
                    case 1
                        T = obj.assembleOuterLeftFluxDirichletBoundaryJumpTerms(d, c);
                    case 2
                        T = obj.assembleOuterRightFluxDirichletBoundaryJumpTerms(d, c);
                    case 3
                        T1 = obj.assembleOuterLeftFluxDirichletBoundaryJumpTerms(d, c/2);
                        T2 = obj.assembleOuterRightFluxDirichletBoundaryJumpTerms(d, c/2);
                        T = [T1; T2];
                end
            end

            %< ------------------------------------------------------------
            %< TRIPLETS TO SPARSE MATRIX
            %< ------------------------------------------------------------
            n = obj.nGlobalDofs;
            A = obj.sparseFromTriplets(T, n, n);
        end
    end

    methods (Access = protected)
        function T = assembleInnerLeftFluxDirichletBoundaryTerms(obj, d, c)
            % ASSEMBLEINNERLEFTFLUXDIRICHLETBOUNDARYTERMS Assemble inner
            % left flux Dirichlet boundary terms.
            %
            %   T = assembleInnerLeftFluxDirichletBoundaryTerms(obj, d, c)
            %   computes the inner contribution to a left flux along the
            %   specified dimension on boundaries with Dirichlet boundary
            %   conditions for auxiliary variables.
            %
            % Inputs:
            %   obj - The AuxiliaryGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix [row, col, value] for boundary contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == 1 | m(:, d) == obj.mesh.resolution(d);
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
            %   conditions. Combines interior and extraction contributions.
            %
            % Inputs:
            %   obj - The AuxiliaryGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Combined triplet matrix for boundary contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            %< part 1: interior
            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.mesh.indexer.multiToLinear(m, 1);
            m(:, d) = m(:, d) - 1;
            l2 = obj.mesh.indexer.multiToLinear(m, 1);

            T1 = zeros(np^2*ne, 3);
            
            F = G.fluxData(2, 2*d-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T1(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*np, i(:)), [], 1);
            T1(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*np, j(:)), [], 1);
            T1(k, 3) = kron(c(l1), v(:));

            %< part 2: extraction
            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 1);
            T2 = zeros(np^2*ne, 3);
            
            F = G.fluxData(1, 2*d-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T2(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, i(:)), [], 1);
            T2(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*np, j(:)), [], 1);
            T2(k, 3) = kron(c(l), v(:));
            
            T = [T1; T2];
        end
    
        function T = assembleInnerRightFluxDirichletBoundaryTerms(obj, d, c)
            % ASSEMBLEINNERRIGHTFLUXDIRICHLETBOUNDARYTERMS Assemble inner
            % right flux Dirichlet boundary terms.
            %
            %   T = assembleInnerRightFluxDirichletBoundaryTerms(obj, d, c)
            %   computes the inner contribution to a right flux along the
            %   specified dimension on boundaries with Dirichlet boundary
            %   conditions for auxiliary variables.
            %
            % Inputs:
            %   obj - The AuxiliaryGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for boundary contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == 1 | m(:, d) == obj.mesh.resolution(d);
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
            %   conditions. Combines interior and extraction contributions.
            %
            % Inputs:
            %   obj - The AuxiliaryGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Combined triplet matrix for boundary contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            %< part 1: interior
            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.mesh.indexer.multiToLinear(m, 1);
            m(:, d) = m(:, d) + 1;
            l2 = obj.mesh.indexer.multiToLinear(m, 1);

            T1 = zeros(np^2*ne, 3);
            
            F = G.fluxData(2*d, 2).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T1(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*np, i(:)), [], 1);
            T1(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*np, j(:)), [], 1);
            T1(k, 3) = kron(c(l1), v(:));

            %< part 2: extraction
            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 1);
            T2 = zeros(np^2*ne, 3);
            
            F = G.fluxData(1, 2*d).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T2(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, i(:)), [], 1);
            T2(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*np, j(:)), [], 1);
            T2(k, 3) = kron(c(l), v(:));
            
            T = [T1; T2];
        end

        function T = assembleOuterLeftFluxDirichletBoundaryJumpTerms(obj, d, c)
            % ASSEMBLEOUTERLEFTFLUXDIRICHLETBOUNDARYJUMPTERMS Assemble
            % outer left flux Dirichlet boundary jump terms.
            %
            %   T = assembleOuterLeftFluxDirichletBoundaryJumpTerms(obj, d,
            %   c) computes the outer jump contribution to a left flux
            %   along the specified dimension on boundaries with Dirichlet
            %   boundary conditions for auxiliary variables.
            %
            % Inputs:
            %   obj - The AuxiliaryGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for boundary jump contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == 1;
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

        function T = assembleOuterRightFluxDirichletBoundaryJumpTerms(obj, d, c)
            % ASSEMBLEOUTERRIGHTFLUXDIRICHLETBOUNDARYJUMPTERMS Assemble
            % outer right flux Dirichlet boundary jump terms.
            %
            %   T = assembleOuterRightFluxDirichletBoundaryJumpTerms(obj,
            %   d, c) computes the outer jump contribution to a right flux
            %   along the specified dimension on boundaries with Dirichlet
            %   boundary conditions for auxiliary variables.
            %
            % Inputs:
            %   obj - The AuxiliaryGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for boundary jump contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 1);

            T = zeros(np^2*ne, 3);
            
            F = G.fluxData(1, 2*d).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l), v(:));
        end
    end
end