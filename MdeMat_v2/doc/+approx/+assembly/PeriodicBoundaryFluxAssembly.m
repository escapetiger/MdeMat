classdef PeriodicBoundaryFluxAssembly < approx.assembly.Assembly
    % PERIODICBOUNDARYFLUXASSEMBLY Assembly for periodic boundary flux
    % terms.
    %
    %   PeriodicBoundaryFluxAssembly handles the assembly of flux terms
    %   at boundaries with periodic boundary conditions in discontinuous
    %   Galerkin methods. Periodic boundary conditions enforce continuity
    %   between opposite boundaries of the computational domain, creating
    %   a seamless connection that simulates an infinite or cyclic domain.
    %
    %   The class implements flux assembly for elements at both the
    %   minimum and maximum boundaries in each coordinate direction,
    %   ensuring proper coupling between these boundaries through
    %   numerical fluxes. This approach maintains the conservative
    %   properties of DG methods while enforcing periodicity.
    %
    %   Periodic boundaries are particularly useful for simulating
    %   flows in periodic domains, wave propagation problems, and
    %   cases where artificial boundary effects need to be minimized.
    %
    % See also:
    %   approx.assembly.Assembly, 
    %   approx.assembly.DirichletBoundaryFluxAssembly

    methods
        function T = assembleInnerLeftFlux(obj, dim, coe)
            % ASSEMBLEINNERLEFTFLUX Assemble inner left flux on periodic
            % boundaries.
            %
            %   T = assembleInnerLeftFlux(obj, dim, coe) assembles the
            %   inner contribution to left-biased flux terms at periodic
            %   boundaries along the specified dimension. This method
            %   handles elements at both minimum and maximum boundaries.
            %
            % Inputs:
            %   obj - The PeriodicBoundaryFluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val]

            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == 1 | m(:, dim) == obj.space.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.space.mesh.indexer.multiToLinear(m, 0);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l), v(:));
        end

        function T = assembleInnerRightFlux(obj, dim, coe)
            % ASSEMBLEINNERRIGHTFLUX Assemble inner right flux on periodic
            % boundaries.
            %
            %   T = assembleInnerRightFlux(obj, dim, coe) assembles the
            %   inner contribution to right-biased flux terms at periodic
            %   boundaries along the specified dimension.
            %
            % Inputs:
            %   obj - The PeriodicBoundaryFluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val]

            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == 1 | m(:, dim) == obj.space.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.space.mesh.indexer.multiToLinear(m, 0);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l), v(:));
        end

        function T = assembleOuterLeftFlux(obj, dim, coe)
            % ASSEMBLEOUTERLEFTFLUX Assemble outer left flux on periodic
            % boundaries.
            %
            %   T = assembleOuterLeftFlux(obj, dim, coe) assembles the
            %   outer contribution to left-biased flux terms at periodic
            %   boundaries. This method implements the periodic coupling
            %   by connecting elements at opposite boundaries.
            %
            % Inputs:
            %   obj - The PeriodicBoundaryFluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val]

            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == 1 | m(:, dim) == obj.space.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.space.mesh.indexer.multiToLinear(m, 0);
            m(:, dim) = m(:, dim) - 1;
            l2 = obj.space.mesh.indexer.multiToLinear(m, 0);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(2, 2*dim-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l1), v(:));
        end

        function T = assembleOuterRightFlux(obj, dim, coe)
            % ASSEMBLEOUTERRIGHTFLUX Assemble outer right flux on periodic
            % boundaries.
            %
            %   T = assembleOuterRightFlux(obj, dim, coe) assembles the
            %   outer contribution to right-biased flux terms at periodic
            %   boundaries. This completes the periodic coupling by
            %   handling the forward direction connections.
            %
            % Inputs:
            %   obj - The PeriodicBoundaryFluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val]

            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == 1 | m(:, dim) == obj.space.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.space.mesh.indexer.multiToLinear(m, 0);
            m(:, dim) = m(:, dim) + 1;
            l2 = obj.space.mesh.indexer.multiToLinear(m, 0);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(2, 2*dim).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l1), v(:));
        end
    end
end