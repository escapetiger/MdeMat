classdef InteriorFluxAssembly < approx.assembly.Assembly
    % INTERIORFLUXASSEMBLY Assembly for interior flux terms in DG methods.
    %
    %   InteriorFluxAssembly handles the assembly of flux terms between
    %   neighboring elements in the interior of the computational domain.
    %   These flux terms are essential for providing inter-element coupling
    %   in discontinuous Galerkin (DG) methods, where the solution is
    %   allowed to be discontinuous across element boundaries.
    %
    %   The class implements both inner and outer flux contributions:
    %   - Inner fluxes: contributions from the current element's boundary
    %   - Outer fluxes: contributions from neighboring elements' boundaries
    %   
    %   Each flux type (left-biased, right-biased) is supported, enabling
    %   the construction of upwind-like or downwind-like numerical schemes.
    %   The flux terms are assembled using triplet format for efficient
    %   sparse matrix construction.
    %
    % See also:
    %   approx.assembly.Assembly, approx.assembly.FluxAssembly,
    %   approx.element.L2ElementOperator

    methods
        function T = assembleInnerLeftFlux(obj, dim, coe)
            % ASSEMBLEINNNERLEFTFLUX Assemble inner left-biased flux terms.
            %
            %   T = assembleInnerLeftFlux(obj, dim, coe) assembles the
            %   inner flux contributions for left-biased numerical fluxes
            %   along the specified dimension. Inner fluxes represent the
            %   contribution from the current element's boundary to the
            %   flux integral.
            %
            % Inputs:
            %   obj - The InteriorFluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val]

            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) > 1 & m(:, dim) < obj.space.mesh.resolution(dim);
            m = m(e, :);
            l = obj.space.mesh.indexer.multiToLinear(m, 0);
            ne = numel(l);

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
            % ASSEMBLEINNERRIGHTFLUX Assemble inner right-biased flux terms.
            %
            %   T = assembleInnerRightFlux(obj, dim, coe) assembles the
            %   inner flux contributions for right-biased numerical fluxes
            %   along the specified dimension. The negative sign in the
            %   flux matrix extraction accounts for the directional
            %   convention in right-biased schemes.
            %
            % Inputs:
            %   obj - The InteriorFluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val]

            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) > 1 & m(:, dim) < obj.space.mesh.resolution(dim);
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
            % ASSEMBLEOUTERLEFTFLUX Assemble outer left-biased flux terms.
            %
            %   T = assembleOuterLeftFlux(obj, dim, coe) assembles the
            %   outer flux contributions for left-biased numerical fluxes
            %   along the specified dimension. Outer fluxes represent the
            %   contribution from neighboring elements' boundaries to the
            %   flux integral, providing inter-element coupling.
            %
            % Inputs:
            %   obj - The InteriorFluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val]

            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) > 1 & m(:, dim) < obj.space.mesh.resolution(dim);
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
            % ASSEMBLEOUTERRIGHTFLUX Assemble outer right-biased flux terms.
            %
            %   T = assembleOuterRightFlux(obj, dim, coe) assembles the
            %   outer flux contributions for right-biased numerical fluxes
            %   along the specified dimension. The method handles the
            %   connectivity between current elements and their neighbors
            %   in the positive coordinate direction.
            %
            % Inputs:
            %   obj - The InteriorFluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val]

            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) > 1 & m(:, dim) < obj.space.mesh.resolution(dim);
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