classdef DirichletBoundaryFluxAssembly < approx.assembly.Assembly
    % DIRICHLETBOUNDARYFLUXASSEMBLY Assembly for Dirichlet boundary flux
    % terms.
    %
    %   DirichletBoundaryFluxAssembly handles the assembly of flux terms
    %   at boundaries with Dirichlet boundary conditions in discontinuous
    %   Galerkin methods. It supports multiple flux modes (upwind,
    %   auxiliary, primal) that provide different approaches for enforcing
    %   Dirichlet conditions and ensuring stability.
    %
    %   The class implements various flux formulations:
    %   - Upwind flux: provides stability through upwind bias
    %   - Auxiliary flux: auxiliary variable in LDG method
    %   - Primal flux: primal variable in LDG method
    %
    %   Each flux mode handles both trace assembly (for prescribed boundary
    %   values) and jump assembly (for boundary condition enforcement) with
    %   appropriate treatment of inner and outer flux contributions.
    %
    % See also:
    %   approx.assembly.Assembly,
    %   approx.assembly.PeriodicBoundaryFluxAssembly

    properties
        fluxMode % Flux mode: {'upwind', 'auxiliary', 'primal'}
    end

    methods
        function obj = DirichletBoundaryFluxAssembly(space, operator, fluxMode)
            % DIRICHLETBOUNDARYFLUXASSEMBLY Constructor for
            % DirichletBoundaryFluxAssembly.
            %
            %   obj = DirichletBoundaryFluxAssembly(space, operator)
            %   creates an assembly object for Dirichlet boundary flux with
            %   default flux mode.
            %
            %   obj = DirichletBoundaryFluxAssembly(space, operator,
            %   fluxMode) creates an assembly object with the specified
            %   flux mode.
            %
            % Inputs:
            %   space - MeshSpace object
            %   operator - ElementOperator object
            %   fluxMode - Flux mode: {'upwind', 'auxiliary', 'primal'} (optional)
            %
            % Outputs:
            %   obj - Constructed DirichletBoundaryFluxAssembly object

            if nargin < 3, fluxMode = []; end

            obj@approx.assembly.Assembly(space, operator);
            obj.setFluxMode(fluxMode);
        end

        function obj = setFluxMode(obj, fluxMode)
            % SETFLUXMODE Set the flux mode.
            %
            %   obj = setFluxMode(obj, fluxMode) sets the flux mode for
            %   Dirichlet boundary condition enforcement. Different modes
            %   provide different approaches for stability and accuracy.
            %
            % Inputs:
            %   obj - The DirichletBoundaryFluxAssembly object
            %   fluxMode - Flux mode: {'upwind', 'auxiliary', 'primal'}
            %
            % Outputs:
            %   obj - The DirichletBoundaryFluxAssembly object

            core.except.assert(isempty(fluxMode) || ...
                ismember(fluxMode, {'upwind', 'auxiliary', 'primal'}), ...
                'InvalidInput', 'Flux mode is not supported.');

            obj.fluxMode = fluxMode;
        end

        function T = assembleInnerLeftFlux(obj, dim, coe)
            % ASSEMBLEINNERLEFTFLUX Assemble inner left flux on boundaries.
            %
            %   T = assembleInnerLeftFlux(obj, dim, coe) assembles the
            %   inner contribution to left-biased flux terms at Dirichlet
            %   boundaries along the specified dimension. The specific
            %   implementation depends on the flux mode.
            %
            % Inputs:
            %   obj - The DirichletBoundaryFluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val] for inner left flux terms

            switch lower(obj.fluxMode)
                case 'upwind'
                    T = obj.assembleUpwindInnerLeftFlux(dim, coe);
                case 'auxiliary'
                    T = obj.assembleAuxiliaryInnerLeftFlux(dim, coe);
                case 'primal'
                    T = obj.assemblePrimalInnerLeftFlux(dim, coe);
            end
        end

        function T = assembleOuterLeftFlux(obj, dim, coe)
            % ASSEMBLEOUTERLEFTFLUX Assemble outer left flux on boundaries.
            %
            %   T = assembleOuterLeftFlux(obj, dim, coe) assembles the
            %   outer contribution to left-biased flux terms at Dirichlet
            %   boundaries along the specified dimension.
            %
            % Inputs:
            %   obj - The DirichletBoundaryFluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val] for outer left flux terms

            switch lower(obj.fluxMode)
                case 'upwind'
                    T = obj.assembleUpwindOuterLeftFlux(dim, coe);
                case 'auxiliary'
                    T = obj.assembleAuxiliaryOuterLeftFlux(dim, coe);
                case 'primal'
                    T = obj.assemblePrimalOuterLeftFlux(dim, coe);
            end
        end

        function T = assembleInnerRightFlux(obj, dim, coe)
            % ASSEMBLEINNERRIGHTFLUX Assemble inner right flux on
            % boundaries.
            %
            %   T = assembleInnerRightFlux(obj, dim, coe) assembles the
            %   inner contribution to right-biased flux terms at Dirichlet
            %   boundaries along the specified dimension.
            %
            % Inputs:
            %   obj - The DirichletBoundaryFluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val]

            switch lower(obj.fluxMode)
                case 'upwind'
                    T = obj.assembleUpwindInnerRightFlux(dim, coe);
                case 'auxiliary'
                    T = obj.assembleAuxiliaryInnerRightFlux(dim, coe);
                case 'primal'
                    T = obj.assemblePrimalInnerRightFlux(dim, coe);
            end
        end

        function T = assembleOuterRightFlux(obj, dim, coe)
            % ASSEMBLEOUTERRIGHTFLUX Assemble outer right flux on
            % boundaries.
            %
            %   T = assembleOuterRightFlux(obj, dim, coe) assembles the
            %   outer contribution to right-biased flux terms at Dirichlet
            %   boundaries along the specified dimension.
            %
            % Inputs:
            %   obj - The DirichletBoundaryFluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val]

            switch lower(obj.fluxMode)
                case 'upwind'
                    T = obj.assembleUpwindOuterRightFlux(dim, coe);
                case 'auxiliary'
                    T = obj.assembleAuxiliaryOuterRightFlux(dim, coe);
                case 'primal'
                    T = obj.assemblePrimalOuterRightFlux(dim, coe);
            end
        end

        function T = assembleTrace(obj, i, coe, f, varargin)
            % ASSEMBLETRACE Assemble Dirichlet boundary condition.
            %
            %   T = assembleTrace(obj, i, coe, f) assembles the Dirichlet
            %   boundary condition along the dimension @a dim with the
            %   coefficient @a coe and boundary function @a f. This creates
            %   the trace terms that enforce prescribed boundary values.
            %
            % Inputs:
            %   obj - The DirichletBoundaryFluxAssembly object
            %   i - Flux index
            %   coe - Coefficient vector (nElements x 1)
            %   f - Boundary function handle
            %   varargin - Additional arguments for boundary function
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val]

            dim = ceil(i / 2);
            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            if mod(i, 2) == 1
                e = m(:, dim) == 1;
            else
                e = m(:, dim) == obj.space.mesh.resolution(dim);
            end
            ne = numel(find(e));
            m = m(e, :);
            l = obj.space.mesh.indexer.multiToLinear(m, 1);
            X = obj.space.mesh.collocate(obj.space.element.fluxData(i).nodes, m);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, i).vector;
            U = f(i, X, varargin{:});
            U = reshape(U, size(F, 2), []);
            Q = P.project(F*U);
            Q = reshape(Q, nl*ne, []);
            nq = size(Q, 2);
            k = 1:numel(Q);
            T(k, 1) = repmat(reshape(bsxfun(@plus, (l(:).' - 1)*nl, (1:nl).'), [], 1), nq, 1);
            T(k, 2) = reshape(repmat(1:nq, ne*nl, 1), [], 1);
            T(k, 3) = reshape(kron(coe(l), ones(nl, 1)).*Q, [], 1);
        end

        function T = assembleJump(obj, i, coe, f, g, varargin)
            % ASSEMBLEJUMP Assemble boundary jump terms.
            %
            %   T = assembleJump(obj, i, coe, f, g) assembles the
            %   jump between the trace and the extraction on the
            %   boundary along the dimension dim.
            %
            % Inputs:
            %   obj - The DirichletBoundaryFluxAssembly object
            %   i - Flux index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %   f - Boundary function handle
            %   g - Interior extraction function handle
            %   varargin - Additional arguments for functions
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val]

            dim = ceil(i / 2);
            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            if mod(i, 2) == 1
                e = m(:, dim) == 1;
                a = 1;
            else
                e = m(:, dim) == obj.space.mesh.resolution(dim);
                a = -1;
            end
            ne = numel(find(e));
            m = m(e, :);
            l = obj.space.mesh.indexer.multiToLinear(m, 1);
            xRef = obj.space.element.fluxData(i).nodes;
            X = obj.space.mesh.collocate(xRef, m);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, i).vector;
            U = f(i, X, varargin{:});
            U = reshape(U, size(F, 2), []);
            V = g(xRef);
            V = reshape(V, size(F, 2), obj.space.nTotalElements, []);
            V = V(:, l, :);
            V = reshape(V, size(F, 2), []);
            Q = P.project(a*F*(V - U));
            Q = reshape(Q, nl * ne, []);
            nq = size(Q, 2);
            k = 1:numel(Q);
            T(k, 1) = repmat(reshape(bsxfun(@plus, (l(:).' - 1)*nl, (1:nl).'), [], 1), nq, 1);
            T(k, 2) = reshape(repmat(1:nq, ne*nl, 1), [], 1);
            T(k, 3) = reshape(kron(coe(l), ones(nl, 1)).*Q, [], 1);
        end

        function T = assembleImplicitLeftJump(obj, dim, coe)
            % ASSEMBLEIMPLICITLEFTJUMP Assemble implicit left boundary
            % jump.
            %
            %   T = assembleImplicitLeftJump(obj, dim, coe) assembles
            %   implicit jump terms for left boundaries in matrix form.
            %
            % Inputs:
            %   obj - The DirichletBoundaryFluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val]

            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l = obj.space.mesh.indexer.multiToLinear(m, 1);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim-1).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l), v(:));
        end

        function T = assembleImplicitRightJump(obj, dim, coe)
            % ASSEMBLEIMPLICITRIGHTJUMP Assemble implicit right boundary
            % jump.
            %
            %   T = assembleImplicitRightJump(obj, dim, coe) assembles
            %   implicit jump terms for right boundaries in matrix form.
            %
            % Inputs:
            %   obj - The DirichletBoundaryFluxAssembly object
            %   dim - Dimension index (positive integer)
            %   coe - Coefficient vector (nElements x 1)
            %
            % Outputs:
            %   T - Triplet matrix [row, col, val]

            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == obj.space.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.space.mesh.indexer.multiToLinear(m, 1);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l), v(:));
        end
    end

    methods (Access = private)
        function T = assembleAuxiliaryInnerLeftFlux(obj, dim, coe)
            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == 1 | m(:, dim) == obj.space.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.space.mesh.indexer.multiToLinear(m, 1);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l), v(:));
        end

        function T = assembleAuxiliaryOuterLeftFlux(obj, dim, coe)
            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            %< part 1: interior
            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == obj.space.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.space.mesh.indexer.multiToLinear(m, 1);
            m(:, dim) = m(:, dim) - 1;
            l2 = obj.space.mesh.indexer.multiToLinear(m, 1);

            T1 = zeros(nl^2*ne, 3);

            F = G.fluxData(2, 2*dim-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T1(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*nl, i(:)), [], 1);
            T1(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*nl, j(:)), [], 1);
            T1(k, 3) = kron(coe(l1), v(:));

            %< part 2: extraction
            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l = obj.space.mesh.indexer.multiToLinear(m, 1);
            T2 = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T2(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, i(:)), [], 1);
            T2(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, j(:)), [], 1);
            T2(k, 3) = kron(coe(l), v(:));

            T = [T1; T2];
        end

        function T = assembleAuxiliaryInnerRightFlux(obj, dim, coe)
            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == 1 | m(:, dim) == obj.space.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.space.mesh.indexer.multiToLinear(m, 1);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l), v(:));
        end

        function T = assembleAuxiliaryOuterRightFlux(obj, dim, coe)
            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            %< part 1: interior
            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.space.mesh.indexer.multiToLinear(m, 1);
            m(:, dim) = m(:, dim) + 1;
            l2 = obj.space.mesh.indexer.multiToLinear(m, 1);

            T1 = zeros(nl^2*ne, 3);

            F = G.fluxData(2, 2*dim).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T1(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*nl, i(:)), [], 1);
            T1(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*nl, j(:)), [], 1);
            T1(k, 3) = kron(coe(l1), v(:));

            %< part 2: extraction
            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == obj.space.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.space.mesh.indexer.multiToLinear(m, 1);
            T2 = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T2(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, i(:)), [], 1);
            T2(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, j(:)), [], 1);
            T2(k, 3) = kron(coe(l), v(:));

            T = [T1; T2];
        end

        function T = assemblePrimalInnerLeftFlux(obj, dim, coe)
            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l = obj.space.mesh.indexer.multiToLinear(m, 1);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l), v(:));
        end

        function T = assemblePrimalOuterLeftFlux(obj, dim, coe)
            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == obj.space.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.space.mesh.indexer.multiToLinear(m, 1);
            m(:, dim) = m(:, dim) - 1;
            l2 = obj.space.mesh.indexer.multiToLinear(m, 1);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(2, 2*dim-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l1), v(:));
        end

        function T = assemblePrimalInnerRightFlux(obj, dim, coe)
            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == obj.space.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.space.mesh.indexer.multiToLinear(m, 1);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l), v(:));
        end

        function T = assemblePrimalOuterRightFlux(obj, dim, coe)
            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.space.mesh.indexer.multiToLinear(m, 1);
            m(:, dim) = m(:, dim) + 1;
            l2 = obj.space.mesh.indexer.multiToLinear(m, 1);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(2, 2*dim).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l1), v(:));
        end

        function T = assembleUpwindInnerLeftFlux(obj, dim, coe)
            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == 1 | m(:, dim) == obj.space.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.space.mesh.indexer.multiToLinear(m, 1);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l), v(:));
        end

        function T = assembleUpwindOuterLeftFlux(obj, dim, coe)
            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == obj.space.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.space.mesh.indexer.multiToLinear(m, 1);
            m(:, dim) = m(:, dim) - 1;
            l2 = obj.space.mesh.indexer.multiToLinear(m, 1);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(2, 2*dim-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l1), v(:));
        end

        function T = assembleUpwindInnerRightFlux(obj, dim, coe)
            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == 1 | m(:, dim) == obj.space.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.space.mesh.indexer.multiToLinear(m, 1);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l), v(:));
        end

        function T = assembleUpwindOuterRightFlux(obj, dim, coe)
            nl = obj.space.nLocalDofs;
            P = obj.space.element.projector;
            G = obj.operator.gradient;

            m = obj.space.mesh.allElementMultiIndices;
            e = m(:, dim) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.space.mesh.indexer.multiToLinear(m, 1);
            m(:, dim) = m(:, dim) + 1;
            l2 = obj.space.mesh.indexer.multiToLinear(m, 1);

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