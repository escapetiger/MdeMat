classdef DirichletBoundaryFluxAssembly < fem.assembly.Assembly
    % DIRICHLETBOUNDARYASSEMBLY Dirichlet assembly on the boundary
    % elements.

    properties
        fluxMode % Flux mode: {'upwind', 'auxiliary', 'primal'}
    end

    methods
        function obj = DirichletBoundaryFluxAssembly(context, fluxMode)
            % DIRICHLETBOUNDARYFLUXASSEMBLY Constructor for
            % DirichletBoundaryFluxAssembly.
            %
            %   obj = DirichletBoundaryFluxAssembly(context) creates an
            %   assembly object for the Dirichelt boundary flux with the
            %   specified flux mode.
            %
            % Inputs:
            %   context - Context object
            %   fluxMode - Flux mode: {'upwind', 'auxiliary', 'primal'}
            %
            % Outputs:
            %   obj - Constructed DirichletBoundaryFluxAssembly object

            if nargin < 2, fluxMode = []; end

            obj@fem.assembly.Assembly(context);
            obj.setFluxMode(fluxMode);
        end

        function obj = setFluxMode(obj, fluxMode)
            % SETFLUXMODE Set the flux mode.
            %
            %   obj = setFluxMode(obj, fluxMode) set the flux mode.
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
            % ASSEMBLEINNERLEFTFLUX Inner contribution to a left flux on
            % the boundary along a specific dimension.

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
            % ASSEMBLEOUTERLEFTFLUX Outer contribution to a left flux on the
            % boundary along a specific dimension.

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
            % ASSEMBLEINNERRIGHTFLUX Inner contribution to a right flux on the
            % boundary along a specific dimension.
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
            % ASSEMBLEOUTERRIGHTFLUX Outer contribution to a right flux on
            % the boundary along a specific dimension.

            switch lower(obj.fluxMode)
                case 'upwind'
                    T = obj.assembleUpwindOuterRightFlux(dim, coe);
                case 'auxiliary'
                    T = obj.assembleAuxiliaryOuterRightFlux(dim, coe);
                case 'primal'
                    T = obj.assemblePrimalOuterRightFlux(dim, coe);
            end
        end

        function T = assembleLeftTrace(obj, dim, coe, f, varargin)
            % ASSEMBLELEFTTRACE Assemble left Dirichlet boundary condition.
            %
            %   T = assembleLeftTrace(obj, dim, coe) assembles the
            %   Dirichlet boundary condition on the left boundary along the
            %   dimenison @a dim with the coefficient @a coe.
            %
            % Inluts:
            %   obj - The DirichletBoundaryFluxAssembly object
            %   dim - Dimension index
            %   coe - Coefficient vector
            %
            % Outputs:
            %   T - Triplet matrix for boundary condition terms

            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l = obj.context.mesh.indexer.multiToLinear(m, 1);
            X = obj.context.mesh.collocate(obj.context.fe.fluxData(2*dim-1).nodes, m);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim-1).vector;
            U = f(X, varargin{:});
            U = reshape(U, size(F, 2), []);
            Q = P.project(-F*U);
            Q = reshape(Q, nl*ne, []);
            nq = size(Q, 2);
            k = 1:numel(Q);
            T(k, 1) = repmat(reshape(bsxfun(@plus, (l(:).' - 1)*nl, (1:nl).'), [], 1), nq, 1);
            T(k, 2) = reshape(repmat(1:nq, ne*nl, 1), [], 1);
            T(k, 3) = reshape(kron(coe(l), ones(nl, 1)).*Q, [], 1);
        end

        function T = assembleRightTrace(obj, dim, coe, f, varargin)
            % ASSEMBLERIGHTTRACE Assemble right Dirichlet condition.
            %
            %   T = assembleRightTrace(obj, dim, coe) assembles the
            %   Dirichlet boundary condition on the right boundary along
            %   the dimenison @a dim with the coefficient @a coe.
            %
            % Inluts:
            %   obj - The DirichletBoundaryFluxAssembly object
            %   dim - Dimension index
            %   coe - Coefficient vector
            %
            % Outputs:
            %   T - Triplet matrix for boundary condition terms

            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == obj.context.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.context.mesh.indexer.multiToLinear(m, 1);
            X = obj.context.mesh.collocate(obj.context.fe.fluxData(2*dim).nodes, m);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim).vector;
            U = f(X, varargin{:});
            U = reshape(U, size(F, 2), []);
            Q = P.project(-F*U);
            Q = reshape(Q, nl*ne, []);
            nq = size(Q, 2);
            k = 1:numel(Q);
            T(k, 1) = repmat(reshape(bsxfun(@plus, (l(:).' - 1)*nl, (1:nl).'), [], 1), nq, 1);
            T(k, 2) = reshape(repmat(1:nq, ne*nl, 1), [], 1);
            T(k, 3) = reshape(kron(coe(l), ones(nl, 1)).*Q, [], 1);
        end

        function T = assembleLeftJump(obj, dim, coe, f, g, varargin)
            % ASSEMBLELEFTJUMP Assemble left boundary jump.
            %
            %   T = assembleLeftJump(obj, dim, coe) assembles the jump
            %   between the trace and the extraction on the left boundary
            %   along the dimenison @a dim with the coefficient @a coe.
            %
            % Inluts:
            %   obj - The DirichletBoundaryFluxAssembly object
            %   dim - Dimension index
            %   coe - Coefficient vector
            %
            % Outputs:
            %   T - Triplet matrix for boundary jump terms

            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l = obj.context.mesh.indexer.multiToLinear(m, 1);
            X = obj.context.mesh.collocate(obj.context.fe.fluxData(2*dim-1).nodes, m);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim-1).vector;
            U = f(X, varargin{:});
            U = reshape(U, size(F, 2), []);
            V = g(X).';
            V = reshape(V, size(F, 2), []);
            V = V(:, l);
            Q = P.project(-F*(V - U));
            nq = numel(Q);
            k = 1:nq;
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, (1:nl).'), [], 1);
            T(k, 2) = ones(size(T(k, 1)));
            T(k, 3) = kron(coe(l), ones(nl, 1)) .* Q(:);
        end

        function T = assembleRightJump(obj, dim, coe, f, g, varargin)
            % ASSEMBLERIGHTJUMP Assemble right boundary jump.
            %
            %   T = assembleRightJump(obj, dim, coe) assembles the jump
            %   between the trace and the extraction on the right boundary
            %   along the dimenison @a dim with the coefficient @a coe.
            %
            % Inluts:
            %   obj - The DirichletBoundaryFluxAssembly object
            %   dim - Dimension index
            %   coe - Coefficient vector
            %
            % Outputs:
            %   T - Triplet matrix for boundary jump terms

            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == obj.context.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.context.mesh.indexer.multiToLinear(m, 1);
            X = obj.context.mesh.collocate(obj.context.fe.fluxData(2*dim).nodes, m);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim).vector;
            U = f(X, varargin{:});
            U = reshape(U, size(F, 2), []);
            V = g(X).';
            V = reshape(V, size(F, 2), []);
            V = V(:, l);
            Q = P.project(F*(U - V));
            nq = numel(Q);
            k = 1:nq;
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, (1:nl).'), [], 1);
            T(k, 2) = ones(size(T(k, 1)));
            T(k, 3) = kron(coe(l), ones(nl, 1)) .* Q(:);
        end

        function T = assembleImplicitLeftJump(obj, dim, coe)
            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l = obj.context.mesh.indexer.multiToLinear(m, 1);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l), v(:));
        end

        function T = assembleImplicitRightJump(obj, dim, coe)
            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == obj.context.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.context.mesh.indexer.multiToLinear(m, 1);

            T = zeros(nl^2*ne, 3);

            F = G.fluxData(1, 2*dim).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*nl, j(:)), [], 1);
            T(k, 3) = kron(coe(l), v(:));
        end
    end

    methods (Access = private)
        function T = assembleAuxiliaryInnerLeftFlux(obj, dim, coe)
            % ASSEMBLEINNERLEFTFLUX Assemble inner contribution to a left
            % boundary flux.
            %
            %   T = assembleInnerLeftFlux(obj, dim, coe) computes the inner
            %   contribution to a left flux along the specified dimension
            %   on boundaries with Dirichlet boundary conditions for
            %   auxiliary variables.
            %
            % Inluts:
            %   obj - The AuxiliaryGridAssembly object
            %   dim - Dimension index
            %   coe - Coefficient vector
            %
            % Outputs:
            %   T - Triplet matrix for boundary contribution

            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == 1 | m(:, dim) == obj.context.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.context.mesh.indexer.multiToLinear(m, 1);

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
            % ASSEMBLEOUTERLEFTFLUX Assemble outer contribution to a left
            % boundary flux.
            %
            %   T = assembleOuterLeftFlux(obj, dim, coe) computes the outer
            %   contribution to a left flux along the specified dimension
            %   on boundaries with Dirichlet boundary conditions. Combines
            %   interior and extraction contributions.
            %
            % Inluts:
            %   obj - The AuxiliaryGridAssembly object
            %   dim - Dimension index
            %   coe - Coefficient vector
            %
            % Outputs:
            %   T - Combined triplet matrix for boundary contribution

            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            %< part 1: interior
            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == obj.context.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.context.mesh.indexer.multiToLinear(m, 1);
            m(:, dim) = m(:, dim) - 1;
            l2 = obj.context.mesh.indexer.multiToLinear(m, 1);

            T1 = zeros(nl^2*ne, 3);

            F = G.fluxData(2, 2*dim-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T1(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*nl, i(:)), [], 1);
            T1(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*nl, j(:)), [], 1);
            T1(k, 3) = kron(coe(l1), v(:));

            %< part 2: extraction
            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l = obj.context.mesh.indexer.multiToLinear(m, 1);
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
            % ASSEMBLEINNERRIGHTFLUX Assemble inner contribution to a right
            % boundary flux.
            %
            %   T = assembleInnerRightFlux(obj, dim, coe) computes the inner
            %   contribution to a right flux along the specified dimension
            %   on boundaries with Dirichlet boundary conditions for
            %   auxiliary variables.
            %
            % Inluts:
            %   obj - The AuxiliaryGridAssembly object
            %   dim - Dimension index
            %   coe - Coefficient vector
            %
            % Outputs:
            %   T - Triplet matrix for boundary contribution

            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == 1 | m(:, dim) == obj.context.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.context.mesh.indexer.multiToLinear(m, 1);

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
            % ASSEMBLEOUTERRIGHTFLUX Assemble outer contribution to a right
            % boundary flux.
            %
            %   T = assembleOuterRightFlux(obj, dim, coe) computes the outer
            %   contribution to a right flux along the specified dimension
            %   on boundaries with Dirichlet boundary conditions. Combines
            %   interior and extraction contributions.
            %
            % Inluts:
            %   obj - The AuxiliaryGridAssembly object
            %   dim - Dimension index
            %   coe - Coefficient vector
            %
            % Outputs:
            %   T - Combined triplet matrix for boundary contribution

            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            %< part 1: interior
            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.context.mesh.indexer.multiToLinear(m, 1);
            m(:, dim) = m(:, dim) + 1;
            l2 = obj.context.mesh.indexer.multiToLinear(m, 1);

            T1 = zeros(nl^2*ne, 3);

            F = G.fluxData(2*dim, 2).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T1(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*nl, i(:)), [], 1);
            T1(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*nl, j(:)), [], 1);
            T1(k, 3) = kron(coe(l1), v(:));

            %< part 2: extraction
            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == obj.context.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.context.mesh.indexer.multiToLinear(m, 1);
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
            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l = obj.context.mesh.indexer.multiToLinear(m, 1);

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
            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == obj.context.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.context.mesh.indexer.multiToLinear(m, 1);
            m(:, dim) = m(:, dim) - 1;
            l2 = obj.context.mesh.indexer.multiToLinear(m, 1);

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
            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == obj.context.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.context.mesh.indexer.multiToLinear(m, 1);

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
            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.context.mesh.indexer.multiToLinear(m, 1);
            m(:, dim) = m(:, dim) + 1;
            l2 = obj.context.mesh.indexer.multiToLinear(m, 1);

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
            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == 1 | m(:, dim) == obj.context.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.context.mesh.indexer.multiToLinear(m, 1);

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
            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == obj.context.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.context.mesh.indexer.multiToLinear(m, 1);
            m(:, dim) = m(:, dim) - 1;
            l2 = obj.context.mesh.indexer.multiToLinear(m, 1);

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
            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == 1 | m(:, dim) == obj.context.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.context.mesh.indexer.multiToLinear(m, 1);

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
            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.context.mesh.indexer.multiToLinear(m, 1);
            m(:, dim) = m(:, dim) + 1;
            l2 = obj.context.mesh.indexer.multiToLinear(m, 1);

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
