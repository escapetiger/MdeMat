classdef GridDirichletBoundaryFluxAssembly < fem.assembly.GridAssembly
    % GRIDDIRICHLETBOUNDARYASSEMBLY Dirichlet assembly on the boundary
    % elements.

    properties (Constant)
        TYPE = 'dirichlet'
    end

    methods
        function T = assembleInnerLeftFlux(obj, d, c, mode)
            % ASSEMBLEINNERLEFTFLUX Inner contribution to a left flux on
            % the boundary along a specific dimension.
            switch mode
                case 1
                    T = obj.assembleUpwindInnerLeftFlux(d, c);
                case 2
                    T = obj.assembleAuxiliaryInnerLeftFlux(d, c);
                case 3
                    T = obj.assemblePrimalInnerLeftFlux(d, c);
            end
        end

        function T = assembleOuterLeftFlux(obj, d, c, mode)
            % ASSEMBLEOUTERLEFTFLUX Outer contribution to a left flux on the
            % boundary along a specific dimension.
            switch mode
                case 1
                    T = obj.assembleUpwindOuterLeftFlux(d, c);
                case 2
                    T = obj.assembleAuxiliaryOuterLeftFlux(d, c);
                case 3
                    T = obj.assemblePrimalOuterLeftFlux(d, c);
            end
        end

        function T = assembleInnerRightFlux(obj, d, c, mode)
            % ASSEMBLEINNERRIGHTFLUX Inner contribution to a right flux on the
            % boundary along a specific dimension.
            switch mode
                case 1
                    T = obj.assembleUpwindInnerRightFlux(d, c);
                case 2
                    T = obj.assembleAuxiliaryInnerRightFlux(d, c);
                case 3
                    T = obj.assemblePrimalInnerRightFlux(d, c);
            end
        end

        function T = assembleOuterRightFlux(obj, d, c, mode)
            % ASSEMBLEOUTERRIGHTFLUX Outer contribution to a right flux on the
            % boundary along a specific dimension.
            switch mode
                case 1
                    T = obj.assembleUpwindOuterRightFlux(d, c);
                case 2
                    T = obj.assembleAuxiliaryOuterRightFlux(d, c);
                case 3
                    T = obj.assemblePrimalOuterRightFlux(d, c);
            end
        end

        function T = assembleLeftTrace(obj, d, c, f)
            % ASSEMBLELEFTTRACE Assemble left Dirichlet boundary condition.
            %
            %   T = assembleLeftTrace(obj, d, c) assembles the Dirichlet
            %   boundary condition on the left boundary along the dimenison
            %   @a d with the coefficient @a c.
            %
            % Inputs:
            %   obj - The GridDirichletBoundaryAssembly object
            %   d - Dimension index
            %   c - Coefficient vector
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
            U = f(X);
            U = reshape(U, np, []);
            Q = P.project(-F*U);
            Q = reshape(Q, np*ne, []);
            nq = size(Q, 2);
            k = 1:numel(Q);
            T(k, 1) = repmat(reshape(bsxfun(@plus, (l(:).' - 1)*np, (1:np).'), [], 1), nq, 1);
            T(k, 2) = reshape(repmat(1:nq, ne*np, 1), [], 1);
            T(k, 3) = reshape(kron(c(l), ones(np, 1)).*Q, [], 1);
        end

        function T = assembleRightTrace(obj, d, c, f)
            % ASSEMBLERIGHTTRACE Assemble right Dirichlet condition.
            %
            %   T = assembleRightTrace(obj, d, c) assembles the Dirichlet
            %   boundary condition on the right boundary along the
            %   dimenison @a d with the coefficient @a c.
            %
            % Inputs:
            %   obj - The GridDirichletBoundaryAssembly object
            %   d - Dimension index
            %   c - Coefficient vector
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
            U = f(X);
            U = reshape(U, np, []);
            Q = P.project(-F*U);
            Q = reshape(Q, np*ne, []);
            nq = size(Q, 2);
            k = 1:numel(Q);
            T(k, 1) = repmat(reshape(bsxfun(@plus, (l(:).' - 1)*np, (1:np).'), [], 1), nq, 1);
            T(k, 2) = reshape(repmat(1:nq, ne*np, 1), [], 1);
            T(k, 3) = reshape(kron(c(l), ones(np, 1)).*Q, [], 1);
        end

        function T = assembleLeftJump(obj, d, c, f, g)
            % ASSEMBLELEFTJUMP Assemble left boundary jump.
            %
            %   T = assembleLeftJump(obj, d, c) assembles the jump between
            %   the trace and the extraction on the left boundary along the
            %   dimenison @a d with the coefficient @a c.
            %
            % Inputs:
            %   obj - The GridDirichletBoundaryAssembly object
            %   d - Dimension index
            %   c - Coefficient vector
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
            U = f(X);
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

        function T = assembleRightJump(obj, d, c, f, g)
            % ASSEMBLERIGHTJUMP Assemble right boundary jump.
            %
            %   T = assembleRightJump(obj, d, c) assembles the jump between
            %   the trace and the extraction on the right boundary along
            %   the dimenison @a d with the coefficient @a c.
            %
            % Inputs:
            %   obj - The GridDirichletBoundaryAssembly object
            %   d - Dimension index
            %   c - Coefficient vector
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
            U = f(X);
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

        function T = assembleImplicitLeftJump(obj, d, c)
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

        function T = assembleImplicitRightJump(obj, d, c)
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

    methods (Access = private)
        function T = assembleAuxiliaryInnerLeftFlux(obj, d, c)
            % ASSEMBLEINNERLEFTFLUX Assemble inner contribution to a left
            % boundary flux.
            %
            %   T = assembleInnerLeftFlux(obj, d, c) computes the inner
            %   contribution to a left flux along the specified dimension
            %   on boundaries with Dirichlet boundary conditions for
            %   auxiliary variables.
            %
            % Inputs:
            %   obj - The AuxiliaryGridAssembly object
            %   d - Dimension index
            %   c - Coefficient vector
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

            F = G.fluxData(1, 2*d).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l), v(:));
        end

        function T = assembleAuxiliaryOuterLeftFlux(obj, d, c)
            % ASSEMBLEOUTERLEFTFLUX Assemble outer contribution to a left
            % boundary flux.
            %
            %   T = assembleOuterLeftFlux(obj, d, c) computes the outer
            %   contribution to a left flux along the specified dimension
            %   on boundaries with Dirichlet boundary conditions. Combines
            %   interior and extraction contributions.
            %
            % Inputs:
            %   obj - The AuxiliaryGridAssembly object
            %   d - Dimension index
            %   c - Coefficient vector
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

        function T = assembleAuxiliaryInnerRightFlux(obj, d, c)
            % ASSEMBLEINNERRIGHTFLUX Assemble inner contribution to a right
            % boundary flux.
            %
            %   T = assembleInnerRightFlux(obj, d, c) computes the inner
            %   contribution to a right flux along the specified dimension
            %   on boundaries with Dirichlet boundary conditions for
            %   auxiliary variables.
            %
            % Inputs:
            %   obj - The AuxiliaryGridAssembly object
            %   d - Dimension index
            %   c - Coefficient vector
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

        function T = assembleAuxiliaryOuterRightFlux(obj, d, c)
            % ASSEMBLEOUTERRIGHTFLUX Assemble outer contribution to a right
            % boundary flux.
            %
            %   T = assembleOuterRightFlux(obj, d, c) computes the outer
            %   contribution to a right flux along the specified dimension
            %   on boundaries with Dirichlet boundary conditions. Combines
            %   interior and extraction contributions.
            %
            % Inputs:
            %   obj - The AuxiliaryGridAssembly object
            %   d - Dimension index
            %   c - Coefficient vector
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
        
        function T = assemblePrimalInnerLeftFlux(obj, d, c)
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

        function T = assemblePrimalOuterLeftFlux(obj, d, c)
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

        function T = assemblePrimalInnerRightFlux(obj, d, c)
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

        function T = assemblePrimalOuterRightFlux(obj, d, c)
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

        function T = assembleUpwindInnerLeftFlux(obj, d, c)
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

        function T = assembleUpwindOuterLeftFlux(obj, d, c)
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

        function T = assembleUpwindInnerRightFlux(obj, d, c)
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

        function T = assembleUpwindOuterRightFlux(obj, d, c)
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
