classdef PeriodicBoundaryFluxAssembly < fem.assembly.Assembly   
    methods
        function T = assembleInnerLeftFlux(obj, dim, coe)
            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == 1 | m(:, dim) == obj.context.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.context.mesh.indexer.multiToLinear(m, 0);

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
            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == 1 | m(:, dim) == obj.context.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.context.mesh.indexer.multiToLinear(m, 0);

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
            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == 1 | m(:, dim) == obj.context.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.context.mesh.indexer.multiToLinear(m, 0);
            m(:, dim) = m(:, dim) - 1;
            l2 = obj.context.mesh.indexer.multiToLinear(m, 0);

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
            nl = obj.context.nLocalDofs;
            P = obj.context.fe.projector;
            G = obj.context.feOp.gradient;

            m = obj.context.mesh.allElementMultiIndices;
            e = m(:, dim) == 1 | m(:, dim) == obj.context.mesh.resolution(dim);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.context.mesh.indexer.multiToLinear(m, 0);
            m(:, dim) = m(:, dim) + 1;
            l2 = obj.context.mesh.indexer.multiToLinear(m, 0);

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
