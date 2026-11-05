classdef GridInteriorFluxAssembly < fem.assembly.GridAssembly
    methods
        function T = assembleInnerLeftFlux(obj, d, c)
            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) > 1 & m(:, d) < obj.mesh.resolution(d);
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 0);
            ne = numel(l);

            T = zeros(np^2*ne, 3);

            F = G.fluxData(1, 2*d).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l), v(:));
        end

        function T = assembleInnerRightFlux(obj, d, c)

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) > 1 & m(:, d) < obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 0);

            T = zeros(np^2*ne, 3);

            F = G.fluxData(1, 2*d-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l), v(:));
        end
        
        function T = assembleOuterLeftFlux(obj, d, c)
            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) > 1 & m(:, d) < obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.mesh.indexer.multiToLinear(m, 0);
            m(:, d) = m(:, d) - 1;
            l2 = obj.mesh.indexer.multiToLinear(m, 0);

            T = zeros(np^2*ne, 3);

            F = G.fluxData(2, 2*d-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l1), v(:));
        end

        function T = assembleOuterRightFlux(obj, d, c)
            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) > 1 & m(:, d) < obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.mesh.indexer.multiToLinear(m, 0);
            m(:, d) = m(:, d) + 1;
            l2 = obj.mesh.indexer.multiToLinear(m, 0);

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
