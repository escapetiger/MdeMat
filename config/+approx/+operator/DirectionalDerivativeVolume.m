classdef DirectionalDerivativeVolume < approx.operator.MeshOperator    
    methods
        function D = addBilinear(obj, coe)
            d = obj.space.nDims;
            if ~isscalar(coe)
                C = reshape(coe, 1, 1, d, []);
            else
                C = coe;
            end

            invJ = obj.space.mesh.computeElementInverseJacobians();
            W = diag(obj.space.element.volume.weights);
            U = obj.space.element.volume.derivatives(:, :, 1:d);
            V = obj.space.element.volume.values;
            A = pagemtimes(pagemtimes(U, W), V.');
            A = reshape(A, [], d) * reshape(invJ, d, []);
            A = reshape(A, size(U, 1), size(V, 1), d, size(invJ, 3));
            A = A .* C;
            D = obj.addBilinear([], [], [], [], A);
        end
    end
end

