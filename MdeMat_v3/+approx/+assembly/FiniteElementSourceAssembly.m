classdef FiniteElementSourceAssembly < approx.assembly.FiniteElementAssembly
    %FINITEELEMENTSOURCEASSEMBLY

    methods
        function G = assemblePointwiseVector(obj, func, varargin)
            % ASSEMBLEPOINTWISEVECTOR Assemble vector contribution from
            % volume projections.

            nl = obj.space.nLocalDofs;
            V = obj.space.element.volume.values;
            W = diag(obj.space.element.volume.weights);
            nq = obj.space.element.volume.nPoints;
            xRef = obj.space.element.volume.nodes;
            U = obj.space.functionEvaluate([], xRef, func, varargin{:});
            U = reshape(U, nq, []);
            L = V * W * U;
            L = obj.space.element.approximator.project(L);
            L = reshape(L, nl, []);
            G = obj.createVector([], [], L);
        end
    end
end
