classdef FiniteElementMultiplierAssembly < approx.assembly.FiniteElementAssembly
    % FINITEELEMENTMULTIPLICATIONASSEMBLY Assembly for finite element
    % multiplication operators.

    methods
        function G = assembleConstant(obj, field)
            % ASSEMBLECONSTANT Assemble diagonal matrix with constant
            % field.
            %
            %   G = assemble(obj, field) creates a diagonal sparse matrix
            %   where each element contributes a constant coefficient
            %   repeated for all its local degrees of freedom.
            %
            % Inputs:
            %   obj - The FiniteElementConstantAssembly object
            %   coe - Coefficients
            %
            % Outputs:
            %   G - Sparse matrix

            ng = obj.space.nGlobalDofs;
            nl = obj.space.nLocalDofs;

            if isscalar(field)
               field = repmat(field, obj.space.nMeshElements, 1); 
            end

            G = spdiags(kron(field(:), ones(nl, 1)), 0, ng, ng);
        end

        function G = assemblePointwise(obj, field, varargin)
            nq = obj.space.element.volume.nPoints;
            ne = obj.space.nMeshElements;
            xRef = obj.space.element.volume.nodes;
            A = obj.space.functionEvaluate([], xRef, field, varargin{:});
            A = reshape(A, nq, []);
            W0 = obj.space.element.volume.weights; % 1 x nq
            W0 = W0(:) .* A; % nq x ne
            W = zeros(nq, nq, ne);
            [i, j] = ndgrid(1:nq, 1:ne);
            I = sub2ind([nq nq ne], i, i, j);
            W(I) = W0(:);
            V = obj.space.element.volume.values; % nl x nq
            U = obj.space.element.volume.values; % nl x nq
            L = pagemtimes(pagemtimes(V, W), U.');
            L = obj.space.element.approximator.project(L);
            L = reshape(L, size(V, 1), size(U, 1), []);
            G = obj.createMatrix([], [], [], [], L);
        end
    end
end