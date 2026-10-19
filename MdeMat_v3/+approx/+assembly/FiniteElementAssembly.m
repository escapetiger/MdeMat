classdef FiniteElementAssembly < approx.assembly.Assembly
    % FINITEELEMENTASSEMBLY Base class for finite element assembly
    % operations.
    %
    %   FiniteElementAssembly provides the foundation for assembling
    %   discrete operators from finite element spaces. The assembly process
    %   involves iterating over mesh elements, computing local data, and
    %   accumulating contributions into global sparse data using
    %   appropriate connectivity information from the finite element space.
    %
    %   The class provides essential methods for creating sparse matrices
    %   and vectors from local element contributions. Local matrices and
    %   vectors computed on individual elements are assembled into global
    %   sparse structures using element connectivity information from the
    %   finite element space.
    %
    %   Key assembly operations include matrix assembly from element-wise
    %   contributions, vector assembly from local vectors, and batched
    %   vector assembly for multiple right-hand sides. The class handles
    %   the mapping between local element degrees of freedom and global
    %   system indices.
    %
    % See also:
    %   approx.assembly.Assembly, approx.space.FiniteElementSpace

    methods (Access = protected)
        function obj = checkValidSpace(obj)
            % CHECKVALIDSPACE Validate that the space is a FiniteElementSpace object.
            %
            %   obj = checkValidSpace(obj) verifies that the space property
            %   is a valid FiniteElementSpace object. Throws an error if the
            %   space is not of the expected type.
            %
            % Outputs:
            %   obj - The FiniteElementAssembly object (for method chaining)
            
            tf = isa(obj.space, 'approx.space.FiniteElementSpace');
            core.except.assert(tf, 'InvalidInput', ...
                'FiniteElementAssembly requires a FiniteElementSpace object.');
        end

        function G = createZeros(obj, nCols)
            % CREATEZEROS Create an all-zero sparse matrix.
            %
            %   G = createZeros(obj) creates a sparse zeros matrix as a
            %   placeholder for assembly where nCols is the number of degrees
            %   of freedom in trial space.
            %
            %   G = createZeros(obj, nCols) creates a sparse zeros matrix as
            %   a placeholder for assembly with the specified number of
            %   columns nCols.
            %
            % Inputs:
            %   obj - The FiniteElementAssembly object
            %   nCols - The number of columns (optional, scalar)
            %
            % Outputs:
            %   G - Sparse matrix with zeros

            if nargin < 2 || isempty(nCols)
                nCols = obj.space.nGlobalDofs;
            end

            nRows = obj.space.nGlobalDofs;
            G = sparse(nRows, nCols);
        end

        function G = createMatrix(obj, IE, IL, JE, JL, L)
            % CREATEMATRIX Assemble matrix contribution from local matrices.
            %
            %   G = createMatrix(obj, IE, IL, JE, JL, L) assembles local
            %   matrices L to a global sparse matrix G with specified
            %   elements (IE, JE) and local DoFs (IL, JL) for test space
            %   and trial space, respectively.
            %
            % Inputs:
            %   obj - The FiniteElementAssembly object
            %   IE - Element linear indices for test space (ne x 1 vector or empty)
            %   IL - Local DoF indices for test space (ni x 1 vector or empty)
            %   JE - Element linear indices for trial space (ne x 1 vector or empty)
            %   JL - Local DoF indices for trial space (nj x 1 vector or empty)
            %   L - Local matrices (nl x nl matrix or nl x nl x ne array)
            %
            % Outputs:
            %   G - Global sparse matrix (ng x ng)

            if isempty(IE) && ~isempty(JE)
                IE = JE;
            elseif ~isempty(IE) && isempty(JE)
                JE = IE;
            elseif isempty(IE) && isempty(JE)
                IE = 1:obj.space.nMeshElements;
                JE = IE;
            else
                core.except.assert(length(IE) == length(JE), ...
                    'InvalidInput', ...
                    'Incomptiable row and columen element indices.');
            end
            ne = length(IE);

            nl = obj.space.nLocalDofs;
            if isempty(IL), IL = 1:nl; end
            if isempty(JL), JL = 1:nl; end

            if ismatrix(L)
                [I, J, V] = find(L(IL, JL));
                nnz = numel(V);
                T = zeros(nnz * ne, 3);
                if nnz == 0
                    G = obj.createZeros();
                    return;
                end
                T(:, 1) = obj.space.localToGlobal(IE, I);
                T(:, 2) = obj.space.localToGlobal(JE, J);
                T(:, 3) = repmat(V(:), ne, 1);
            else
                L = reshape(L, nl, nl, []);
                core.except.assert(size(L, 3) == ne, ...
                    'InvalidInput', 'Invalid local matrices.');

                ni = length(IL);
                nj = length(JL);
                T = zeros(ni*nj*ne, 3);

                [IL, JL] = ndgrid(IL, JL);
                IL = IL(:);
                JL = JL(:);
                T(:, 1) = obj.space.localToGlobal(IE, IL);
                T(:, 2) = obj.space.localToGlobal(JE, JL);
                L = reshape(L, nl^2, []);
                KL = IL + (JL - 1) * nl;

                L = L(KL, :);
                T(:, 3) = L(:);
            end

            ng = obj.space.nGlobalDofs;
            G = core.linalg.sparseFromTriplets(T, ng, ng);
        end

        function G = createVector(obj, IE, IL, L)
            % CREATEVECTOR Assemble vector contribution from local vectors.
            %
            %   G = createVector(obj, IE, IL, L) assembles local vectors L
            %   to a global sparse vector G with specified elements IE and
            %   local DoFs IL.
            %
            % Inputs:
            %   obj - The FiniteElementAssembly object
            %   IE - Element linear indices (ne x 1 vector or empty)
            %   IL - Local DoF indices (ni x 1 vector or empty)
            %   L - Local vectors (nl x 1 vector or nl x ne matrix)
            %
            % Outputs:
            %   G - Global sparse vector (ng x 1)

            if isempty(IE), IE = 1:obj.space.nMeshElements; end
            ne = length(IE);

            nl = obj.space.nLocalDofs;
            if isempty(IL), IL = 1:nl; end

            if size(L, 2) == 1
                L = reshape(L, nl, []);
                L = repmat(L, 1, ne);
            else
                L = reshape(L, nl, ne);
                core.except.assert(size(L, 2) == ne, ...
                    'InvalidInput', 'Invalid local vectors.');
            end
            L = L(IL, :);
            ni = length(IL);

            T = zeros(ni*ne, 3);
            T(:, 1) = obj.space.localToGlobal(IE, []);
            T(:, 2) = 1;
            T(:, 3) = L(:);

            ng = obj.space.nGlobalDofs;
            G = core.linalg.sparseFromTriplets(T, ng, 1);
        end

        function G = createVectorBatch(obj, IE, IL, L)
            % CREATEVECTORBATCH Assemble batched vector contribution.
            %
            %   G = createVectorBatch(obj, IE, IL, L) assembles batched
            %   local vectors L to a global sparse matrix G with specified
            %   elements IE and local DoFs IL. Each column represents a
            %   different vector in the batch.
            %
            % Inputs:
            %   obj - The FiniteElementAssembly object
            %   IE - Element linear indices (ne x 1 vector or empty)
            %   IL - Local DoF indices (ni x 1 vector or empty)
            %   L - Batched local vectors (nl x nc matrix or nl x ne x nc array)
            %
            % Outputs:
            %   G - Global sparse matrix (ng x nc)

            if isempty(IE), IE = 1:obj.space.nMeshElements; end
            ne = length(IE);

            nl = obj.space.nLocalDofs;
            if isempty(IL), IL = 1:nl; end

            if ismatrix(L)
                L = reshape(L, nl, []);
                L = repmat(reshape(L, nl, 1, []), 1, ne, 1);
            else
                L = reshape(L, nl, ne, []);
                core.except.assert(size(L, 2) == ne, ...
                    'InvalidInput', 'Invalid local vectors.');
            end
            L = L(IL, :, :);
            ni = length(IL);
            nc = size(L, 3);

            T = zeros(ni*ne*nc, 3);
            T(:, 1) = repmat(obj.space.localToGlobal(IE, []), nc, 1);
            T(:, 2) = reshape(repmat(1:nc, ni*ne, 1), [], 1);
            T(:, 3) = L(:);

            ng = obj.space.nGlobalDofs;
            G = core.linalg.sparseFromTriplets(T, ng, nc);
        end
    end
end