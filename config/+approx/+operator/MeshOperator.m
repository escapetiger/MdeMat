classdef MeshOperator < handle
    % MESHOPERATOR Base class for operators between mesh spaces.

    properties
        space % MeshSpace object
    end

    methods
        function obj = MeshOperator(space)
            % OPERATOR Constructor for Operator.
            %
            %   obj = MeshOperator(space) creates an operator with the
            %   specified mesh space.
            %
            % Inputs:
            %   space - MeshSpace object
            %
            % Outputs:
            %   obj - Constructed Operator object

            core.except.assert(isa(space, 'approx.space.MeshSpace'), ...
                'InvalidInput', 'BilinearForm requires a mesh space.');
            obj.space = space;
        end

        function D = addBilinear(obj, EI, LI, EJ, LJ, A)
            % ADDBILINEAR Add bilinear contribution.
            %
            %   D = addBilinear(obj, EI, LI, EJ, LJ, A) converts local
            %   matrices @a A to a global block diagonal matrix @a D with
            %   specified elements for test space @a EI and trial space @a
            %   EJ.
            %
            % Inputs:
            %   obj - The MeshOperator object
            %   EI - Element linear indices for test space (ne x 1 vector)
            %   LI - Local DoF indices for test space (ni x 1 vector)
            %   EJ - Element linear indices for trial space (ne x 1 vector)
            %   LJ - Local DoF indices for trial space (nj x 1 vector)
            %   A - Local matrices (nl x nl matrix or nl × nl × ne array)
            %
            % Outputs:
            %   D - Block diagonal matrix (ng x ng sparse)

            if isempty(EI) && ~isempty(EJ)
                EI = EJ;
            elseif ~isempty(EI) && isempty(EJ)
                EJ = EI;
            elseif isempty(EI) && isempty(EJ)
                EI = 1:obj.space.nMeshElements;
                EJ = EI;
            else
                core.except.assert(length(EI) == length(EJ), ...
                    'InvalidInput', ...
                    'Incomptiable row and columen element indices.');
            end
            ne = length(EI);

            nl = obj.space.nLocalDofs;
            if isempty(LI), LI = 1:nl; end
            if isempty(LJ), LJ = 1:nl; end

            if ismatrix(A)
                [I, J, V] = find(A(LI, LJ));
                nnz = numel(V);
                T = zeros(nnz * ne, 3);
                T(:, 1) = obj.space.localToGlobal(EI, I);
                T(:, 2) = obj.space.localToGlobal(EJ, J);
                T(:, 3) = V(:);
            else
                A = reshape(A, nl, nl, []);
                core.except.assert(size(A, 3) == ne, ...
                    'InvalidInput', 'Invalid local matrices.');

                ni = length(LI);
                nj = length(LJ);
                T = zeros(ni*nj*ne, 3);

                [LI, LJ] = ndgrid(LI, LJ);
                LI = LI(:);
                LJ = LJ(:);
                T(:, 1) = obj.space.localToGlobal(EI, LI);
                T(:, 2) = obj.space.localToGlobal(EJ, LJ);

                A = A(LI, LJ, :);
                T(:, 3) = A(:);
            end

            ng = obj.space.nGlobalDofs;
            D = core.linalg.sparseFromTriplets(T, ng, ng);
        end

        function D = addLinear(obj, E, L, A)
            % ADDLINEAR Add linear contribution.
            %
            %   D = addLinear(obj, E, L, A) converts local vectors @a A to
            %   a global block diagonal matrix @a D with specified elements
            %   for test space @a E.
            %
            % Inputs:
            %   obj - The MeshOperator object
            %   E - Element linear indices (ne x 1 vector)
            %   L - Local DoF indices (ni x 1 vector)
            %   A - Local vectors (nl x nc matrix or nl x ne x nc array)
            %
            % Outputs:
            %   D - Block diagonal matrix (ng x nc sparse)

            if isempty(E), E = 1:obj.space.nMeshElements; end
            ne = length(E);

            nl = obj.space.nLocalDofs;
            if isempty(L), L = 1:nl; end

            if nargin < 5 || isempty(C), C = ones(1, ne); end

            if ismatrix(A)
                A = reshape(A, nl, []); % nl × nc
                A = repmat(reshape(A, nl, 1, []), 1, ne, 1); % nl × ne × nc
            else
                A = reshape(A, nl, ne, []); % nl × ne × nc
                core.except.assert(size(A, 2) == ne, ...
                    'InvalidInput', 'Invalid local vectors.');
            end
            A = A(L, :, :); % ni × ne × nc
            ni = length(L);
            nc = size(A, 3);

            T = zeros(ni*ne*nc, 3);
            T(:, 1) = repmat(obj.space.localToGlobal(E, []), nc, 1);
            T(:, 2) = reshape(repmat(1:nc, ni*ne, 1), [], 1);
            C = reshape(repmat(C(:).', ni, 1), ni, ne, 1); % ni × ne × 1
            T(:, 3) = reshape(C.*A, [], 1);

            ng = obj.space.nGlobalDofs;
            D = core.linalg.sparseFromTriplets(T, ng, nc);
        end
    end
end
