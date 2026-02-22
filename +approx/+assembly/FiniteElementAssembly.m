classdef FiniteElementAssembly < handle
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
    %   approx.space.FiniteElementSpace

    properties
        Space % Finite element space
    end

    methods
        function obj = FiniteElementAssembly(space)
            % FINITEELEMENTASSEMBLY Constructor for FiniteElementAssembly
            % class.
            %
            %   obj = FiniteElementAssembly(space) creates a finite element
            %   assembly object with the specified finite element space.
            %   The constructor validates that the provided space is a
            %   valid FiniteElementSpace object.

            arguments
                space approx.space.FiniteElementSpace
            end

            obj.Space = space;
        end
    end
    
    methods (Access = protected)
        function G = createZeros(obj, options)
            % CREATEZEROS Create an all-zero sparse matrix.
            %
            %   G = createZeros(obj) creates a sparse zeros matrix as a
            %   placeholder for assembly with square dimensions.
            %
            %   G = createZeros(obj, nCols=nCols) creates a sparse zeros matrix
            %   as a placeholder for assembly with the specified number of
            %   columns @a nCols.
            
            arguments
                obj approx.assembly.FiniteElementAssembly
                options.nCols {mustBePositive, mustBeInteger} = []
            end
            
            nRows = obj.Space.NGlobalDofs;
            if isempty(options.nCols)
                nCols = obj.Space.NGlobalDofs;
            else
                nCols = options.nCols;
            end
            
            G = sparse(nRows, nCols);
        end
        
        function G = createMatrix(obj, L, options)
            % CREATEMATRIX Create global sparse matrix from local matrices.
            %
            %   G = createMatrix(obj, L) assembles local matrices @a L to a
            %   global sparse matrix @a G using all elements and DoFs.
            %
            %   G = createMatrix(obj, L, TEI=TEI, TLI=TLI, SEI=SEI, SLI=SLI)
            %   assembles with specified elements (@a TEI, @a SEI) and local
            %   DoFs (@a TLI, @a SLI) for test space and trial space.
            
            arguments
                obj approx.assembly.FiniteElementAssembly
                L {mustBeNumeric}
                options.TEI {mustBeInteger} = []
                options.TLI {mustBeInteger} = []
                options.SEI {mustBeInteger} = []
                options.SLI {mustBeInteger} = []
            end
            
            if isempty(options.TEI) && ~isempty(options.SEI)
                TEI = options.SEI;
                SEI = options.SEI;
            elseif ~isempty(options.TEI) && isempty(options.SEI)
                TEI = options.TEI;
                SEI = options.TEI;
            elseif isempty(options.TEI) && isempty(options.SEI)
                TEI = 1:obj.Space.NMeshElements;
                SEI = TEI;
            else
                TEI = options.TEI;
                SEI = options.SEI;
                core.except.assert(length(TEI) == length(SEI), ...
                    'InvalidInput', ...
                    'Incompatible row and column element indices.');
            end
            ne = length(TEI);
            
            nl = obj.Space.NLocalDofs;
            if isempty(options.TLI)
                TLI = 1:nl;
            else
                TLI = options.TLI;
            end
            if isempty(options.SLI)
                SLI = 1:nl;
            else
                SLI = options.SLI;
            end
            
            if ismatrix(L)
                [I, J, V] = find(L(TLI, SLI));
                nnz = numel(V);
                T = zeros(nnz * ne, 3);
                if nnz == 0
                    G = obj.createZeros();
                    return;
                end
                T(:, 1) = obj.Space.localToGlobal(TEI, I);
                T(:, 2) = obj.Space.localToGlobal(SEI, J);
                T(:, 3) = repmat(V(:), ne, 1);
            else
                L = reshape(L, nl, nl, []);
                core.except.assert(size(L, 3) == ne, ...
                    'InvalidInput', 'Invalid local matrices.');
                
                ni = length(TLI);
                nj = length(SLI);
                T = zeros(ni*nj*ne, 3);
                
                [TLI, SLI] = ndgrid(TLI, SLI);
                TLI = TLI(:);
                SLI = SLI(:);
                T(:, 1) = obj.Space.localToGlobal(TEI, TLI);
                T(:, 2) = obj.Space.localToGlobal(SEI, SLI);
                L = reshape(L, nl^2, []);
                KL = TLI + (SLI - 1) * nl;
                
                L = L(KL, :);
                T(:, 3) = L(:);
            end
            
            ng = obj.Space.NGlobalDofs;
            G = core.linalg.sparseFromTriplets(T, ng, ng);
        end
        
        % function G = createVector(obj, L, options)
        %     % CREATEVECTOR Create global sparse vector from local vectors.
        %     %
        %     %   G = createVector(obj, L) assembles local vectors @a L to a
        %     %   global sparse vector @a G using all elements and DoFs.
        %     %
        %     %   G = createVector(obj, L, TEI=TEI, TLI=TLI) assembles with
        %     %   specified elements @a TEI and local DoFs @a TLI.
        % 
        %     arguments
        %         obj approx.assembly.FiniteElementAssembly
        %         L {mustBeNumeric}
        %         options.TEI {mustBeInteger} = []
        %         options.TLI {mustBeInteger} = []
        %     end
        % 
        %     if isempty(options.TEI)
        %         TEI = 1:obj.Space.NMeshElements;
        %     else
        %         TEI = options.TEI;
        %     end
        %     ne = length(TEI);
        % 
        %     nl = obj.Space.NLocalDofs;
        %     if isempty(options.TLI)
        %         TLI = 1:nl;
        %     else
        %         TLI = options.TLI;
        %     end
        % 
        %     if size(L, 2) == 1
        %         L = reshape(L, nl, []);
        %         L = repmat(L, 1, ne);
        %     else
        %         L = reshape(L, nl, ne);
        %         core.except.assert(size(L, 2) == ne, ...
        %             'InvalidInput', 'Invalid local vectors.');
        %     end
        %     L = L(TLI, :);
        %     ni = length(TLI);
        % 
        %     T = zeros(ni*ne, 3);
        %     T(:, 1) = obj.Space.localToGlobal(TEI, []);
        %     T(:, 2) = 1;
        %     T(:, 3) = L(:);
        % 
        %     ng = obj.Space.NGlobalDofs;
        %     G = core.linalg.sparseFromTriplets(T, ng, 1);
        % end
        
        function G = createVector(obj, L, options)
            % CREATEVECTOR Create global sparse matrix from a batch of
            % local vectors.
            %
            %   G = createVector(obj, L) assembles batched local vectors
            %   @a L to a global sparse matrix @a G using all elements and DoFs.
            %
            %   G = createVector(obj, L, TEI=TEI, TLI=TLI) assembles with
            %   specified elements @a TEI and local DoFs @a TLI. Each column
            %   represents a different vector in the batch.
            
            arguments
                obj approx.assembly.FiniteElementAssembly
                L {mustBeNumeric}
                options.TEI {mustBeInteger} = []
                options.TLI {mustBeInteger} = []
            end
            
            if isempty(options.TEI)
                TEI = 1:obj.Space.NMeshElements;
            else
                TEI = options.TEI;
            end
            ne = length(TEI);
            
            nl = obj.Space.NLocalDofs;
            if isempty(options.TLI)
                TLI = 1:nl;
            else
                TLI = options.TLI;
            end
            
            % if ismatrix(L)
            %     L = reshape(L, nl, []);
            %     L = repmat(reshape(L, nl, 1, []), 1, ne, 1);
            % else
            %     L = reshape(L, nl, ne, []);
            %     core.except.assert(size(L, 2) == ne, ...
            %         'InvalidInput', 'Invalid local vectors.');
            % end
            L = reshape(L, nl, ne, []);
            L = L(TLI, :, :);
            ni = length(TLI);
            nc = size(L, 3);
            
            T = zeros(ni*ne*nc, 3);
            T(:, 1) = repmat(obj.Space.localToGlobal(TEI, []), nc, 1);
            T(:, 2) = reshape(repmat(1:nc, ni*ne, 1), [], 1);
            T(:, 3) = L(:);
            
            ng = obj.Space.NGlobalDofs;
            G = core.linalg.sparseFromTriplets(T, ng, nc);
        end
    end
end