classdef FiniteElementUpwindAssembly < approx.assembly.FiniteElementAssembly
    % FINITEELEMENTUPWINDASSEMBLY Assembly for finite element upwind
    % operators.
    %
    %   FiniteElementUpwindAssembly handles the assembly of upwind finite
    %   element operators. These operators are mainly used for describing
    %   advection procedure.

    properties (Constant)
        VALID_BC_TYPES = {'periodic', 'dirichlet'}
    end

    properties
        bcType % Boundary condition type
    end

    methods
        function obj = FiniteElementUpwindAssembly(space, bcType)
            obj@approx.assembly.FiniteElementAssembly(space);
            obj.setBcType(bcType);
        end

        function obj = setBcType(obj, bcType)
            obj.checkValidBcType(bcType);
            obj.bcType = bcType;
        end

        function G = assembleMatrix(obj, field, direction)
            % ASSEMBLEVOLUMEMATRIX Assemble the matrix contribution
            % for the divergence of piecewise constant advection field.

            d = obj.space.nDims;

            %< Advection field
            field = reshape(field, d, []);

            %< Direction
            if nargin < 3 || isempty(direction)
                direction = field;
            else
                direction = reshape(direction, d, []);
            end

            %< Assembly
            GV = obj.assembleVolumeMatrix(field);
            GF = obj.assembleFluxMatrix(field, direction);
            G = GF - GV;
        end

        function G = assembleVector(obj, field, direction, func, varargin)
            % ASSEMBLEVECTOR Assemble the vector contribution for
            % the divergence of piecewise constant advection field.

            %< Skip for periodic boundary condition
            if strcmpi(obj.bcType, 'periodic')
                G = obj.createZeros(1);
                return;
            end

            d = obj.space.nDims;

            %< Advection field
            field = reshape(field, d, []);

            %< Direction
            if nargin < 3 || isempty(direction)
                direction = field;
            else
                direction = reshape(direction, d, []);
            end

            G = obj.assembleFluxVector(field, direction, func, varargin{:});
        end

    end

    methods (Access = protected)
        function obj = checkValidBcType(obj, bcType)
            tf = ismember(bcType, obj.VALID_BC_TYPES);
            core.except.assert(tf, 'InvalidInput', ...
                'BC type %s is not supported.', bcType);
        end
    end

    methods (Access = private)
        function G = assembleVolumeMatrix(obj, field)
            % ASSEMBLEVOLUMEMATRIX Assemble the matrix contribution
            % from volume integrals for the divergence of piecewise
            % constant advection field.

            %< Sizes
            nd = obj.space.nDims;

            %< Inverse of Jacobian
            invJ = obj.space.mesh.computeElementInverseJacobians(); % nd x nd x ne

            %< Integration data
            W = diag(obj.space.element.volume.weights); % nq x nq
            V = obj.space.element.volume.derivatives{1}; % nl x nq x nd
            U = obj.space.element.volume.values; % nl x nq

            %< Reference matrices
            R = pagemtimes(pagemtimes(V, W), U.'); % nl x nl x nd
            R = obj.space.element.approximator.project(R); % nl x nl x nd

            %< Physical marices
            L = reshape(R, [], nd) * reshape(invJ, nd, []); % nl^2 x nd x ne
            L = reshape(L, size(V, 1), size(U, 1), nd, []);
            L = sum(L.*reshape(field, 1, 1, nd, []), 3);
            L = reshape(L, size(V, 1), size(U, 1), []);
            G = obj.createMatrix([], [], [], [], L);
        end

        function G = assembleFluxMatrix(obj, field, direction)
            % ASSEMBLEFLUXMATRIX Assemble the matrix contribution from flux
            % integrals for the divergence of piecewise constant advection
            % field.

            %< Initialization
            G = obj.createZeros();

            %< Add contribution from internal elements
            IEI = obj.space.mesh.getInternalElements(); % [] x 1
            for IFE = 1:obj.space.element.nFluxes
                GF = obj.assembleInternalFluxMatrix(IEI, IFE, field, direction);
                G = G + GF;
            end

            %< Add contribution from boundary elements
            IEB = obj.space.mesh.getBoundaryElements(); % [] x 1
            for IFE = 1:obj.space.element.nFluxes
                %< Find elements of which IFE-th face is not on boundary
                IE = obj.space.mesh.findInternalElements(IEB, IFE);
                if ~isempty(IE)
                    GF = obj.assembleInternalFluxMatrix(IE, IFE, field, direction);
                    G = G + GF;
                end

                %< Find elements of which IFE-th face is on boundary
                IE = obj.space.mesh.findBoundaryElements(IEB, IFE);
                if ~isempty(IE)
                    GF = obj.assembleInternalFluxMatrix(IE, IFE, field, direction);
                    G = G + GF;
                end
            end
        end

        function G = assembleInternalFluxMatrix(obj, IE, IFE, field, direction)
            %< Sizes
            nl = obj.space.nLocalDofs;

            %< Outward normal vectors
            normal = obj.space.mesh.computeOutwardNormals(IE, IFE); % nd x []

            %< Stabilizer
            eta = sum(direction.*normal, 1);
            mask = eta < 0;
            if size(direction, 2) == 1 && size(normal, 2) == 1
                mask = repmat(mask, 1, length(IE));
            end

            %< Test basis values and weights
            V = obj.space.element.flux(IFE).values; % nl x nq
            W = diag(obj.space.element.flux(IFE).weights); % nq x nq

            %< Face scaling factors
            ratio = obj.space.mesh.computeFaceScalingFactors(IE, IFE); % 1 x []
            ratio = reshape(ratio, 1, []);

            %< Face coefficients
            c = ratio .* sum(field.*normal, 1);

            %< Initialization
            G = obj.createZeros();

            %< Interior contribution
            if ~any(~mask)
                A = obj.createZeros();
            else
                U = obj.space.element.flux(IFE).values; % nl x nq
                L = V * W * U.'; % nl x nl
                L = obj.space.element.approximator.project(L); % nl x nl
                L = reshape(L, nl, nl, 1);
                C = reshape(c(~mask), 1, 1, []);
                A = obj.createMatrix(IE(~mask), [], IE(~mask), [], L.*C);
            end
            G = G + A;

            %< Exterior contribution
            if ~any(mask)
                A = obj.createZeros();
            else
                [JE, JFE] = obj.space.mesh.findNeighborElements(IE(mask), IFE, obj.bcType); % ne x 1
                U = arrayfun(@(j) obj.space.element.flux(j).values, JFE, 'Un', 0);
                U = cat(3, U{:}); % nl x nq x ne
                L = pagemtimes(V*W, U.'); % nl x nl x ne
                L = obj.space.element.approximator.project(L); % nl x nl x ne
                C = reshape(c(mask), 1, 1, []);
                A = obj.createMatrix(IE(mask), [], JE, [], L.*C);
            end
            G = G + A;
        end

        %         function G = assemblePointwiseVolume(obj, field)
        %             % ASSEMBLEPOINTWISEVOLUME Assemble the volume integrals with
        %             % spatially varying advection field.
        %
        %             d = obj.space.nDims;
        %             invJ = obj.space.mesh.computeElementInverseJacobians(); % d x d x ne
        %             X = obj.space.element.volume.nodes;
        %             X = obj.space.mesh.collocate(X, []);
        %             W = obj.space.element.volume.weights(:).';
        %             W = W .* reshape(field(X).', [], d);
        %             W = arrayfun(@(i) diag(W(:, i)), 1:d, 'Un', 0);
        %             W = cat(3, W{:}); % nq x nq x d
        %             V = obj.space.element.volume.derivatives{1}; % nl x nq x d
        %             U = obj.space.element.volume.values; % nl x nq
        %             L = pagemtimes(pagemtimes(V, W), U.'); % nl x nl x d
        %             L = reshape(L, [], d) * reshape(invJ, d, []); % nl^2 x d x ne
        %             L = reshape(L, size(V, 1), size(U, 1), d, size(invJ, 3));
        %             G = obj.createMatrix([], [], [], [], L);
        %         end
        %

        function G = assembleFluxVector(obj, field, direction, func, varargin)
            % ASSEMBLEFLUXVECTOR Assemble the vector contribution
            % from flux integrals for the divergence of piecewise constant
            % advection field.

            %< Initialization
            G = obj.createZeros(1);

            %< Contribution from boundary elements
            IEB = obj.space.mesh.getBoundaryElements(); % [] x 1
            for IFE = 1:obj.space.element.nFluxes
                %< Find elements of which IFE-th face is on boundary
                IE = obj.space.mesh.findBoundaryElements(IEB, IFE);
                if ~isempty(IE)
                    GB = obj.assembleDirichletBoundaryFluxVector(IE, IFE, field, direction, func, varargin{:});
                    G = G + GB;
                end
            end
        end

        function G = assembleDirichletBoundaryFluxVector(obj, IE, IFE, field, direction, func, varargin)
            %< Outward normal vectors
            normal = obj.space.mesh.computeOutwardNormals(IE, IFE); % d x []

            %< Stabilizer
            eta = sum(direction.*normal, 1);
            mask = eta < 0;
            if size(direction, 2) == 1 && size(normal, 2) == 1
                mask = repmat(mask, 1, length(IE));
            end

            if ~any(mask)
                G = obj.createZeros(1);
                return;
            end

            %< Face scaling factors
            ratio = obj.space.mesh.computeFaceScalingFactors(IE, IFE); % 1 x []
            ratio = reshape(ratio, 1, []);

            %< Face coefficients
            c = ratio .* sum(field.*normal, 1); % 1 x []

            %< Test basis values and weights
            V = obj.space.element.flux(IFE).values; % nl x nq
            W = diag(obj.space.element.flux(IFE).weights); % nq x nq

            %< Face integration points
            nq = obj.space.element.flux(IFE).nPoints;
            xRef = obj.space.element.flux(IFE).nodes; % d x nq
            U = obj.space.functionEvaluate(IE(mask), xRef, func, varargin{:});
            U = reshape(U, nq, []); % nq x ne

            nl = obj.space.nLocalDofs;
            L = V * W * U; % nl x ne
            L = reshape(L, nl, []); % nl x ne
            L = obj.space.element.approximator.project(L); % nl x ne
            C = reshape(c(mask), 1, []); % 1 x ne
            G = obj.createVector(IE(mask), [], L.*C);
        end
    end
end
