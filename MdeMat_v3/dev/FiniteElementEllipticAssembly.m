classdef FiniteElementEllipticAssembly < approx.assembly.FiniteElementAssembly
    % FINITEELEMENTELLIPTICASSEMBLY Assembly for finite element elliptic
    % operators.

    properties (Constant)
        VALID_BC_TYPES = {'periodic', 'dirichlet'}
        VALID_PENALTY_TYPES = { ...
            {'boundary', 'all'}, {'central', 'right', 'left'}; ...
            {'central', 'right', 'left'}, {}}
    end

    properties
        bcType % Boundary condition type
        penaltyType % Penalty type
        cache % Cached data
    end

    methods
        function obj = FiniteElementEllipticAssembly(space, bcType, penaltyType)
            obj@approx.assembly.FiniteElementAssembly(space);
            obj.setBcType(bcType);
            obj.setPenaltyType(penaltyType);
            obj.reset();
        end

        function obj = setBcType(obj, bcType)
            obj.checkValidBcType(bcType);
            obj.bcType = bcType;
        end

        function obj = setPenaltyType(obj, penaltyType)
            obj.checkValidPenaltyType(penaltyType);
            obj.penaltyType = penaltyType;
        end

        function obj = reset(obj)
            obj.cache = struct();
        end

        function [A, b] = forceCompletion(obj, A, b)
            V = obj.space.element.volume.values; % nl x nq
            W = diag(obj.space.element.volume.weights); % nq x nq
            L = V * W * ones(length(W), 1); % nl x 1
            c = repmat(L, obj.space.nMeshElements, 1);
            A = [A, c; c.', 0];
            b = [b; 0];
        end

        function G = assembleConstantMatrix(obj, field)
            %< Sizes
            nd = obj.space.nDims;

            %< Diffusion field
            field = reshape(field, nd, nd, []);

            %< Assembly
            GV = obj.assembleConstantVolumeMatrix(field);
            GF = obj.assembleConstantFluxMatrix(field);
            G = obj.zerosLike();
            if strcmpi(obj.bcType, 'dirichlet')
                obj.cache.G12 = cell(1, nd);
            end
            for d = 1:nd
                G12 = GF{1, 2}{d} - GV{1, 2}{d};
                G21 = GF{2, 1}{d} - GV{2, 1}{d};
                G11 = GF{1, 1}{d};
                G = G + G12 * G21 - G11;
                if strcmpi(obj.bcType, 'dirichlet')
                    obj.cache.G12{d} = G12;
                end
            end
        end

        function G = assembleConstantVolumeMatrix(obj, field)
            % ASSEMBLECONSTANTVOLUMEMATRIX Assemble the matrix contribution
            % from volume integrals for the divergence of piecewise
            % constant diffusion field.

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

            %< Initialization
            G = repmat({repmat({obj.zerosLike()}, 1, nd)}, 2, 2);

            %< Physical matrices for (1, 2) component
            L = reshape(R, [], nd) * reshape(invJ, nd, []); % nl^2 x nd x ne
            L = reshape(L, size(V, 1), size(U, 1), nd, []); % nl x nl x nd x ne

            for d = 1:nd
                Ld = squeeze(L(:, :, d, :));
                G{1, 2}{d} = obj.assembleMatrix([], [], [], [], Ld);
            end

            %< Physical matrices for (2, 1) component
            field = permute(field, [2, 1, 3]);
            C = reshape(pagemtimes(invJ, field), nd, nd, []); % nd x nd x ne
            L = reshape(R, [], nd) * reshape(C, nd, []); % nl^2 x nd x ne
            L = reshape(L, size(V, 1), size(U, 1), nd, []); % nl x nl x nd x ne

            for d = 1:nd
                Ld = squeeze(L(:, :, d, :));
                G{2, 1}{d} = obj.assembleMatrix([], [], [], [], Ld);
            end
        end

        function G = assembleConstantFluxMatrix(obj, field)
            %< Initialization
            nd = obj.space.nDims;
            G = repmat({repmat({obj.zerosLike()}, 1, nd)}, 2, 2);

            %< Add contribution from internal elements
            IEI = obj.space.mesh.getInternalElements(); % [] x 1
            for IFE = 1:obj.space.element.nFluxes
                GF = obj.assembleConstantInternalFluxMatrix(IEI, IFE, field);
                for d = 1:nd
                    G{1, 2}{d} = G{1, 2}{d} + GF{1, 2}{d};
                    G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                    G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                end
            end

            %< Add contribution from boundary elements
            IEB = obj.space.mesh.getBoundaryElements(); % [] x 1
            for IFE = 1:obj.space.element.nFluxes
                %< Find elements of which IFE-th face is not on boundary
                IE = obj.space.mesh.findInternalElements(IEB, IFE);
                if ~isempty(IE)
                    GF = obj.assembleConstantInternalFluxMatrix(IE, IFE, field);
                    for d = 1:nd
                        G{1, 2}{d} = G{1, 2}{d} + GF{1, 2}{d};
                        G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                        G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                    end
                end

                %< Find elements of which IFE-th face is on boundary
                IE = obj.space.mesh.findBoundaryElements(IEB, IFE);
                if ~isempty(IE)
                    switch obj.bcType
                        case 'periodic'
                            GF = obj.assembleConstantInternalFluxMatrix(IE, IFE, field);
                        case 'dirichlet'
                            GF = obj.assembleConstantDirichletBoundaryFluxMatrix(IE, IFE, field);
                    end
                    for d = 1:nd
                        G{1, 2}{d} = G{1, 2}{d} + GF{1, 2}{d};
                        G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                        G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                    end
                end
            end
        end

        function G = assembleConstantInternalFluxMatrix(obj, IE, IFE, field)
            %< Sizes
            nd = obj.space.nDims;
            nl = obj.space.nLocalDofs;

            %< Outward normal vectors
            normal = obj.space.mesh.computeOutwardNormals(IE, IFE); % nd x []

            %< Face scaling factors
            ratio = obj.space.mesh.computeFaceScalingFactors(IE, IFE); % [] x 1
            ratio = reshape(ratio, 1, []);

            %< Determine side
            direction = obj.space.mesh.computeDirection(); % nd x 1
            side = sign(sum(direction.*normal, 1)); % 1 x []

            %< Penalty parameters
            P = obj.computeInternalFluxPenalty(side, ratio);

            %< Face coefficients
            field = reshape(field, nd, nd, []);
            FC1 = ratio .* normal; % nd x []
            normal = reshape(normal, nd, 1, []); % nd x 1 x []
            FC2 = ratio .* reshape(pagemtimes(field, normal), nd, []);

            %< Test basis values and weights
            V = obj.space.element.flux(IFE).values; % nl x nq
            W = diag(obj.space.element.flux(IFE).weights); % nq x nq

            %< Initialization
            G = repmat({repmat({obj.zerosLike()}, 1, nd)}, 2, 2);

            %< Interior contribution
            U = obj.space.element.flux(IFE).values; % nl x nq
            L = V * W * U.'; % nl x nl
            L = obj.space.element.approximator.project(L); % nl x nl
            L = reshape(L, nl, nl, 1);
            for d = 1:nd
                %< (1, 2) component
                C = reshape(FC1(d, :), 1, 1, []);
                A = obj.assembleMatrix(IE, [], IE, [], L.*C);
                G{1, 2}{d} = G{1, 2}{d} + A / 2;
                C = reshape(FC1(d, :).*P{1, 2}, 1, 1, []);
                A = obj.assembleMatrix(IE, [], IE, [], L.*C);
                G{1, 2}{d} = G{1, 2}{d} - A;
                %< (2, 1) component
                C = reshape(FC2(d, :), 1, 1, []);
                A = obj.assembleMatrix(IE, [], IE, [], L.*C);
                G{2, 1}{d} = G{2, 1}{d} + A / 2;
                C = reshape(FC2(d, :).*P{2, 1}, 1, 1, []);
                A = obj.assembleMatrix(IE, [], IE, [], L.*C);
                G{2, 1}{d} = G{2, 1}{d} - A;
                %< (1, 1) component
                C = reshape(FC1(d, :).*P{1, 1}, 1, 1, []);
                A = obj.assembleMatrix(IE, [], IE, [], L.*C);
                G{1, 1}{d} = G{1, 1}{d} - A;
            end

            %< Exterior contribution
            [JE, JFE] = obj.space.mesh.findNeighborElements(IE, IFE, obj.bcType); % ne x 1
            U = arrayfun(@(j) obj.space.element.flux(j).values, JFE, 'Un', 0);
            U = cat(3, U{:}); % nl x nq x ne
            L = pagemtimes(V*W, U.'); % nl x nl x ne
            L = obj.space.element.approximator.project(L); % nl x nl x ne
            for d = 1:nd
                %< (1, 2) component
                C = reshape(FC1(d, :), 1, 1, []);
                A = obj.assembleMatrix(IE, [], JE, [], L.*C);
                G{1, 2}{d} = G{1, 2}{d} + A / 2;
                C = reshape(FC1(d, :).*P{1, 2}, 1, 1, []);
                A = obj.assembleMatrix(IE, [], JE, [], L.*C);
                G{1, 2}{d} = G{1, 2}{d} + A;
                %< (2, 1) component
                C = reshape(FC2(d, :), 1, 1, []);
                A = obj.assembleMatrix(IE, [], JE, [], L.*C);
                G{2, 1}{d} = G{2, 1}{d} + A / 2;
                C = reshape(FC2(d, :).*P{2, 1}, 1, 1, []);
                A = obj.assembleMatrix(IE, [], JE, [], L.*C);
                G{2, 1}{d} = G{2, 1}{d} + A;
                %< (1, 1) component
                C = reshape(FC1(d, :).*P{1, 1}, 1, 1, []);
                A = obj.assembleMatrix(IE, [], JE, [], L.*C);
                G{1, 1}{d} = G{1, 1}{d} + A;
            end
        end

        function G = assembleConstantDirichletBoundaryFluxMatrix(obj, IE, IFE, field)
            %< Sizes
            nd = obj.space.nDims;
            nl = obj.space.nLocalDofs;

            %< Outward normal vectors
            normal = obj.space.mesh.computeOutwardNormals(IE, IFE); % nd x []

            %< Face scaling factors
            ratio = obj.space.mesh.computeFaceScalingFactors(IE, IFE); % [] x 1
            ratio = reshape(ratio, 1, []);

            %< Determine side
            direction = obj.space.mesh.computeDirection(); % nd x 1
            side = sign(sum(direction.*normal, 1)); % 1 x []

            %< Penalty parameters
            P = obj.computeDirichletBoundaryFluxPenalty(side, ratio);

            %< Face coefficients
            field = reshape(field, nd, nd, []);
            FC1 = ratio .* normal; % nd x []
            normal = reshape(normal, nd, 1, []); % nd x 1 x []
            FC2 = ratio .* reshape(pagemtimes(field, normal), nd, []);

            %< Test basis values and weights
            V = obj.space.element.flux(IFE).values; % nl x nq
            W = diag(obj.space.element.flux(IFE).weights); % nq x nq

            %< Initialization
            G = repmat({repmat({obj.zerosLike()}, 1, nd)}, 2, 2);

            %< Interior contribution
            U = obj.space.element.flux(IFE).values; % nl x nq
            L = V * W * U.'; % nl x nl
            L = obj.space.element.approximator.project(L); % nl x nl
            L = reshape(L, nl, nl, 1);
            for d = 1:nd
                %< (1, 2) component
                C = reshape(FC1(d, :), 1, 1, []);
                A = obj.assembleMatrix(IE, [], IE, [], L.*C); % average
                G{1, 2}{d} = G{1, 2}{d} + A;
                %< (1, 1) component
                C = reshape(FC2(d, :).*P{1, 1}, 1, 1, []);
                A = obj.assembleMatrix(IE, [], IE, [], L.*C); % jump
                G{1, 1}{d} = G{1, 1}{d} - A;
            end
        end

        function G = assembleConstantVector(obj, field, func, varargin)
            %< Skip for periodic boundary condition
            if strcmpi(obj.bcType, 'periodic')
                G = obj.zerosLike(1);
                return;
            end

            nd = obj.space.nDims;

            field = reshape(field, nd, nd, []);

            GF = obj.assembleConstantFluxVector(field, func, varargin{:});

            G = obj.zerosLike(1);
            for d = 1:nd
                G = G + obj.cache.G12{d} * GF{2, 1}{d} - GF{1, 1}{d};
            end
        end

        function G = assembleConstantFluxVector(obj, field, func, varargin)
            %< Initialization
            nd = obj.space.nDims;
            G = repmat({repmat({obj.zerosLike(1)}, 1, nd)}, 2, 2);

            %< Contribution from boundary elements
            IEB = obj.space.mesh.getBoundaryElements(); % [] x 1
            for IFE = 1:obj.space.element.nFluxes
                %< Find elements of which IFE-th face is on boundary
                IE = obj.space.mesh.findBoundaryElements(IEB, IFE);
                if ~isempty(IE)
                    GF = obj.assembleConstantDirichletBoundaryFluxVector(IE, IFE, field, func, varargin{:});
                    for d = 1:nd
                        G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                        G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                    end
                end
            end
        end

        function G = assembleConstantDirichletBoundaryFluxVector(obj, IE, IFE, field, func, varargin)
            %< Sizes
            nd = obj.space.nDims;
            nl = obj.space.nLocalDofs;

            %< Outward normal vectors
            normal = obj.space.mesh.computeOutwardNormals(IE, IFE); % d x []

            %< Face scaling factors
            ratio = obj.space.mesh.computeFaceScalingFactors(IE, IFE); % [] x 1
            ratio = reshape(ratio, 1, []);

            %< Determine side
            direction = obj.space.mesh.computeDirection(); % nd x 1
            side = sign(sum(direction.*normal, 1)); % 1 x []

            %< Penalty parameters
            P = obj.computeDirichletBoundaryFluxPenalty(side, ratio);

            %< Face coefficients
            field = reshape(field, nd, nd, []);
            normal = reshape(normal, nd, 1, []); % nd x 1 x []
            FC2 = ratio .* reshape(pagemtimes(field, normal), nd, []);

            %< Test basis values and weights
            V = obj.space.element.flux(IFE).values; % nl x nq
            W = diag(obj.space.element.flux(IFE).weights); % nq x nq

            %< Face integration points
            nq = obj.space.element.flux(IFE).nPoints;
            xRef = obj.space.element.flux(IFE).nodes; % d x nq
            U = obj.space.functionEvaluate(IE, xRef, func, varargin{:});
            U = reshape(U, nq, []); % nq x ne

            G = repmat({repmat({obj.zerosLike(1)}, 1, nd)}, 2, 2);
            for d = 1:nd
                L = V * W * U; % nl x ne
                L = reshape(L, nl, []); % nl x ne
                L = obj.space.element.approximator.project(L); % nl x ne
                C = reshape(FC2(d, :), 1, []); % 1 x ne
                G{2, 1}{d} = obj.assembleVector(IE, [], L.*C);
                C = reshape(FC2(d, :).*P{1, 1}, 1, []); % 1 x ne
                G{1, 1}{d} = obj.assembleVector(IE, [], L.*C);
            end
        end

    end

    methods (Access = protected)
        function obj = checkValidBcType(obj, bcType)
            tf = ismember(bcType, obj.VALID_BC_TYPES);
            core.except.assert(tf, 'InvalidInput', ...
                'BC type %s is not supported.', bcType);
        end

        function obj = checkValidPenaltyType(obj, penaltyType)
            core.except.assert(isequal(size(penaltyType), [2, 2]), ...
                'InvalidInput', 'Penalty types are incomplete.');
            for i = 1:numel(penaltyType)
                tf1 = ismember(penaltyType{i}, obj.VALID_PENALTY_TYPES{i});
                tf2 = isempty(penaltyType{i});
                tf = tf1 | tf2;
                core.except.assert(tf, 'InvalidInput', ...
                    'Penalty type %s is not supported.', penaltyType{i});
            end
        end

        function P = computeInternalFluxPenalty(obj, side, ratio)
            PT = obj.penaltyType;

            switch PT{1, 1}
                case 'all'
                    switch PT{1, 2}
                        case 'central'
                            beta = ratio / 2;
                        case 'right'
                            beta = ratio .* (side > 0);
                        case 'left'
                            beta = ratio .* (side < 0);
                        otherwise
                            beta = 0;
                    end
                otherwise
                    beta = 0;
            end
            eta = {beta, 1; 1, 0};
            P = obj.computePenalty(side, eta); % 2 x 2 cell
        end

        function P = computeDirichletBoundaryFluxPenalty(obj, side, ratio)
            PT = obj.penaltyType;

            switch PT{1,1}
                case 'boundary'
                    switch PT{1, 2}
                        case 'central'
                            beta = -ratio / 2;
                        case 'right'
                            beta = -ratio .* (side > 0);
                        case 'left'
                            beta = -ratio .* (side < 0);
                        otherwise
                            beta = 0;
                    end
                case 'all'
                    switch PT{1, 2}
                        case 'central'
                            beta = ratio / 2;
                        case 'right'
                            beta = ratio .* (side > 0);
                        case 'left'
                            beta = ratio .* (side < 0);
                        otherwise
                            beta = 0;
                    end
                otherwise
                    beta = 0;
            end
            eta = {beta, 0; 0, 0};
            P = obj.computePenalty(side, eta);
        end

        function P = computePenalty(obj, side, eta)
            PT = obj.penaltyType;
            P = cell(2, 2);

            switch PT{1, 2}
                case 'right'
                    beta = 1 / 2;
                case 'left'
                    beta = -1 / 2;
                otherwise
                    beta = 0;
            end

            P{1, 2} = side .* eta{1, 2} .* beta;

            switch PT{1, 1}
                case {'boundary', 'all'}
                    alpha = 1;
                otherwise
                    alpha = 0;
            end

            P{1, 1} = side .* eta{1, 1} .* alpha;

            switch PT{2, 1}
                case 'right'
                    beta = 1 / 2;
                case 'left'
                    beta = -1 / 2;
                otherwise
                    beta = 0;
            end

            P{2, 1} = side .* eta{2, 1} .* beta;
        end
    end
end
