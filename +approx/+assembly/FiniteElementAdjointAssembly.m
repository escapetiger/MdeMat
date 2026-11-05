classdef FiniteElementAdjointAssembly < approx.assembly.FiniteElementAssembly
    % FINITEELEMENTADJOINTASSEMBLY Assembly for finite element adjoint operator pairs.
    %
    %   FiniteElementAdjointAssembly implements assembly operations for
    %   finite element adjoint operator pairs. These operators arise in
    %   discontinuous Galerkin methods where adjoint consistency is
    %   required for optimal convergence rates and stability.
    %
    %   The class handles both periodic and Dirichlet boundary conditions
    %   with various penalty formulations. Penalty types are specified as
    %   a 2x2 cell array where the first row controls penalty application
    %   (boundary vs all elements, weighted vs unweighted) and the second
    %   row controls penalty direction (central, right, left bias). Elements
    %   can be empty to disable specific penalty terms.
    %
    %   The assembly process combines volume integrals from the divergence
    %   theorem with flux integrals across element interfaces, incorporating
    %   appropriate penalty terms to ensure weak enforcement of boundary
    %   conditions and inter-element continuity.
    %
    % See also:
    %   FiniteElementAssembly, FiniteElementEllipticAssembly

    properties (Constant)
        VALID_BC_TYPES = {'periodic', 'dirichlet'} % Valid boundary condition types
        VALID_PENALTY_TYPES = { ...  % Valid penalty type combinations
            {'boundary', 'boundaryw', 'all', 'allw'}, {'central', 'right', 'left'}; ...
            {'central', 'right', 'left'}, {}}
    end

    properties
        bcType % Boundary condition type
        penaltyType % Penalty type specification (2x2 cell array)
    end

    methods
        function obj = FiniteElementAdjointAssembly(space, bcType, penaltyType)
            % FINITEELEMENTADJOINTASSEMBLY Constructor for adjoint assembly.
            %
            %   Creates an adjoint assembly object with specified boundary
            %   conditions and penalty formulation.
            %
            % Inputs:
            %   space - FiniteElementSpace object
            %   bcType - Boundary condition type
            %   penaltyType - Penalty type specification
            %
            % Outputs:
            %   obj - Constructed FiniteElementAdjointAssembly object

            obj@approx.assembly.FiniteElementAssembly(space);
            obj.setBcType(bcType);
            obj.setPenaltyType(penaltyType);
        end

        function obj = setBcType(obj, bcType)
            % SETBCTYPE Set the boundary condition type.
            %
            %   Validates and sets the boundary condition type for the
            %   assembly operation.
            %
            % Inputs:
            %   obj - The FiniteElementAdjointAssembly object
            %   bcType - Boundary condition type
            %
            % Outputs:
            %   obj - The FiniteElementAdjointAssembly object

            obj.checkValidBcType(bcType);
            obj.bcType = bcType;
        end

        function obj = setPenaltyType(obj, penaltyType)
            % SETPENALTYTYPE Set the penalty type specification.
            %
            %   Validates and sets the penalty formulation for flux
            %   terms in the discontinuous Galerkin method.
            %
            % Inputs:
            %   obj - The FiniteElementAdjointAssembly object
            %   penaltyType - Penalty type specification (2x2 cell array)
            %
            % Outputs:
            %   obj - The FiniteElementAdjointAssembly object

            obj.checkValidPenaltyType(penaltyType);
            obj.penaltyType = penaltyType;
        end

        function G = assembleMatrix(obj)
            % ASSEMBLEMATRIX Assemble the adjoint operator matrix.
            %
            %   Constructs the global adjoint operator matrix by combining
            %   volume and flux contributions with appropriate penalty terms.
            %
            % Outputs:
            %   G - Cell array of assembled matrices (2x2 cell, each containing nd cells)

            nd = obj.space.nDims;

            GV = obj.assembleVolumeMatrix();
            GF = obj.assembleFluxMatrix();
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);
            for d = 1:nd
                G{1, 2}{d} = GF{1, 2}{d} - GV{1, 2}{d};
                G{2, 1}{d} = GF{2, 1}{d} - GV{2, 1}{d};
                G{1, 1}{d} = GF{1, 1}{d};
            end
        end

        function G = assembleVector(obj, func, varargin)
            % ASSEMBLEVECTOR Assemble the adjoint operator vector.
            %
            %   Constructs the right-hand side vector for adjoint problems
            %   with Dirichlet boundary conditions.
            %
            % Inputs:
            %   obj - The FiniteElementAdjointAssembly object
            %   func - Boundary function handle
            %   varargin - Additional arguments for function evaluation
            %
            % Outputs:
            %   G - Assembled vector (ng x 1 sparse)

            if strcmpi(obj.bcType, 'periodic')
                G = obj.createZeros(1);
                return;
            end

            G = obj.assembleFluxVector(func, varargin{:});
        end
    
        function G = assemblePointwiseMatrix(obj, field)
            nd = obj.space.nDims;

            GV = obj.assemblePointwiseVolumeMatrix(field);
            GF = obj.assemblePointwiseFluxMatrix(field);
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);
            for d = 1:nd
                G{1, 2}{d} = GF{1, 2}{d} - GV{1, 2}{d};
                G{2, 1}{d} = GF{2, 1}{d} - GV{2, 1}{d};
                G{1, 1}{d} = GF{1, 1}{d};
            end
        end
    end

    methods (Access = protected)
        function obj = checkValidBcType(obj, bcType)
            % CHECKVALIDBCTYPE Validate boundary condition type.
            %
            %   Ensures the boundary condition type is supported by the
            %   adjoint assembly implementation.
            %
            % Inputs:
            %   obj - The FiniteElementAdjointAssembly object
            %   bcType - Boundary condition type to validate
            %
            % Outputs:
            %   obj - The FiniteElementAdjointAssembly object

            tf = ismember(bcType, obj.VALID_BC_TYPES);
            core.except.assert(tf, 'InvalidInput', ...
                'BC type %s is not supported.', bcType);
        end

        function obj = checkValidPenaltyType(obj, penaltyType)
            % CHECKVALIDPENALTYTYPE Validate penalty type specification.
            %
            %   Ensures the penalty type specification is complete and
            %   contains only supported penalty formulations.
            %
            % Inputs:
            %   obj - The FiniteElementAdjointAssembly object
            %   penaltyType - Penalty type specification to validate
            %
            % Outputs:
            %   obj - The FiniteElementAdjointAssembly object

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
            % COMPUTEINTERNALFLUXPENALTY Compute penalty parameters for internal fluxes.
            %
            %   Calculates penalty coefficients for internal flux terms
            %   based on the specified penalty formulation and face orientation.
            %
            % Inputs:
            %   obj - The FiniteElementAdjointAssembly object
            %   side - Face orientation indicator
            %   ratio - Face scaling factor
            %
            % Outputs:
            %   P - Penalty parameter structure (2x2 cell)

            PT = obj.penaltyType;

            switch PT{1, 1}
                case 'all'
                    beta = 1;
                case 'allw'
                    beta = ratio;
                otherwise
                    beta = 0;
            end
            eta = {beta, 1; 1, 0};
            P = obj.computePenalty(side, eta);
        end

        function P = computeDirichletBoundaryFluxPenalty(obj, side, ratio)
            % COMPUTEDIRICHLETBOUNDARYFLUXPENALTY Compute penalty
            % parameters for Dirichlet boundaries.
            %
            %   P = computeDirichletBoundaryFluxPenalty(obj, side, ratio)
            %   calculates penalty coefficients for Dirichlet boundary flux
            %   terms based on the specified penalty formulation.
            %
            % Inputs:
            %   obj - The FiniteElementAdjointAssembly object
            %   side - Face orientation indicator
            %   ratio - Face scaling factor
            %
            % Outputs:
            %   P - Penalty parameter structure (2x2 cell)

            PT = obj.penaltyType;

            switch PT{1, 1}
                case 'boundary'
                    switch PT{1, 2}
                        case 'central'
                            beta = 1 / 2;
                        case 'right'
                            beta = (side > 0);
                        case 'left'
                            beta = (side < 0);
                        otherwise
                            beta = 0;
                    end
                case 'boundaryw'
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
                case 'all'
                    switch PT{1, 2}
                        case 'central'
                            beta = 1 / 2;
                        case 'right'
                            beta = (side > 0);
                        case 'left'
                            beta = (side < 0);
                        otherwise
                            beta = 0;
                    end
                case 'allw'
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
            % COMPUTEPENALTY Compute penalty coefficients from base
            % parameters.
            %
            %   P = computePanalty(obj, side, eta) transforms base penalty
            %   parameters into the final penalty coefficients used in the
            %   assembly process.
            %
            % Inputs:
            %   obj - The FiniteElementAdjointAssembly object
            %   side - Face orientation indicator
            %   eta - Base penalty parameters (2x2 cell)
            %
            % Outputs:
            %   P - Final penalty coefficients (2x2 cell)

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
                case {'boundary', 'all', 'boundaryw', 'allw'}
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

    methods (Access = private)
        function G = assemblePointwiseVolumeMatrix(obj, field)
            nd = obj.space.nDims;
            ne = obj.space.nMeshElements;
            nq = obj.space.element.volume.nPoints;
            invJ = obj.space.mesh.computeElementInverseJacobians(); % d x d x ne
            xRef = obj.space.element.volume.nodes;
            A = obj.space.functionEvaluate([], xRef, field);
            A = reshape(A, nq, []);
            W0 = obj.space.element.volume.weights;
            W0 = W0(:) .* A;
            W = zeros(nq, nq, ne);
            [i, j] = ndgrid(1:nq, 1:ne);
            I = sub2ind([nq nq ne], i, i, j);
            W(I) = W0(:);
            V = obj.space.element.volume.derivatives{1}; % nl x nq x d
            U = obj.space.element.volume.values; % nl x nq
            R = pagemtimes(pagemtimes(V, W), U.'); % nl x nl x d
            R = obj.space.element.approximator.project(R);
            L = reshape(R, size(V,1), size(U,1), nd, 1, ne) .* reshape(invJ, 1, 1, nd, nd, ne); % nl^2 x d x ne
            L = sum(L, 3);
            L = reshape(L, size(V, 1), size(U, 1), nd, []);
            for d = 1:nd
                Ld = squeeze(L(:, :, d, :));
                G{1, 2}{d} = obj.createMatrix([], [], [], [], Ld);
                G{2, 1}{d} = G{1, 2}{d};
            end
        end

        function G = assemblePointwiseFluxMatrix(obj, field)
            % ASSEMBLEFLUXMATRIX Assemble flux contribution to adjoint
            % matrix.
            
            nd = obj.space.nDims;
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);

            IEI = obj.space.mesh.getInternalElements();
            for IFE = 1:obj.space.element.nFluxes
                GF = obj.assemblePointwiseInternalFluxMatrix(IEI, IFE, field);
                for d = 1:nd
                    G{1, 2}{d} = G{1, 2}{d} + GF{1, 2}{d};
                    G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                    G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                end
            end

            IEB = obj.space.mesh.getBoundaryElements();
            for IFE = 1:obj.space.element.nFluxes
                IE = obj.space.mesh.findInternalElements(IEB, IFE);
                if ~isempty(IE)
                    GF = obj.assemblePointwiseInternalFluxMatrix(IE, IFE, field);
                    for d = 1:nd
                        G{1, 2}{d} = G{1, 2}{d} + GF{1, 2}{d};
                        G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                        G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                    end
                end

                IE = obj.space.mesh.findBoundaryElements(IEB, IFE);
                if ~isempty(IE)
                    switch obj.bcType
                        case 'periodic'
                            GF = obj.assemblePointwiseInternalFluxMatrix(IE, IFE, field);
                        case 'dirichlet'
                            % TBA
                    end
                    for d = 1:nd
                        G{1, 2}{d} = G{1, 2}{d} + GF{1, 2}{d};
                        G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                        G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                    end
                end
            end
        end
        
        function G = assemblePointwiseInternalFluxMatrix(obj, IE, IFE, field)
            % ASSEMBLEINTERNALFLUXMATRIX Assemble internal flux
            % contributions.
            
            nd = obj.space.nDims;
            nl = obj.space.nLocalDofs;
            ne = length(IE);
            nq = obj.space.element.flux(IFE).nPoints;
            normal = obj.space.mesh.computeOutwardNormals(IE, IFE);
            ratio = obj.space.mesh.computeFaceScalingFactors(IE, IFE);
            ratio = reshape(ratio, 1, []);
            direction = obj.space.mesh.computeDirection();
            side = sign(sum(direction.*normal, 1));
            P = obj.computeInternalFluxPenalty(side, ratio);
            FC = ratio .* normal;
            V = obj.space.element.flux(IFE).values;
            xRef = obj.space.element.flux(IFE).nodes;
            A = obj.space.functionEvaluate(IE, xRef, field);
            A = reshape(A, nq, []);
            W0 = obj.space.element.flux(IFE).weights;
            W0 = W0(:) .* A;
            W = zeros(nq, nq, ne);
            [i, j] = ndgrid(1:nq, 1:ne);
            I = sub2ind([nq nq ne], i, i, j);
            W(I) = W0(:);

            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);

            U = obj.space.element.flux(IFE).values;
            L = pagemtimes(pagemtimes(V, W), U.');
            L = obj.space.element.approximator.project(L);
            L = reshape(L, nl, nl, []);
            for d = 1:nd
                C = reshape(FC(d, :), 1, 1, []);
                A = obj.createMatrix(IE, [], IE, [], L.*C);
                G{1, 2}{d} = G{1, 2}{d} + A / 2;
                C = reshape(FC(d, :).*P{1, 2}, 1, 1, []);
                A = obj.createMatrix(IE, [], IE, [], L.*C);
                G{1, 2}{d} = G{1, 2}{d} - A;
                C = reshape(FC(d, :), 1, 1, []);
                A = obj.createMatrix(IE, [], IE, [], L.*C);
                G{2, 1}{d} = G{2, 1}{d} + A / 2;
                C = reshape(FC(d, :).*P{2, 1}, 1, 1, []);
                A = obj.createMatrix(IE, [], IE, [], L.*C);
                G{2, 1}{d} = G{2, 1}{d} - A;
                C = reshape(FC(d, :).*P{1, 1}, 1, 1, []);
                A = obj.createMatrix(IE, [], IE, [], L.*C);
                G{1, 1}{d} = G{1, 1}{d} - A;
            end

            [JE, JFE] = obj.space.mesh.findNeighborElements(IE, IFE, obj.bcType);
            U = arrayfun(@(j) obj.space.element.flux(j).values, JFE, 'Un', 0);
            U = cat(3, U{:});
            L = pagemtimes(pagemtimes(V, W), U.');
            L = obj.space.element.approximator.project(L);
            for d = 1:nd
                C = reshape(FC(d, :), 1, 1, []);
                A = obj.createMatrix(IE, [], JE, [], L.*C);
                G{1, 2}{d} = G{1, 2}{d} + A / 2;
                C = reshape(FC(d, :).*P{1, 2}, 1, 1, []);
                A = obj.createMatrix(IE, [], JE, [], L.*C);
                G{1, 2}{d} = G{1, 2}{d} + A;
                C = reshape(FC(d, :), 1, 1, []);
                A = obj.createMatrix(IE, [], JE, [], L.*C);
                G{2, 1}{d} = G{2, 1}{d} + A / 2;
                C = reshape(FC(d, :).*P{2, 1}, 1, 1, []);
                A = obj.createMatrix(IE, [], JE, [], L.*C);
                G{2, 1}{d} = G{2, 1}{d} + A;
                C = reshape(FC(d, :).*P{1, 1}, 1, 1, []);
                A = obj.createMatrix(IE, [], JE, [], L.*C);
                G{1, 1}{d} = G{1, 1}{d} + A;
            end
        end

        function G = assembleVolumeMatrix(obj)
            % ASSEMBLEVOLUMEMATRIX Assemble volume contribution to adjoint
            % matrix.
            
            nd = obj.space.nDims;
            invJ = obj.space.mesh.computeElementInverseJacobians();
            W = diag(obj.space.element.volume.weights);
            V = obj.space.element.volume.derivatives{1};
            U = obj.space.element.volume.values;
            R = pagemtimes(pagemtimes(V, W), U.');
            R = obj.space.element.approximator.project(R);
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);
            L = reshape(R, [], nd) * reshape(invJ, nd, []);
            L = reshape(L, size(V, 1), size(U, 1), nd, []);

            for d = 1:nd
                Ld = squeeze(L(:, :, d, :));
                G{1, 2}{d} = obj.createMatrix([], [], [], [], Ld);
                G{2, 1}{d} = G{1, 2}{d};
            end
        end

        function G = assembleFluxMatrix(obj)
            % ASSEMBLEFLUXMATRIX Assemble flux contribution to adjoint
            % matrix.
            
            nd = obj.space.nDims;
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);

            IEI = obj.space.mesh.getInternalElements();
            for IFE = 1:obj.space.element.nFluxes
                GF = obj.assembleInternalFluxMatrix(IEI, IFE);
                for d = 1:nd
                    G{1, 2}{d} = G{1, 2}{d} + GF{1, 2}{d};
                    G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                    G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                end
            end

            IEB = obj.space.mesh.getBoundaryElements();
            for IFE = 1:obj.space.element.nFluxes
                IE = obj.space.mesh.findInternalElements(IEB, IFE);
                if ~isempty(IE)
                    GF = obj.assembleInternalFluxMatrix(IE, IFE);
                    for d = 1:nd
                        G{1, 2}{d} = G{1, 2}{d} + GF{1, 2}{d};
                        G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                        G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                    end
                end

                IE = obj.space.mesh.findBoundaryElements(IEB, IFE);
                if ~isempty(IE)
                    switch obj.bcType
                        case 'periodic'
                            GF = obj.assembleInternalFluxMatrix(IE, IFE);
                        case 'dirichlet'
                            GF = obj.assembleDirichletBoundaryFluxMatrix(IE, IFE);
                    end
                    for d = 1:nd
                        G{1, 2}{d} = G{1, 2}{d} + GF{1, 2}{d};
                        G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                        G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                    end
                end
            end
        end

        function G = assembleInternalFluxMatrix(obj, IE, IFE)
            % ASSEMBLEINTERNALFLUXMATRIX Assemble internal flux
            % contributions.
            
            nd = obj.space.nDims;
            nl = obj.space.nLocalDofs;
            normal = obj.space.mesh.computeOutwardNormals(IE, IFE);
            ratio = obj.space.mesh.computeFaceScalingFactors(IE, IFE);
            ratio = reshape(ratio, 1, []);
            direction = obj.space.mesh.computeDirection();
            side = sign(sum(direction.*normal, 1));
            P = obj.computeInternalFluxPenalty(side, ratio);
            FC = ratio .* normal;
            V = obj.space.element.flux(IFE).values;
            W = diag(obj.space.element.flux(IFE).weights);
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);

            U = obj.space.element.flux(IFE).values;
            L = V * W * U.';
            L = obj.space.element.approximator.project(L);
            L = reshape(L, nl, nl, 1);
            for d = 1:nd
                C = reshape(FC(d, :), 1, 1, []);
                A = obj.createMatrix(IE, [], IE, [], L.*C);
                G{1, 2}{d} = G{1, 2}{d} + A / 2;
                C = reshape(FC(d, :).*P{1, 2}, 1, 1, []);
                A = obj.createMatrix(IE, [], IE, [], L.*C);
                G{1, 2}{d} = G{1, 2}{d} - A;
                C = reshape(FC(d, :), 1, 1, []);
                A = obj.createMatrix(IE, [], IE, [], L.*C);
                G{2, 1}{d} = G{2, 1}{d} + A / 2;
                C = reshape(FC(d, :).*P{2, 1}, 1, 1, []);
                A = obj.createMatrix(IE, [], IE, [], L.*C);
                G{2, 1}{d} = G{2, 1}{d} - A;
                C = reshape(FC(d, :).*P{1, 1}, 1, 1, []);
                A = obj.createMatrix(IE, [], IE, [], L.*C);
                G{1, 1}{d} = G{1, 1}{d} - A;
            end

            [JE, JFE] = obj.space.mesh.findNeighborElements(IE, IFE, obj.bcType);
            U = arrayfun(@(j) obj.space.element.flux(j).values, JFE, 'Un', 0);
            U = cat(3, U{:});
            L = pagemtimes(V*W, U.');
            L = obj.space.element.approximator.project(L);
            for d = 1:nd
                C = reshape(FC(d, :), 1, 1, []);
                A = obj.createMatrix(IE, [], JE, [], L.*C);
                G{1, 2}{d} = G{1, 2}{d} + A / 2;
                C = reshape(FC(d, :).*P{1, 2}, 1, 1, []);
                A = obj.createMatrix(IE, [], JE, [], L.*C);
                G{1, 2}{d} = G{1, 2}{d} + A;
                C = reshape(FC(d, :), 1, 1, []);
                A = obj.createMatrix(IE, [], JE, [], L.*C);
                G{2, 1}{d} = G{2, 1}{d} + A / 2;
                C = reshape(FC(d, :).*P{2, 1}, 1, 1, []);
                A = obj.createMatrix(IE, [], JE, [], L.*C);
                G{2, 1}{d} = G{2, 1}{d} + A;
                C = reshape(FC(d, :).*P{1, 1}, 1, 1, []);
                A = obj.createMatrix(IE, [], JE, [], L.*C);
                G{1, 1}{d} = G{1, 1}{d} + A;
            end
        end

        function G = assembleDirichletBoundaryFluxMatrix(obj, IE, IFE)
            % ASSEMBLEDIRICHLETBOUNDARYFLUXMATRIX Assemble Dirichlet
            % boundary flux contributions.
            
            nd = obj.space.nDims;
            nl = obj.space.nLocalDofs;
            normal = obj.space.mesh.computeOutwardNormals(IE, IFE);
            ratio = obj.space.mesh.computeFaceScalingFactors(IE, IFE);
            ratio = reshape(ratio, 1, []);
            direction = obj.space.mesh.computeDirection();
            side = sign(sum(direction.*normal, 1));
            P = obj.computeDirichletBoundaryFluxPenalty(side, ratio);
            FC = ratio .* normal;
            V = obj.space.element.flux(IFE).values;
            W = diag(obj.space.element.flux(IFE).weights);
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);

            U = obj.space.element.flux(IFE).values;
            L = V * W * U.';
            L = obj.space.element.approximator.project(L);
            L = reshape(L, nl, nl, 1);
            for d = 1:nd
                C = reshape(FC(d, :), 1, 1, []);
                A = obj.createMatrix(IE, [], IE, [], L.*C);
                G{1, 2}{d} = G{1, 2}{d} + A;
                C = reshape(FC(d, :).*P{1, 1}, 1, 1, []);
                A = obj.createMatrix(IE, [], IE, [], L.*C);
                G{1, 1}{d} = G{1, 1}{d} - A;
            end
        end

        function G = assembleFluxVector(obj, func, varargin)
            % ASSEMBLEFLUXVECTOR Assemble flux contribution to adjoint
            % vector.
            
            nd = obj.space.nDims;
            G = repmat({repmat({obj.createZeros(1)}, 1, nd)}, 2, 2);

            IEB = obj.space.mesh.getBoundaryElements();
            for IFE = 1:obj.space.element.nFluxes
                IE = obj.space.mesh.findBoundaryElements(IEB, IFE);
                if ~isempty(IE)
                    GF = obj.assembleDirichletBoundaryFluxVector(IE, IFE, func, varargin{:});
                    for d = 1:nd
                        G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                        G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                    end
                end
            end
        end

        function G = assembleDirichletBoundaryFluxVector(obj, IE, IFE, func, varargin)
            % ASSEMBLEDIRICHLETBOUNDARYFLUXVECTOR Assemble Dirichlet
            % boundary flux vector.
            
            nd = obj.space.nDims;
            nl = obj.space.nLocalDofs;
            normal = obj.space.mesh.computeOutwardNormals(IE, IFE);
            ratio = obj.space.mesh.computeFaceScalingFactors(IE, IFE);
            ratio = reshape(ratio, 1, []);
            direction = obj.space.mesh.computeDirection();
            side = sign(sum(direction.*normal, 1));
            P = obj.computeDirichletBoundaryFluxPenalty(side, ratio);
            FC = ratio .* normal;
            V = obj.space.element.flux(IFE).values;
            W = diag(obj.space.element.flux(IFE).weights);
            nq = obj.space.element.flux(IFE).nPoints;
            xRef = obj.space.element.flux(IFE).nodes;
            U = obj.space.functionEvaluate(IE, xRef, func, varargin{:});
            U = reshape(U, nq, []);

            G = repmat({repmat({obj.createZeros(1)}, 1, nd)}, 2, 2);
            for d = 1:nd
                L = V * W * U;
                L = reshape(L, nl, []);
                L = obj.space.element.approximator.project(L);
                C = reshape(FC(d, :), 1, []);
                G{2, 1}{d} = obj.createVector(IE, [], L.*C);
                C = reshape(FC(d, :).*P{1, 1}, 1, []);
                G{1, 1}{d} = obj.createVector(IE, [], L.*C);
            end
        end
    end
end