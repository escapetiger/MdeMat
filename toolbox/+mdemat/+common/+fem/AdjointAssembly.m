classdef AdjointAssembly < mdemat.common.fem.FiniteElementAssembly
    % FINITEELEMENTADJOINTASSEMBLY Assembly for finite element adjoint
    % operator pairs.
    %
    %   AdjointAssembly implements assembly operations for
    %   finite element adjoint operator pairs. These operators arise in
    %   discontinuous Galerkin methods where adjoint consistency is
    %   required for optimal convergence rates and stability.
    %
    %   The class handles both periodic and Dirichlet boundary conditions
    %   with various penalty formulations. Penalty types are specified as a
    %   2x2 cell array where the first row controls penalty application
    %   (boundary vs all elements, weighted vs unweighted) and the second
    %   row controls penalty direction (central, right, left bias).
    %   Elements can be empty to disable specific penalty terms.
    %
    %   The assembly process combines volume integrals from the divergence
    %   theorem with flux integrals across element interfaces,
    %   incorporating appropriate penalty terms to ensure weak enforcement
    %   of boundary conditions and inter-element continuity.
    %
    % See also:
    %   mdemat.common.fem.FiniteElementAssembly, 
    %   mdemat.common.fem.EllipticAssembly

    properties
        BcType % Boundary condition type
        PenaltyType % Penalty type specification (2x2 cell array)
    end

    methods
        function obj = AdjointAssembly(space, bcType, penaltyType)
            % FINITEELEMENTADJOINTASSEMBLY Constructor for adjoint
            % assembly.
            %
            %   Creates an adjoint assembly object with specified boundary
            %   conditions and penalty formulation.

            arguments
                space mdemat.common.fem.FiniteElementSpace
                bcType{mustBeMember(bcType, {'periodic', 'dirichlet'})}
                penaltyType{mustBePenaltyType(penaltyType)}
            end

            obj@mdemat.common.fem.FiniteElementAssembly(space);
            obj.setBcType(bcType);
            obj.setPenaltyType(penaltyType);
        end

        function obj = setBcType(obj, bcType)
            % SETBCTYPE Set the boundary condition type.
            %
            %   obj = setBcType(obj, bcType) validates and sets the
            %   boundary condition type for the assembly operation.

            arguments
                obj mdemat.common.fem.AdjointAssembly
                bcType{mustBeMember(bcType, {'periodic', 'dirichlet'})}
            end

            obj.BcType = bcType;
        end

        function obj = setPenaltyType(obj, penaltyType)
            % SETPENALTYTYPE Set the penalty type specification.
            %
            %   obj = setPenaltyType(obj, penaltyType) validates and sets
            %   the penalty formulation for flux terms in the discontinuous
            %   Galerkin method.

            arguments
                obj mdemat.common.fem.AdjointAssembly
                penaltyType{mustBePenaltyType(penaltyType)}
            end

            obj.PenaltyType = penaltyType;
        end

        function G = assembleMatrix(obj, field, options)
            % ASSEMBLEMATRIX Assemble the adjoint operator matrix.
            %
            %   G = assembleMatrix(obj, field) constructs the global
            %   adjoint operator matrix by combining volume and flux
            %   contributions with appropriate penalty terms.
            %
            %   G = assembleMatrix(obj, field, args=args) constructs 
            %   the global adjoint operator with spatially varying 
            %   field and @a args passed to @a field.

            arguments
                obj mdemat.common.fem.AdjointAssembly
                field{mustBeA(field, {'numeric', 'function_handle'})}
                options.args = {}
            end

            nd = obj.Space.NDims;

            if isa(field, 'function_handle')
                GV = obj.assemblePointwiseVolumeMatrix(field, options.args{:});
                GF = obj.assemblePointwiseFluxMatrix(field, options.args{:});
            else
                GV = obj.assembleConstantVolumeMatrix(field);
                GF = obj.assembleConstantFluxMatrix(field);
            end

            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);
            for d = 1:nd
                G{1, 2}{d} = GF{1, 2}{d} - GV{1, 2}{d};
                G{2, 1}{d} = GF{2, 1}{d} - GV{2, 1}{d};
                G{1, 1}{d} = GF{1, 1}{d};
            end
        end

        function G = assembleVector(obj, field, func, options)
            % ASSEMBLEVECTOR Assemble the adjoint operator vector.
            %
            %   G = assembleVector(obj, func) constructs the right-hand
            %   side vector for adjoint problems with Dirichlet boundary
            %   conditions using boundary function @a func.
            %
            %   G = assembleVector(obj, func, options) allows additional
            %   arguments to be passed to function evaluation through @a
            %   options.args.

            arguments
                obj mdemat.common.fem.AdjointAssembly
                field{mustBeA(field, {'numeric', 'function_handle'})}
                func{mustBeA(func, 'function_handle')}
                options.mode = ''
                options.args = {}
            end

            if strcmpi(obj.BcType, 'periodic')
                G = obj.createZeros(nCols=1);
                return;
            end

            nd = obj.Space.NDims;

            if isa(field, 'function_handle')
                G = obj.assemblePointwiseFluxVector(field, func, options.mode, options.args{:});
            else
                field = reshape(field, nd, nd, []);
                G = obj.assembleConstantFluxVector(field, func, options.mode, options.args{:});
            end
        end
    end

    methods (Access = protected)
        function P = computeInternalFluxPenalty(obj, side, ratio)
            % COMPUTEINTERNALFLUXPENALTY Compute penalty parameters for
            % internal fluxes.
            %
            %   P = computeInternalFluxPenalty(obj, side, ratio) computes
            %   penalty coefficients for internal flux terms based on the
            %   flow direction @a side and face scaling factors @a ratio.
            %   Returns penalty structure @a P for discontinuous Galerkin
            %   stabilization.

            PT = obj.PenaltyType;

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
            %   terms based on the flow direction @a side and face scaling
            %   factors @a ratio. Returns penalty structure @a P for
            %   boundary stabilization in discontinuous Galerkin methods.

            T = obj.PenaltyType;

            switch T{1, 1}
                case 'boundary'
                    switch T{1, 2}
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
                    switch T{1, 2}
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
                    switch T{1, 2}
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
                    switch T{1, 2}
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
            %   P = computePenalty(obj, side, eta) transforms base penalty
            %   parameters @a eta and flow direction @a side into the final
            %   penalty coefficients @a P used in the assembly process. The
            %   penalty formulation depends on the configured penalty type.

            T = obj.PenaltyType;
            P = cell(2, 2);

            switch T{1, 2}
                case 'right'
                    beta = 1 / 2;
                case 'left'
                    beta = -1 / 2;
                otherwise
                    beta = 0;
            end

            P{1, 2} = side .* eta{1, 2} .* beta;

            switch T{1, 1}
                case {'boundary', 'all', 'boundaryw', 'allw'}
                    alpha = 1;
                otherwise
                    alpha = 0;
            end

            P{1, 1} = side .* eta{1, 1} .* alpha;

            switch T{2, 1}
                case 'right'
                    beta = 1 / 2;
                case 'left'
                    beta = -1 / 2;
                otherwise
                    beta = 0;
            end

            P{2, 1} = side .* eta{2, 1} .* beta;
        end

        function G = assemblePointwiseVolumeMatrix(obj, field, varargin)
            % The divergence block is unweighted; the adjoint block carries D.
            nd = obj.Space.NDims;
            ne = obj.Space.NMeshElements;
            quadrature = obj.Space.Element.Volume;
            values = quadrature.Values;
            partials = conj(quadrature.Partials{1});
            weights = quadrature.Weights;
            tensor = obj.samplePointwiseTensor(field, quadrature.Nodes, 1:ne, varargin{:});
            inverse = obj.Space.Mesh.computeElementInverseJacobians();
            if size(inverse, 3) == 1, inverse = repmat(inverse, 1, 1, ne); end
            G = obj.assembleConstantVolumeMatrix(eye(nd));
            for d = 1:nd
                G{2, 1}{d} = obj.createZeros();
                for element = 1:ne
                    local = zeros(obj.Space.NLocalDofs);
                    for k = 1:nd
                        coefficient = reshape(tensor(d, k, :, element), 1, []);
                        for alpha = 1:nd
                            local = local + inverse(alpha, k, element) ...
                                * (partials(:, :, alpha).*(weights.*coefficient))*values.';
                        end
                    end
                    local = obj.Space.Element.Approximator.project(local);
                    G{2, 1}{d} = G{2, 1}{d}+obj.createMatrix(local, TEI=element, SEI=element);
                end
            end
        end

        function G = assemblePointwiseFluxMatrix(obj, field, varargin)
            % ASSEMBLEPOINTWISEFLUXMATRIX Assemble flux contribution to
            % adjoint matrix.

            nd = obj.Space.NDims;
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);

            EII = obj.Space.Mesh.getInternalElements();
            for LFI = 1:obj.Space.Element.NFluxes
                % Empty element selections otherwise mean all elements in
                % mesh/assembly helpers; an empty interior must add nothing.
                if isempty(EII), break; end
                GF = obj.assemblePointwiseInternalFluxMatrix(field, EII, LFI, varargin{:});
                for d = 1:nd
                    G{1, 2}{d} = G{1, 2}{d} + GF{1, 2}{d};
                    G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                    G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                end
            end

            EIB = obj.Space.Mesh.getBoundaryElements();
            for LFI = 1:obj.Space.Element.NFluxes
                EI = obj.Space.Mesh.findInternalElements(EIB, LFI);
                if ~isempty(EI)
                    GF = obj.assemblePointwiseInternalFluxMatrix(field, EI, LFI, varargin{:});
                    for d = 1:nd
                        G{1, 2}{d} = G{1, 2}{d} + GF{1, 2}{d};
                        G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                        G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                    end
                end

                EI = obj.Space.Mesh.findBoundaryElements(EIB, LFI);
                if ~isempty(EI)
                    switch obj.BcType
                        case 'periodic'
                            GF = obj.assemblePointwiseInternalFluxMatrix(field, EI, LFI, varargin{:});
                        case 'dirichlet'
                            GF = obj.assemblePointwiseFaceMatrix(field, EI, LFI, true, varargin{:});
                    end
                    for d = 1:nd
                        G{1, 2}{d} = G{1, 2}{d} + GF{1, 2}{d};
                        G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                        G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                    end
                end
            end
        end

        function G = assemblePointwiseInternalFluxMatrix(obj, field, EI, LFI, varargin)
            G = obj.assemblePointwiseFaceMatrix(field, EI, LFI, false, varargin{:});
        end

        function G = assemblePointwiseFaceMatrix(obj, field, EI, LFI, boundary, varargin)
            nd = obj.Space.NDims;
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);
            quadrature = obj.Space.Element.Flux(LFI);
            values = quadrature.Values;
            tensor = obj.samplePointwiseTensor(field, quadrature.Nodes, EI, varargin{:});
            normal = obj.Space.Mesh.computeOutwardNormals(EI, LFI);
            ratio = obj.Space.Mesh.computeFaceScalingFactors(EI, LFI);
            if size(normal, 2) == 1, normal = repmat(normal, 1, numel(EI)); end
            if isscalar(ratio), ratio = repmat(ratio, numel(EI), 1); end
            side = sign(sum(obj.Space.Mesh.computeDirection().*normal, 1));
            if boundary
                penalty = obj.computeDirichletBoundaryFluxPenalty(side, ratio(:).');
            else
                penalty = obj.computeInternalFluxPenalty(side, ratio(:).');
                [neighbors, neighborFaces] = obj.Space.Mesh.findNeighborElements(EI, LFI, obj.BcType);
            end
            for element = 1:numel(EI)
                for d = 1:nd
                    coefficient12 = ratio(element)*normal(d, element);
                    coefficient21 = zeros(1, quadrature.NPoints);
                    for k = 1:nd
                        coefficient21 = coefficient21+ratio(element)*normal(k, element) ...
                            *reshape(tensor(d, k, :, element), 1, []);
                    end
                    local12 = obj.Space.Element.Approximator.project( ...
                        (conj(values).*(quadrature.Weights*coefficient12))*values.');
                    local21 = obj.Space.Element.Approximator.project( ...
                        (conj(values).*(quadrature.Weights.*coefficient21))*values.');
                    if boundary
                        G{1, 2}{d} = G{1, 2}{d}+obj.createMatrix(local12, TEI=EI(element), SEI=EI(element));
                        G{1, 1}{d} = G{1, 1}{d}-obj.createMatrix( ...
                            local21*penalty{1, 1}(element), TEI=EI(element), SEI=EI(element));
                    else
                        G{1, 2}{d} = G{1, 2}{d}+obj.createMatrix( ...
                            local12*(0.5-penalty{1, 2}(element)), TEI=EI(element), SEI=EI(element));
                        G{2, 1}{d} = G{2, 1}{d}+obj.createMatrix( ...
                            local21*(0.5-penalty{2, 1}(element)), TEI=EI(element), SEI=EI(element));
                        G{1, 1}{d} = G{1, 1}{d}-obj.createMatrix( ...
                            local12*penalty{1, 1}(element), TEI=EI(element), SEI=EI(element));
                        neighborFace = neighborFaces(min(element, numel(neighborFaces)));
                        neighborValues = obj.Space.Element.Flux(neighborFace).Values;
                        external12 = obj.Space.Element.Approximator.project( ...
                            (conj(values).*(quadrature.Weights*coefficient12))*neighborValues.');
                        external21 = obj.Space.Element.Approximator.project( ...
                            (conj(values).*(quadrature.Weights.*coefficient21))*neighborValues.');
                        G{1, 2}{d} = G{1, 2}{d}+obj.createMatrix( ...
                            external12*(0.5+penalty{1, 2}(element)), TEI=EI(element), SEI=neighbors(element));
                        G{2, 1}{d} = G{2, 1}{d}+obj.createMatrix( ...
                            external21*(0.5+penalty{2, 1}(element)), TEI=EI(element), SEI=neighbors(element));
                        G{1, 1}{d} = G{1, 1}{d}+obj.createMatrix( ...
                            external12*penalty{1, 1}(element), TEI=EI(element), SEI=neighbors(element));
                    end
                end
            end
        end

        function tensor = samplePointwiseTensor(obj, field, nodes, EI, varargin)
            nd = obj.Space.NDims;
            points = obj.Space.Mesh.collocate(nodes, EI);
            sampled = field(points, varargin{:});
            count = size(points, 2);
            if nd > 1 && isequal(size(sampled), [nd, nd])
                tensor = repmat(sampled, 1, 1, size(nodes, 2), numel(EI));
                return;
            end
            if isscalar(sampled)
                sampled = repmat(sampled, 1, count);
            elseif size(sampled, 2) ~= count && size(sampled, 1) == count
                sampled = sampled.';
            end
            if isequal(size(sampled), [1, count])
                tensor = reshape(eye(nd), nd*nd, 1).*sampled;
            elseif isequal(size(sampled), [nd*nd, count])
                tensor = sampled;
            else
                error('mdemat:fem:InvalidTensorField', ...
                    'Tensor samples require scalar or d^2 components by physical points.');
            end
            if any(~isfinite(tensor(:)))
                error('mdemat:fem:InvalidTensorField', 'Tensor samples must be finite.');
            end
            tensor = reshape(tensor, nd, nd, size(nodes, 2), numel(EI));
        end

        function G = assemblePointwiseFluxVector(obj, field, func, mode, varargin)
            nd = obj.Space.NDims;
            nf = obj.Space.Element.NFluxes;
            if strcmp(mode, 'face_full'), count = nf; else, count = nd; end
            G = repmat({repmat({obj.createZeros(nCols=1)}, 1, count)}, 2, 2);
            for face = 1:nf
                selected = obj.Space.Mesh.getBoundaryElements(face);
                if isempty(selected), continue; end
                quadrature = obj.Space.Element.Flux(face);
                values = quadrature.Values;
                tensor = obj.samplePointwiseTensor(field, quadrature.Nodes, selected, varargin{:});
                normal = obj.Space.Mesh.computeOutwardNormals(selected, face);
                ratio = obj.Space.Mesh.computeFaceScalingFactors(selected, face);
                if size(normal, 2) == 1, normal = repmat(normal, 1, numel(selected)); end
                if isscalar(ratio), ratio = repmat(ratio, numel(selected), 1); end
                side = sign(sum(obj.Space.Mesh.computeDirection().*normal, 1));
                penalty = obj.computeDirichletBoundaryFluxPenalty(side, ratio(:).');
                faceParts = repmat({repmat({obj.createZeros(nCols=1)}, 1, nd)}, 2, 2);
                for element = 1:numel(selected)
                    if ismember(mode, {'face', 'face_full'})
                        args = [{face}, varargin];
                    else
                        args = varargin;
                    end
                    trace = obj.Space.feval(func, quadrature.Nodes, EI=selected(element), args=args);
                    for d = 1:nd
                        coefficient = zeros(1, quadrature.NPoints);
                        for k = 1:nd
                            coefficient = coefficient+ratio(element)*normal(k, element) ...
                                *reshape(tensor(d, k, :, element), 1, []);
                        end
                        local = obj.Space.Element.Approximator.project( ...
                            (conj(values).*(quadrature.Weights.*coefficient))*trace);
                        vector = obj.createVector(local, TEI=selected(element));
                        faceParts{2, 1}{d} = faceParts{2, 1}{d}+vector;
                        faceParts{1, 1}{d} = faceParts{1, 1}{d}+penalty{1, 1}(element)*vector;
                    end
                end
                if strcmp(mode, 'face_full')
                    G{1, 1}{face} = faceParts{1, 1};
                    G{2, 1}{face} = faceParts{2, 1};
                else
                    for d = 1:nd
                        G{1, 1}{d} = G{1, 1}{d}+faceParts{1, 1}{d};
                        G{2, 1}{d} = G{2, 1}{d}+faceParts{2, 1}{d};
                    end
                end
            end
        end

        function G = assembleConstantVolumeMatrix(obj, field)
            % ASSEMBLEVOLUMEMATRIX Assemble volume contribution to elliptic
            % matrix.

            %< Sizes
            nd = obj.Space.NDims;

            %< Inverse of Jacobian
            invJ = obj.Space.Mesh.computeElementInverseJacobians(); % (nd, nd, ne)

            %< Integration data
            W = diag(obj.Space.Element.Volume.Weights); % (nq, nq)
            V = conj(obj.Space.Element.Volume.Partials{1}); % (nl, nq, nd)
            U = obj.Space.Element.Volume.Values; % (nl, nq)

            %< Reference matrices
            R = pagemtimes(pagemtimes(V, W), U.'); % (nl, nl, nd)
            R = obj.Space.Element.Approximator.project(R); % (nl, nl, nd)

            %< Initialization
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);

            %< Physical matrices for (1, 2) component
            L = reshape(R, [], nd) * reshape(invJ, nd, []); % (nl^2, nd, ne)
            L = reshape(L, size(V, 1), size(U, 1), nd, []); % (nl, nl, nd, ne)

            for d = 1:nd
                Ld = squeeze(L(:, :, d, :));
                G{1, 2}{d} = obj.createMatrix(Ld);
            end

            %< Physical matrices for (2, 1) component
            field = permute(field, [2, 1, 3]);
            C = reshape(pagemtimes(invJ, field), nd, nd, []); % (nd, nd, ne)
            L = reshape(R, [], nd) * reshape(C, nd, []); % (nl^2, nd, ne)
            L = reshape(L, size(V, 1), size(U, 1), nd, []); % (nl, nl, nd, ne)

            for d = 1:nd
                Ld = squeeze(L(:, :, d, :));
                G{2, 1}{d} = obj.createMatrix(Ld);
            end
        end

        function G = assembleConstantFluxMatrix(obj, field)
            % ASSEMBLEFLUXMATRIX Assemble flux contribution to elliptic
            % matrix.

            % Initialization
            nd = obj.Space.NDims;
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);

            %< Add contribution from internal elements
            EII = obj.Space.Mesh.getInternalElements(); % (:, 1)
            for LFI = 1:obj.Space.Element.NFluxes
                if isempty(EII), break; end
                GF = obj.assembleConstantInternalFluxMatrix(field, EII, LFI);
                for d = 1:nd
                    G{1, 2}{d} = G{1, 2}{d} + GF{1, 2}{d};
                    G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                    G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                end
            end

            %< Add contribution from boundary elements
            EIB = obj.Space.Mesh.getBoundaryElements(); % (:, 1)
            for LFI = 1:obj.Space.Element.NFluxes
                %< Find elements of which LFI-th face is not on boundary
                EI = obj.Space.Mesh.findInternalElements(EIB, LFI);
                if ~isempty(EI)
                    GF = obj.assembleConstantInternalFluxMatrix(field, EI, LFI);
                    for d = 1:nd
                        G{1, 2}{d} = G{1, 2}{d} + GF{1, 2}{d};
                        G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                        G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                    end
                end

                %< Find elements of which LFI-th face is on boundary
                EI = obj.Space.Mesh.findBoundaryElements(EIB, LFI);
                if ~isempty(EI)
                    switch obj.BcType
                        case 'periodic'
                            GF = obj.assembleConstantInternalFluxMatrix(field, EI, LFI);
                        case 'dirichlet'
                            GF = obj.assembleConstantBoundaryFluxMatrix(field, EI, LFI);
                    end
                    for d = 1:nd
                        G{1, 2}{d} = G{1, 2}{d} + GF{1, 2}{d};
                        G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                        G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                    end
                end
            end
        end

        function G = assembleConstantInternalFluxMatrix(obj, field, EI, LFI)
            % ASSEMBLEINTERNALFLUXMATRIX Assemble internal flux
            % contributions.

            %< Sizes
            nd = obj.Space.NDims;
            nl = obj.Space.NLocalDofs;

            %< Outward normal vectors
            normal = obj.Space.Mesh.computeOutwardNormals(EI, LFI); % (nd, :)

            %< Face scaling factors
            ratio = obj.Space.Mesh.computeFaceScalingFactors(EI, LFI); % (:, 1)
            ratio = reshape(ratio, 1, []);

            %< Determine side
            direction = obj.Space.Mesh.computeDirection(); % (nd, 1)
            side = sign(sum(direction.*normal, 1)); % (1, :)

            %< Penalty parameters
            P = obj.computeInternalFluxPenalty(side, ratio);

            %< Face coefficients
            field = reshape(field, nd, nd, []);
            FC1 = ratio .* normal; % (nd, :)
            normal = reshape(normal, nd, 1, []); % (nd, 1, :)
            FC2 = ratio .* reshape(pagemtimes(field, normal), nd, []);

            %< Test basis values and weights
            V = conj(obj.Space.Element.Flux(LFI).Values); % (nl, nq)
            W = diag(obj.Space.Element.Flux(LFI).Weights); % (nq, nq)

            %< Initialization
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);

            %< Interior contribution
            U = obj.Space.Element.Flux(LFI).Values; % (nl, nq)
            L = V * W * U.'; % nl x nl
            L = obj.Space.Element.Approximator.project(L); % (nl, nl)
            L = reshape(L, nl, nl, 1);
            for d = 1:nd
                %< (1, 2) component
                C = reshape(FC1(d, :), 1, 1, []);
                A = obj.createMatrix(L.*C, TEI=EI, SEI=EI);
                G{1, 2}{d} = G{1, 2}{d} + A / 2;
                C = reshape(FC1(d, :).*P{1, 2}, 1, 1, []);
                A = obj.createMatrix(L.*C, TEI=EI, SEI=EI);
                G{1, 2}{d} = G{1, 2}{d} - A;
                %< (2, 1) component
                C = reshape(FC2(d, :), 1, 1, []);
                A = obj.createMatrix(L.*C, TEI=EI, SEI=EI);
                G{2, 1}{d} = G{2, 1}{d} + A / 2;
                C = reshape(FC2(d, :).*P{2, 1}, 1, 1, []);
                A = obj.createMatrix(L.*C, TEI=EI, SEI=EI);
                G{2, 1}{d} = G{2, 1}{d} - A;
                %< (1, 1) component
                C = reshape(FC1(d, :).*P{1, 1}, 1, 1, []);
                A = obj.createMatrix(L.*C, TEI=EI, SEI=EI);
                G{1, 1}{d} = G{1, 1}{d} - A;
            end

            %< Exterior contribution
            [NEI, NLFI] = obj.Space.Mesh.findNeighborElements(EI, LFI, obj.BcType); % (ne, 1)
            U = arrayfun(@(j) obj.Space.Element.Flux(j).Values, NLFI, 'Un', 0);
            U = cat(3, U{:}); % (nl, nq, ne)
            L = pagemtimes(V*W, U.'); % (nl, nl, ne)
            L = obj.Space.Element.Approximator.project(L); % (nl, nl, ne)
            for d = 1:nd
                %< (1, 2) component
                C = reshape(FC1(d, :), 1, 1, []);
                A = obj.createMatrix(L.*C, TEI=EI, SEI=NEI);
                G{1, 2}{d} = G{1, 2}{d} + A / 2;
                C = reshape(FC1(d, :).*P{1, 2}, 1, 1, []);
                A = obj.createMatrix(L.*C, TEI=EI, SEI=NEI);
                G{1, 2}{d} = G{1, 2}{d} + A;
                %< (2, 1) component
                C = reshape(FC2(d, :), 1, 1, []);
                A = obj.createMatrix(L.*C, TEI=EI, SEI=NEI);
                G{2, 1}{d} = G{2, 1}{d} + A / 2;
                C = reshape(FC2(d, :).*P{2, 1}, 1, 1, []);
                A = obj.createMatrix(L.*C, TEI=EI, SEI=NEI);
                G{2, 1}{d} = G{2, 1}{d} + A;
                %< (1, 1) component
                C = reshape(FC1(d, :).*P{1, 1}, 1, 1, []);
                A = obj.createMatrix(L.*C, TEI=EI, SEI=NEI);
                G{1, 1}{d} = G{1, 1}{d} + A;
            end
        end

        function G = assembleConstantBoundaryFluxMatrix(obj, field, EI, LFI)
            % ASSEMBLEBOUNDARYFLUXMATRIX Assemble boundary flux matrix.

            %< Sizes
            nd = obj.Space.NDims;
            nl = obj.Space.NLocalDofs;

            %< Outward normal vectors
            normal = obj.Space.Mesh.computeOutwardNormals(EI, LFI); % (nd, :)

            %< Face scaling factors
            ratio = obj.Space.Mesh.computeFaceScalingFactors(EI, LFI); % (:, 1)
            ratio = reshape(ratio, 1, []);

            %< Determine side
            direction = obj.Space.Mesh.computeDirection(); % (nd, 1)
            side = sign(sum(direction.*normal, 1)); % (1, :)

            %< Penalty parameters
            P = obj.computeDirichletBoundaryFluxPenalty(side, ratio);

            %< Face coefficients
            field = reshape(field, nd, nd, []);
            FC1 = ratio .* normal; % (nd, :)
            normal = reshape(normal, nd, 1, []); % (nd, 1, :)
            FC2 = ratio .* reshape(pagemtimes(field, normal), nd, []);

            %< Test basis values and weights
            V = conj(obj.Space.Element.Flux(LFI).Values); % (nl, nq)
            W = diag(obj.Space.Element.Flux(LFI).Weights); % (nq, nq)

            %< Initialization
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);

            %< Interior contribution
            U = obj.Space.Element.Flux(LFI).Values; % (nl, nq)
            L = V * W * U.'; % (nl, nl)
            L = obj.Space.Element.Approximator.project(L); % (nl, nl)
            L = reshape(L, nl, nl, 1);
            for d = 1:nd
                %< (1, 2) component
                C = reshape(FC1(d, :), 1, 1, []);
                A = obj.createMatrix(L.*C, TEI=EI, SEI=EI); % average
                G{1, 2}{d} = G{1, 2}{d} + A;
                %< (1, 1) component
                C = reshape(FC2(d, :).*P{1, 1}, 1, 1, []);
                A = obj.createMatrix(L.*C, TEI=EI, SEI=EI); % jump
                G{1, 1}{d} = G{1, 1}{d} - A;
            end
        end

        function G = assembleConstantFluxVector(obj, field, func, mode, varargin)
            % ASSEMBLEFLUXVECTOR Assemble flux contribution to elliptic vector.

            % Initialization
            nd = obj.Space.NDims;
            nf = obj.Space.Element.NFluxes;
            if strcmpi(mode, 'face_full')
                G = repmat({repmat({obj.createZeros(nCols=1)}, 1, nf)}, 2, 2);    
            else
                G = repmat({repmat({obj.createZeros(nCols=1)}, 1, nd)}, 2, 2);
            end

            % Contribution from boundary elements
            EIB = obj.Space.Mesh.getBoundaryElements(); % (:, 1)
            for LFI = 1:nf
                % Find elements of which LFI-th face is on boundary
                EI = obj.Space.Mesh.findBoundaryElements(EIB, LFI);
                if ~isempty(EI)
                    GF = obj.assembleConstantBoundaryFluxVector(field, func, EI, LFI, mode, varargin{:});
                    if strcmpi(mode, 'face_full')
                        G{1, 1}{LFI} = GF{1, 1};
                        G{2, 1}{LFI} = GF{2, 1};
                    else
                        for d = 1:nd
                            G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                            G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                        end
                    end
                end
            end
        end

        function G = assembleConstantBoundaryFluxVector(obj, field, func, EI, LFI, mode, varargin)
            % ASSEMBLEBOUNDARYFLUXVECTOR Assemble boundary flux vector.

            % Sizes
            nd = obj.Space.NDims;
            nl = obj.Space.NLocalDofs;
            ne = length(EI);

            % Outward normal vectors
            normal = obj.Space.Mesh.computeOutwardNormals(EI, LFI); % (nd, :)

            % Face scaling factors
            ratio = obj.Space.Mesh.computeFaceScalingFactors(EI, LFI); % (:, 1)
            ratio = reshape(ratio, 1, []); % (1, :)

            % Determine side
            direction = obj.Space.Mesh.computeDirection(); % (nd, 1)
            side = sign(sum(direction.*normal, 1)); % (1, :)

            % Penalty parameters
            P = obj.computeDirichletBoundaryFluxPenalty(side, ratio);

            % Face coefficients
            field = reshape(field, nd, nd, []);
            normal = reshape(normal, nd, 1, []); % (nd, 1, :)
            FC2 = ratio .* reshape(pagemtimes(field, normal), nd, []);

            % Test basis values and weights
            V = conj(obj.Space.Element.Flux(LFI).Values); % (nl, nq)
            W = diag(obj.Space.Element.Flux(LFI).Weights); % (nq, nq)

            % Face integration points
            nq = obj.Space.Element.Flux(LFI).NPoints;
            xRef = obj.Space.Element.Flux(LFI).Nodes; % (nd, nq)
            switch mode
                case ''
                    U = obj.Space.feval(func, xRef, EI=EI, args = varargin); % (nq*ne, nc)
                case {'face', 'face_full'}
                    U = obj.Space.feval(func, xRef, EI=EI, args = [LFI, varargin]); % (nq*ne, nc)
            end
            U = reshape(U, nq, []); % (nq, ne*nc)

            G = repmat({repmat({obj.createZeros(nCols=1)}, 1, nd)}, 2, 2);
            for d = 1:nd
                L = V * W * U; % (nl, ne*nc)
                L = reshape(L, nl, []); % (nl, ne*nc)
                L = obj.Space.Element.Approximator.project(L); % (nl, ne*nc)
                C = reshape(FC2(d, :), 1, []); % (1, ne)
                L = reshape(L, nl, ne, []);
                G{2, 1}{d} = obj.createVector(L.*C, TEI=EI);
                C = reshape(FC2(d, :).*P{1, 1}, 1, []); % (1, ne)
                G{1, 1}{d} = obj.createVector(L.*C, TEI=EI);
            end
        end
    end
end
