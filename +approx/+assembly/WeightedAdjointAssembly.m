classdef WeightedAdjointAssembly < approx.assembly.AdjointAssembly
    % FINITEELEMENTWEIGHTEDADJOINTASSEMBLY Assembly for weighted finite
    % element adjoint operator pairs.
    %
    %   WeightedAdjointAssembly implements assembly operations
    %   for weighted finite element adjoint operator pairs. These operators
    %   generalize the standard adjoint assembly to support spatially
    %   varying weight functions. The forward operator is defined as
    %   \f$ Au = ca(x)\partial_x(u/a(x)) \f$ and the adjoint operator as
    %   \f$ A^*u = -c\partial_x(a(x)u)/a(x) \f$, where @a a(x) is a
    %   spatially varying weight function and @a c is a constant.
    %
    %   The class handles both periodic and Dirichlet boundary conditions
    %   with various penalty formulations inherited from the parent class.
    %   The weighted formulation preserves adjoint consistency while
    %   allowing for spatially varying coefficients that arise in problems
    %   with non-uniform material properties or coordinate transformations.
    %
    %   The assembly process combines volume integrals from the divergence
    %   theorem with flux integrals across element interfaces, incorporating
    %   the weight function @a a(x) and its derivatives appropriately.
    %
    % See also:
    %   approx.assembly.AdjointAssembly,
    %   approx.assembly.FiniteElementAssembly

    methods
        function obj = WeightedAdjointAssembly(space, bcType, penaltyType)
            % FINITEELEMENTWEIGHTEDADJOINTASSEMBLY Constructor for weighted
            % adjoint assembly.
            %
            %   Creates a weighted adjoint assembly object with specified
            %   boundary conditions and penalty formulation. Inherits all
            %   configuration from the parent AdjointAssembly
            %   class.

            arguments
                space approx.space.FiniteElementSpace
                bcType{mustBeMember(bcType, {'periodic', 'dirichlet'})}
                penaltyType{mustBePenaltyType(penaltyType)}
            end

            obj@approx.assembly.AdjointAssembly(space, bcType, penaltyType);
        end

        function G = assembleMatrix(obj, field, weight, options)
            % ASSEMBLEMATRIX Assemble the weighted adjoint operator matrix.
            %
            %   G = assembleMatrix(obj, field, weight) constructs the global
            %   weighted adjoint operator matrix by combining volume and
            %   flux contributions with the spatially varying weight
            %   function @a weight(x). The @a field parameter specifies the
            %   constant @a c in the operator definition.
            %
            %   G = assembleMatrix(obj, field, weight, args=args) constructs
            %   the global weighted adjoint operator with additional
            %   arguments @a args passed to the weight function evaluation.
            %
            %   G = assembleMatrix(obj, field, weight, dWeight=dWeight)
            %   uses the analytical derivative @a dWeight of the weight
            %   function for accurate volume integral assembly via the
            %   product rule.

            arguments
                obj approx.assembly.WeightedAdjointAssembly
                field{mustBeA(field, 'numeric')}
                weight{mustBeA(weight, 'function_handle')}
                options.dWeight{mustBeA(options.dWeight, 'function_handle')} = []
                options.args = {}
            end

            nd = obj.Space.NDims;

            GV = obj.assembleConstantVolumeMatrix(field, weight, options.dWeight, options.args{:});
            GF = obj.assembleConstantFluxMatrix(field, weight, options.args{:});

            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);
            for d = 1:nd
                G{1, 2}{d} = GF{1, 2}{d} - GV{1, 2}{d};
                G{2, 1}{d} = GF{2, 1}{d} - GV{2, 1}{d};
                G{1, 1}{d} = GF{1, 1}{d};
            end
        end

        function G = assembleVector(obj, field, weight, func, options)
            % ASSEMBLEVECTOR Assemble the weighted adjoint operator vector.
            %
            %   G = assembleVector(obj, field, weight, func) constructs the
            %   right-hand side vector for weighted adjoint problems with
            %   Dirichlet boundary conditions using boundary function @a
            %   func and weight function @a weight(x).
            %
            %   G = assembleVector(obj, field, weight, func, options) allows
            %   additional arguments to be passed to function evaluation
            %   through @a options.args.

            arguments
                obj approx.assembly.WeightedAdjointAssembly
                field{mustBeA(field, 'numeric')}
                weight{mustBeA(weight, 'function_handle')}
                func{mustBeA(func, 'function_handle')}
                options.args = {}
            end

            if strcmpi(obj.BcType, 'periodic')
                G = obj.createZeros(nCols=1);
                return;
            end

            nd = obj.Space.NDims;

            GF = obj.assembleConstantFluxVector(field, weight, func, options.args{:});

            G = repmat({repmat({obj.createZeros(nCols=1)}, 1, nd)}, 2, 2);
            for d = 1:nd
                G{1, 1}{d} = GF{1, 1}{d};
                G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
            end
        end
    end

    methods (Access = protected)
        function G = assembleConstantVolumeMatrix(obj, field, weight, dWeight, varargin)
            % ASSEMBLECONSTANTVOLUMEMATRIX Assemble volume matrix for
            % weighted adjoint operator with constant field.
            %
            %   G = assembleConstantVolumeMatrix(obj, field, weight, dWeight)
            %   assembles the volume contribution to the weighted adjoint
            %   operator matrix with constant @a field and spatially varying
            %   @a weight function. Uses integration by parts with the product
            %   rule to correctly handle ∂ₓ(au) and ∂ₓ(u/a) terms.
            %
            %   For (1,2): -c (∫ (∂ₓv)·u dx - ∫ v·u·(∂ₓa)/a dx)
            %   For (2,1): +c (∫ (∂ₓv)·u dx + ∫ v·u·(∂ₓa)/a dx)

            %< Sizes
            nd = obj.Space.NDims;
            ne = obj.Space.NMeshElements;
            nl = obj.Space.NLocalDofs;
            nq = obj.Space.Element.Volume.NPoints;

            %< Inverse of Jacobian and determinant
            invJ = obj.Space.Mesh.computeElementInverseJacobians(); % (nd, nd, ne)

            %< Evaluate weight function a(x) at quadrature points
            xRef = obj.Space.Element.Volume.Nodes;
            a = obj.Space.feval(weight, xRef, args=varargin);
            a = reshape(a, nq, ne); % (nq, ne)

            %< Integration data
            W0 = obj.Space.Element.Volume.Weights; % (nq,)
            V = obj.Space.Element.Volume.Partials{1}; % (nl, nq, nd)
            U = obj.Space.Element.Volume.Values; % (nl, nq)

            %< Initialization
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);

            %< Assemble (1, 2) component: A*u = -c∂_x(au)/a
            % Weak form: -c(∫ (∂ₓv)·u dx - ∫ v·u·(∂ₓa)/a dx)

            % Standard term: ∫ (∂ₓv)·u dx
            % Vectorized: R1(:,:,d) = Vd * diag(W0) * U'
            R1 = pagemtimes(pagemtimes(V, diag(W0)), U'); % (nl, nl, nd)
            R1 = obj.Space.Element.Approximator.project(R1); % (nl, nl, nd)

            % Correction term with ∂ₓa if provided
            if ~isempty(dWeight)
                dAdx = obj.Space.feval(dWeight, xRef, args=varargin);
                dAdx = reshape(dAdx, nq, ne); % (nq, ne)

                % Vectorized: ∫ v·u·(∂ₓa)/a dx
                % weights: (nq, ne) - include detJ for volume integral
                weights = W0(:) .* (dAdx ./ a);
                % U: (nl, nq), weights: (nq, ne) -> broadcast and compute
                % R2(:,:,ei) = U * diag(weights(:,ei)) * U'
                % = sum_q U(:,q) * weights(q,ei) * U(:,q)'
                Uw = reshape(U, nl, nq, 1) .* reshape(weights, 1, nq, ne); % (nl, nq, ne)
                R2 = pagemtimes(Uw, U'); % (nl, nl, ne)
                R2 = obj.Space.Element.Approximator.project(R2);
            else
                R2 = zeros(nl, nl, ne);
            end

            % Apply inverse Jacobian transformation to both R1 and R2
            % Following the pattern from AdjointAssembly

            % Transform R1: L = reshape(R, [], nd) * reshape(invJ, nd, [])
            L1 = reshape(R1, [], nd) * reshape(invJ, nd, []); % (nl^2, nd*ne)
            L1 = reshape(L1, nl, nl, nd, ne); % (nl, nl, nd, ne)

            % Transform R2 if dWeight provided
            if ~isempty(dWeight)
                L = L1 - reshape(R2, nl, nl, 1, ne);
            else
                L = L1;
            end

            for d = 1:nd
                Ld = squeeze(L(:, :, d, :));
                G{1, 2}{d} = obj.createMatrix(Ld) * (-field);
            end

            %< Assemble (2, 1) component: Au = ca∂_x(u/a)
            % Weak form: +c (∫ (∂ₓv)·u dx + ∫ v·u·(∂ₓa)/a dx)
            if ~isempty(dWeight)
                L = L1 + reshape(R2, nl, nl, 1, ne);
            else
                L = L1;
            end

            % Reuse L from (1,2) component (already computed above)
            for d = 1:nd
                Ld = squeeze(L(:, :, d, :));
                G{2, 1}{d} = obj.createMatrix(Ld) * field;
            end
        end

        function G = assembleConstantFluxMatrix(obj, field, weight, varargin)
            % ASSEMBLECONSTANTFLUXMATRIX Assemble flux contribution to
            % weighted adjoint matrix with constant field.
            %
            %   G = assembleConstantFluxMatrix(obj, field, weight)
            %   assembles flux contributions for the weighted adjoint
            %   operator with constant @a field, incorporating the weight
            %   function @a weight(x) in
            %   the flux terms.

            nd = obj.Space.NDims;
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);

            EII = obj.Space.Mesh.getInternalElements();
            for LFI = 1:obj.Space.Element.NFluxes
                GF = obj.assembleConstantInternalFluxMatrix(weight, EII, LFI, varargin{:});
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
                    GF = obj.assembleConstantInternalFluxMatrix(weight, EI, LFI, varargin{:});
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
                            GF = obj.assembleConstantInternalFluxMatrix(weight, EI, LFI, varargin{:});
                        case 'dirichlet'
                            GF = obj.assembleConstantBoundaryFluxMatrix(weight, EI, LFI, varargin{:});
                    end
                    for d = 1:nd
                        G{1, 2}{d} = G{1, 2}{d} + GF{1, 2}{d};
                        G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                        G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                    end
                end
            end

            %< Apply constant field coefficient
            % A* has negative sign: A*u = -c∂_x(au)/a
            % A has positive sign: Au = ca∂_x(u/a)
            for d = 1:nd
                G{1, 2}{d} = G{1, 2}{d} * (-field);
                G{2, 1}{d} = G{2, 1}{d} * field;
                G{1, 1}{d} = G{1, 1}{d} * (-field);
            end
        end

        function G = assembleConstantInternalFluxMatrix(obj, weight, EI, LFI, varargin)
            % ASSEMBLECONSTANTINTERNALFLUXMATRIX Assemble internal flux
            % contributions for weighted adjoint operator.
            %
            %   (1,2) component assembles A*: A*u = -c∂_x(au)/a
            %   Flux term: (v/a)·n{{au}} (field coefficient applied at final stage)
            %
            %   (2,1) component assembles A: Au = ca∂_x(u/a)
            %   Flux term: v·a·n{{u/a}} (field coefficient applied at final stage)

            %< Sizes
            nd = obj.Space.NDims;
            nl = obj.Space.NLocalDofs;
            ne = length(EI);
            nq = obj.Space.Element.Flux(LFI).NPoints;

            %< Outward normal vectors
            normal = obj.Space.Mesh.computeOutwardNormals(EI, LFI); % (nd, ne)

            %< Face scaling factors
            ratio = obj.Space.Mesh.computeFaceScalingFactors(EI, LFI); % (ne,)
            ratio = reshape(ratio, 1, []);

            %< Determine side
            direction = obj.Space.Mesh.computeDirection(); % (nd, 1)
            side = sign(sum(direction.*normal, 1)); % (1, ne)

            %< Penalty parameters
            P = obj.computeInternalFluxPenalty(side, ratio);

            %< Evaluate weight at flux quadrature points
            xRef = obj.Space.Element.Flux(LFI).Nodes;
            a = obj.Space.feval(weight, xRef, EI=EI, args=varargin);
            a = reshape(a, nq, ne); % (nq, ne)

            %< Test basis values and weights
            V = obj.Space.Element.Flux(LFI).Values; % (nl, nq)
            W0 = obj.Space.Element.Flux(LFI).Weights; % (nq,)

            %< Initialization
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);

            %< Interior contribution
            U = obj.Space.Element.Flux(LFI).Values; % (nl, nq)

            % (1,2) A*: Test by 1/a, trial by a
            Vw = V ./ reshape(a, 1, nq, ne); % (nl, nq, ne)
            Uw = U .* reshape(a, 1, nq, ne); % (nl, nq, ne)
            L12 = pagemtimes(Vw * diag(W0), pagetranspose(Uw)); % (nl, nl, ne)
            L12 = obj.Space.Element.Approximator.project(L12); % (nl, nl, ne)

            % (2,1) A: Test by a, trial by 1/a
            Vw = V .* reshape(a, 1, nq, ne); % (nl, nq, ne)
            Uw = U ./ reshape(a, 1, nq, ne); % (nl, nq, ne)
            L21 = pagemtimes(Vw * diag(W0), pagetranspose(Uw)); % (nl, nl, ne)
            L21 = obj.Space.Element.Approximator.project(L21); % (nl, nl, ne)

            % (1,1) A* penalty: Both by 1/a
            Vw = V ./ reshape(a, 1, nq, ne); % (nl, nq, ne)
            Uw = U ./ reshape(a, 1, nq, ne); % (nl, nq, ne)
            L11 = pagemtimes(Vw * diag(W0), pagetranspose(Uw)); % (nl, nl, ne)
            L11 = obj.Space.Element.Approximator.project(L11); % (nl, nl, ne)

            for d = 1:nd
                %< (1, 2) component
                C = reshape(ratio .* normal(d, :), 1, 1, []);
                A = obj.createMatrix(L12.*C, TEI=EI, SEI=EI);
                G{1, 2}{d} = G{1, 2}{d} + A / 2;
                C = reshape(ratio .* normal(d, :) .* P{1, 2}, 1, 1, []);
                A = obj.createMatrix(L12.*C, TEI=EI, SEI=EI);
                G{1, 2}{d} = G{1, 2}{d} - A;

                %< (2, 1) component
                C = reshape(ratio .* normal(d, :), 1, 1, []);
                A = obj.createMatrix(L21.*C, TEI=EI, SEI=EI);
                G{2, 1}{d} = G{2, 1}{d} + A / 2;
                C = reshape(ratio .* normal(d, :) .* P{2, 1}, 1, 1, []);
                A = obj.createMatrix(L21.*C, TEI=EI, SEI=EI);
                G{2, 1}{d} = G{2, 1}{d} - A;

                %< (1, 1) component
                C = reshape(ratio .* normal(d, :) .* P{1, 1}, 1, 1, []);
                A = obj.createMatrix(L11.*C, TEI=EI, SEI=EI);
                G{1, 1}{d} = G{1, 1}{d} - A;
            end

            %< Exterior contribution
            [NEI, NLFI] = obj.Space.Mesh.findNeighborElements(EI, LFI, obj.BcType);
            U = arrayfun(@(j) obj.Space.Element.Flux(j).Values, NLFI, 'Un', 0);
            U = cat(3, U{:}); % (nl, nq, ne)
            
            % Evaluate weight at neighbor elements for trial functions
            xRefNeighbor = arrayfun(@(j) obj.Space.Element.Flux(j).Nodes, NLFI, 'Un', 0);
            xRefNeighbor = cat(2, xRefNeighbor{:}); % (nd, nq*ne)
            aNeighbor = obj.Space.feval(weight, xRefNeighbor, EI=NEI, args=varargin);
            aNeighbor = reshape(aNeighbor, nq, ne); % (nq, ne)

            %  (current) (neighbor)
            % (1,2) A*: Test by 1/a (current), trial by a (neighbor)
            Vw = V ./ reshape(a, 1, nq, ne); % (nl, nq, ne)
            Uw = U .* reshape(aNeighbor, 1, nq, ne); % (nl, nq, ne)
            L12 = pagemtimes(Vw * diag(W0), pagetranspose(Uw)); % (nl, nl, ne)
            L12 = obj.Space.Element.Approximator.project(L12); % (nl, nl, ne)
            % (2,1) A: Test by a (current), trial by 1/a (neighbor)
            Vw = V .* reshape(a, 1, nq, ne); % (nl, nq, ne)
            Uw = U ./ reshape(aNeighbor, 1, nq, ne); % (nl, nq, ne)
            L21 = pagemtimes(Vw * diag(W0), pagetranspose(Uw)); % (nl, nl, ne)
            L21 = obj.Space.Element.Approximator.project(L21); % (nl, nl, ne)
 
            % (1,1) A* penalty: Both by 1/a (current for test, neighbor for trial)
            Vw = V ./ reshape(a, 1, nq, ne); % (nl, nq, ne)
            Uw = U ./ reshape(aNeighbor, 1, nq, ne); % (nl, nq, ne)
            L11 = pagemtimes(Vw * diag(W0), pagetranspose(Uw)); % (nl, nl, ne)
            L11 = obj.Space.Element.Approximator.project(L11); % (nl, nl, ne)

            for d = 1:nd
                %< (1, 2) component
                C = reshape(ratio .* normal(d, :), 1, 1, []);
                A = obj.createMatrix(L12.*C, TEI=EI, SEI=NEI);
                G{1, 2}{d} = G{1, 2}{d} + A / 2;
                C = reshape(ratio .* normal(d, :) .* P{1, 2}, 1, 1, []);
                A = obj.createMatrix(L12.*C, TEI=EI, SEI=NEI);
                G{1, 2}{d} = G{1, 2}{d} + A;

                %< (2, 1) component
                C = reshape(ratio .* normal(d, :), 1, 1, []);
                A = obj.createMatrix(L21.*C, TEI=EI, SEI=NEI);
                G{2, 1}{d} = G{2, 1}{d} + A / 2;
                C = reshape(ratio .* normal(d, :) .* P{2, 1}, 1, 1, []);
                A = obj.createMatrix(L21.*C, TEI=EI, SEI=NEI);
                G{2, 1}{d} = G{2, 1}{d} + A;

                %< (1, 1) component
                C = reshape(ratio .* normal(d, :) .* P{1, 1}, 1, 1, []);
                A = obj.createMatrix(L11.*C, TEI=EI, SEI=NEI);
                G{1, 1}{d} = G{1, 1}{d} + A;
            end
        end

        function G = assembleConstantBoundaryFluxMatrix(obj, weight, EI, LFI, varargin)
            % ASSEMBLECONSTANTBOUNDARYFLUXMATRIX Assemble boundary flux
            % matrix for weighted adjoint operator.
            %
            %   (1,2) for A*: boundary flux with weight 1/a for test, a for trial
            %   Field coefficient applied at final assembly stage.

            %< Sizes
            nd = obj.Space.NDims;
            nl = obj.Space.NLocalDofs;
            ne = length(EI);
            nq = obj.Space.Element.Flux(LFI).NPoints;

            %< Outward normal vectors
            normal = obj.Space.Mesh.computeOutwardNormals(EI, LFI); % (nd, ne)

            %< Face scaling factors
            ratio = obj.Space.Mesh.computeFaceScalingFactors(EI, LFI); % (ne,)
            ratio = reshape(ratio, 1, []);

            %< Determine side
            direction = obj.Space.Mesh.computeDirection(); % (nd, 1)
            side = sign(sum(direction.*normal, 1)); % (1, ne)

            %< Penalty parameters
            P = obj.computeDirichletBoundaryFluxPenalty(side, ratio);

            %< Evaluate weight at flux quadrature points
            xRef = obj.Space.Element.Flux(LFI).Nodes;
            a = obj.Space.feval(weight, xRef, EI=EI, args=varargin);
            a = reshape(a, nq, ne); % (nq, ne)

            %< Test basis values and weights
            V = obj.Space.Element.Flux(LFI).Values; % (nl, nq)
            W0 = obj.Space.Element.Flux(LFI).Weights; % (nq,)

            %< Initialization
            G = repmat({repmat({obj.createZeros()}, 1, nd)}, 2, 2);

            %< Boundary contribution
            U = obj.Space.Element.Flux(LFI).Values; % (nl, nq)

            % (1,2) A*: Test by 1/a, trial by a
            Vw = V ./ reshape(a, 1, nq, ne); % (nl, nq, ne)
            Uw = U .* reshape(a, 1, nq, ne); % (nl, nq, ne)
            L12 = pagemtimes(Vw * diag(W0), pagetranspose(Uw)); % (nl, nl, ne)
            L12 = obj.Space.Element.Approximator.project(L12); % (nl, nl, ne)

            % (1,1) A* penalty: Both by 1/a
            Vw = V ./ reshape(a, 1, nq, ne); % (nl, nq, ne)
            Uw = U ./ reshape(a, 1, nq, ne); % (nl, nq, ne)
            L11 = pagemtimes(Vw * diag(W0), pagetranspose(Uw)); % (nl, nl, ne)
            L11 = obj.Space.Element.Approximator.project(L11); % (nl, nl, ne)

            for d = 1:nd
                %< (1, 2) component
                C = reshape(ratio .* normal(d, :), 1, 1, []);
                A = obj.createMatrix(L12.*C, TEI=EI, SEI=EI);
                G{1, 2}{d} = G{1, 2}{d} + A;

                %< (1, 1) component
                C = reshape(ratio .* normal(d, :) .* P{1, 1}, 1, 1, []);
                A = obj.createMatrix(L11.*C, TEI=EI, SEI=EI);
                G{1, 1}{d} = G{1, 1}{d} - A;
            end
        end

        function G = assembleConstantFluxVector(obj, field, weight, func, varargin)
            % ASSEMBLECONSTANTFLUXVECTOR Assemble flux vector for weighted
            % adjoint operator.

            nd = obj.Space.NDims;
            G = repmat({repmat({obj.createZeros(nCols=1)}, 1, nd)}, 2, 2);

            EIB = obj.Space.Mesh.getBoundaryElements();
            for LFI = 1:obj.Space.Element.NFluxes
                EI = obj.Space.Mesh.findBoundaryElements(EIB, LFI);
                if ~isempty(EI)
                    GF = obj.assembleConstantBoundaryFluxVector(weight, func, EI, LFI, varargin{:});
                    for d = 1:nd
                        G{1, 1}{d} = G{1, 1}{d} + GF{1, 1}{d};
                        G{2, 1}{d} = G{2, 1}{d} + GF{2, 1}{d};
                    end
                end
            end

            %< Apply constant field coefficient
            % A* has negative sign: A*u = -c∂_x(au)/a
            % A has positive sign: Au = ca∂_x(u/a)
            for d = 1:nd
                G{1, 1}{d} = G{1, 1}{d} * field;
                G{2, 1}{d} = G{2, 1}{d} * field;
            end
        end

        function G = assembleConstantBoundaryFluxVector(obj, weight, func, EI, LFI, varargin)
            % ASSEMBLECONSTANTBOUNDARYFLUXVECTOR Assemble boundary flux
            % vector for weighted adjoint operator.
            %
            %   (2,1) for A: boundary flux with weight a
            %   (1,1) for A*: boundary flux with standard weights

            nd = obj.Space.NDims;
            nl = obj.Space.NLocalDofs;
            ne = length(EI);
            nq = obj.Space.Element.Flux(LFI).NPoints;
            normal = obj.Space.Mesh.computeOutwardNormals(EI, LFI);
            ratio = obj.Space.Mesh.computeFaceScalingFactors(EI, LFI);
            ratio = reshape(ratio, 1, []);
            direction = obj.Space.Mesh.computeDirection();
            side = sign(sum(direction.*normal, 1));
            P = obj.computeDirichletBoundaryFluxPenalty(side, ratio);

            V = obj.Space.Element.Flux(LFI).Values;
            xRef = obj.Space.Element.Flux(LFI).Nodes;

            % Evaluate weight at flux quadrature points
            a = obj.Space.feval(weight, xRef, EI=EI, args=varargin);
            a = reshape(a, nq, ne);

            W0 = obj.Space.Element.Flux(LFI).Weights;

            % Evaluate boundary function
            U = obj.Space.feval(func, xRef, EI=EI, args=varargin);
            U = reshape(U, nq, ne);

            G = repmat({repmat({obj.createZeros(nCols=1)}, 1, nd)}, 2, 2);
            for d = 1:nd
                % For (2,1) component (A): integrates v·a·n(u_bc/a)
                % Weight test function by a, trial function by 1/a
                L = zeros(nl, ne);
                for ei = 1:ne
                    Vw = V .* reshape(a(:, ei), 1, nq); % (nl, nq) - weight by a
                    Uw = U(:, ei) ./ a(:, ei); % (nq, 1) - weight by 1/a
                    W = diag(W0);
                    L(:, ei) = Vw * W * Uw;
                end
                L = obj.Space.Element.Approximator.project(L);

                C = reshape(ratio .* normal(d, :), 1, []);
                G{2, 1}{d} = obj.createVector(L.*C, TEI=EI);

                % For (1,1) component (A*): integrates (v/a)·n(au_bc)
                % Weight test function by 1/a, trial function by a
                L = zeros(nl, ne);
                for ei = 1:ne
                    Vw = V ./ reshape(a(:, ei), 1, nq); % (nl, nq) - weight by 1/a
                    Uw = U(:, ei) .* a(:, ei); % (nq, 1) - weight by a
                    W = diag(W0);
                    L(:, ei) = Vw * W * Uw;
                end
                L = obj.Space.Element.Approximator.project(L);

                C = reshape(ratio .* normal(d, :) .* P{1, 1}, 1, []);
                G{1, 1}{d} = obj.createVector(L.*C, TEI=EI);
            end
        end
    end
end
