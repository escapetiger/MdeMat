classdef UpwindGridAssembly < fem.assembly.EulerianGridAssembly
    % UPWINDGRIDASSEMBLY Grid-based assembly for upwind schemes.
    %
    %   UpwindGridAssembly automatically selects left or right flux based
    %   on velocity direction for upwind discretizations. This class is
    %   specifically designed for convection-dominated problems where
    %   numerical stability requires upwind flux selection to prevent
    %   spurious oscillations.
    %
    %   The assembly chooses the flux direction based on the sign of the
    %   velocity component: positive velocities use left-biased fluxes
    %   (information flows from left to right), while negative velocities
    %   use right-biased fluxes (information flows from right to left).
    %
    % Examples:
    %   % Create upwind assembly with Dirichlet boundaries
    %   assembly = UpwindGridAssembly(fe, mesh, op, 1);
    %   
    %   % Assemble upwind divergence operator
    %   coeffs = [1.0, 0.5];
    %   velocity = [2.0, -1.0]; % positive x, negative y velocity
    %   matrix = assembly.divergence(coeffs, velocity);
    %
    %   % Assemble with boundary conditions
    %   bcVector = assembly.divergenceBc(coeffs, velocity, @bcFunction);
    %
    % Notes:
    %   The upwind choice ensures numerical stability for
    %   convection-dominated flows by respecting the characteristic
    %   directions of information propagation.
    %
    % See also:
    %   fem.assembly.EulerianGridAssembly, fem.assembly.GridAssembly

    methods
        function A = divergence(obj, c, a)
            % DIVERGENCE Assemble upwind divergence operator.
            %
            %   A = divergence(obj, c, a) creates the divergence operator
            %   matrix using upwind flux selection based on velocity
            %   direction. Each spatial dimension uses left flux for
            %   positive velocity and right flux for negative velocity.
            %
            % Inputs:
            %   obj - The UpwindGridAssembly object
            %   c - Coefficient vector, one per spatial dimension
            %   a - Velocity vector, one per spatial dimension (optional)
            %
            % Outputs:
            %   A - Assembled upwind divergence sparse matrix

            if nargin < 3 || isempty(a)
                a = c;
            end

            n = obj.nGlobalDofs;
            A = sparse(n, n);
            for d = 1:obj.nDims
                if c(d) == 0, continue; end
                T = obj.partial(d, c(d), a(d));
                A = A + T;
            end
        end

        function A = partial(obj, d, c, a)
            % PARTIAL Assemble upwind partial derivative operator.
            %
            %   A = partial(obj, d, c, a) creates the partial derivative
            %   operator matrix for the specified dimension using upwind
            %   flux selection. Positive velocity uses left flux, negative
            %   velocity uses right flux for numerical stability.
            %
            % Inputs:
            %   obj - The UpwindGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient value for this dimension
            %   a - Velocity component for this dimension (optional)
            %
            % Outputs:
            %   A - Assembled upwind partial derivative sparse matrix

            if nargin < 4 || isempty(a)
                a = c;
            end

            ne = obj.nTotalElements;

            %< ------------------------------------------------------------
            %< COEFFICIENTS
            %< ------------------------------------------------------------
            if isa(obj.mesh, 'approx.mesh.UniformGrid')
                h = obj.mesh.spacings{d};
                c = ones(1, ne) .* c ./ h;
            end

            if isa(obj.mesh, 'approx.mesh.NonuniformGrid')
                h = obj.mesh.spacings;
                [h{:}] = ndgrid(h{:});
                c = c ./ h{d}(:);
            end

            c = c(:);

            %< ------------------------------------------------------------
            %< VOLUME CONTRIBUTION
            %< ------------------------------------------------------------
            T1 = obj.assembleVolumeTerms(d, c);

            %< ------------------------------------------------------------
            %< FLUX CONTRIBUTION
            %< ------------------------------------------------------------
            if isa(obj.fe, 'fem.element.L2FiniteElement')
                if a > 0
                    T2 = obj.assembleLeftFluxTerms(d, c);
                else
                    T2 = obj.assembleRightFluxTerms(d, c);
                end
            end

            %< ------------------------------------------------------------
            %< TRIPLETS TO SPARSE MATRIX
            %< ------------------------------------------------------------
            n = obj.nGlobalDofs;
            T = [T1; T2];
            A = obj.sparseFromTriplets(T, n, n);
        end

        function b = divergenceBc(obj, c, a, f, varargin)
            % DIVERGENCEBC Assemble upwind divergence boundary conditions.
            %
            %   b = divergenceBc(obj, c, a, f) assembles the boundary
            %   condition vector for the upwind divergence operator using
            %   velocity-based flux selection for each spatial dimension.
            %
            % Inputs:
            %   obj - The UpwindGridAssembly object
            %   c - Coefficient vector, one per spatial dimension
            %   a - Velocity vector, one per spatial dimension
            %   f - Boundary condition function handle
            %   varargin - Additional arguments for boundary function
            %
            % Outputs:
            %   b - Assembled boundary condition sparse vector

            n = obj.nGlobalDofs;
            b = sparse(n, 1);
            for d = 1:obj.nDims
                if c(d) == 0, continue; end
                b = b + obj.partialBc(d, c(d), a(d), f, varargin{:});
            end
        end

        function b = partialBc(obj, d, c, a, f, varargin)
            % PARTIALBC Assemble upwind partial boundary condition.
            %
            %   b = partialBc(obj, d, c, a, f) assembles the boundary
            %   condition contribution for a specific dimension using
            %   upwind flux selection based on velocity direction.
            %
            % Inputs:
            %   obj - The UpwindGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient value for this dimension
            %   a - Velocity component for this dimension
            %   f - Boundary condition function handle
            %   varargin - Additional arguments for boundary function
            %
            % Outputs:
            %   b - Boundary condition contribution sparse vector

            ne = obj.nTotalElements;

            %< ------------------------------------------------------------
            %< COEFFICIENTS
            %< ------------------------------------------------------------
            if isa(obj.mesh, 'approx.mesh.UniformGrid')
                h = obj.mesh.spacings{d};
                c = ones(1, ne) .* c ./ h;
            end

            if isa(obj.mesh, 'approx.mesh.NonuniformGrid')
                h = obj.mesh.spacings;
                [h{:}] = ndgrid(h{:});
                c = c ./ h{d}(:);
            end

            c = c(:);

            %< ------------------------------------------------------------
            %< FLUX CONTRIBUTION
            %< ------------------------------------------------------------
            if isa(obj.fe, 'fem.element.L2FiniteElement')
                if a > 0
                    T = obj.assembleLeftDirichletBoundaryConditionTerms(d, c, f, varargin{:});
                else
                    T = obj.assembleRightDirichletBoundaryConditionTerms(d, c, f, varargin{:});
                end
            end

            %< ------------------------------------------------------------
            %< TRIPLETS TO SPARSE MATRIX
            %< ------------------------------------------------------------
            n = obj.nGlobalDofs;
            b = obj.sparseFromTriplets(T, n, 1);
        end
    end

    methods (Access = protected)
        function T = assembleInnerLeftFluxDirichletBoundaryTerms(obj, d, c)
            % ASSEMBLEINNERLEFTFLUXDIRICHLETBOUNDARYTERMS Assemble inner
            % left flux Dirichlet boundary terms.
            %
            %   T = assembleInnerLeftFluxDirichletBoundaryTerms(obj, d, c)
            %   computes the inner contribution to a left flux along the
            %   specified dimension on boundaries with Dirichlet boundary
            %   conditions for upwind schemes.
            %
            % Inputs:
            %   obj - The UpwindGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix [row, col, value] for boundary contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == 1 | m(:, d) == obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 1);

            T = zeros(np^2*ne, 3);

            F = G.fluxData(1, 2*d).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l), v(:));
        end

        function T = assembleOuterLeftFluxDirichletBoundaryTerms(obj, d, c)
            % ASSEMBLEOUTERLEFTFLUXDIRICHLETBOUNDARYTERMS Assemble outer
            % left flux Dirichlet boundary terms.
            %
            %   T = assembleOuterLeftFluxDirichletBoundaryTerms(obj, d, c)
            %   computes the outer contribution to a left flux along the
            %   specified dimension on boundaries with Dirichlet boundary
            %   conditions for upwind schemes.
            %
            % Inputs:
            %   obj - The UpwindGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for boundary contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.mesh.indexer.multiToLinear(m, 1);
            m(:, d) = m(:, d) - 1;
            l2 = obj.mesh.indexer.multiToLinear(m, 1);

            T = zeros(np^2*ne, 3);

            F = G.fluxData(2, 2*d-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l1), v(:));
        end

        function T = assembleInnerRightFluxDirichletBoundaryTerms(obj, d, c)
            % ASSEMBLEINNERRIGHTFLUXDIRICHLETBOUNDARYTERMS Assemble inner
            % right flux Dirichlet boundary terms.
            %
            %   T = assembleInnerRightFluxDirichletBoundaryTerms(obj, d, c)
            %   computes the inner contribution to a right flux along the
            %   specified dimension on boundaries with Dirichlet boundary
            %   conditions for upwind schemes.
            %
            % Inputs:
            %   obj - The UpwindGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for boundary contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == 1 | m(:, d) == obj.mesh.resolution(d);
            ne = numel(find(e));
            m = m(e, :);
            l = obj.mesh.indexer.multiToLinear(m, 1);

            T = zeros(np^2*ne, 3);

            F = G.fluxData(1, 2*d-1).matrix;
            [i, j, v] = find(-P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l), v(:));
        end

        function T = assembleOuterRightFluxDirichletBoundaryTerms(obj, d, c)
            % ASSEMBLEOUTERRIGHTFLUXDIRICHLETBOUNDARYTERMS Assemble outer
            % right flux Dirichlet boundary terms.
            %
            %   T = assembleOuterRightFluxDirichletBoundaryTerms(obj, d, c)
            %   computes the outer contribution to a right flux along the
            %   specified dimension on boundaries with Dirichlet boundary
            %   conditions for upwind schemes.
            %
            % Inputs:
            %   obj - The UpwindGridAssembly object
            %   d - Spatial dimension index (1, 2, or 3)
            %   c - Coefficient vector, one per element
            %
            % Outputs:
            %   T - Triplet matrix for boundary contribution

            np = obj.nLocalDofs;
            P = obj.fe.projector;
            G = obj.op.gradient;

            m = obj.mesh.allElementMultiIndices;
            e = m(:, d) == 1;
            ne = numel(find(e));
            m = m(e, :);
            l1 = obj.mesh.indexer.multiToLinear(m, 1);
            m(:, d) = m(:, d) + 1;
            l2 = obj.mesh.indexer.multiToLinear(m, 1);

            T = zeros(np^2*ne, 3);

            F = G.fluxData(2, 2*d).matrix;
            [i, j, v] = find(P.project(F));
            nv = numel(v);
            k = 1:(nv * ne);
            T(k, 1) = reshape(bsxfun(@plus, (l1(:).' - 1)*np, i(:)), [], 1);
            T(k, 2) = reshape(bsxfun(@plus, (l2(:).' - 1)*np, j(:)), [], 1);
            T(k, 3) = kron(c(l1), v(:));
        end
    end
end